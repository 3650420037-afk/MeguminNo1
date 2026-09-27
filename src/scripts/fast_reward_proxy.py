# -*- coding: utf-8 -*-
"""快速奖励代理: 从**已有三维位姿**直接算口袋接触特征, 目标 <50 ms/分子。

为什么需要它(目标路线 ②)
------------------------
用 Vina 对接做训练内循环的奖励太慢(实测 ~16 s/分子)。要在 ReST 式
"生成 -> 按奖励筛选 -> 再微调"里当筛选器, 必须有一个**毫秒级**的代理。
本模块不做对接、不做搜索: 它只读分子**已经生成好的位姿**与口袋原子,
算几何接触特征并给出一个接触分。

它是什么 / 不是什么
-------------------
- **是**: 一个几何接触代理 —— 埋藏度、疏水接触、极性接触、穿模惩罚。
- **不是**: 亲和力预测器。它的价值必须**用与真实 Vina 分的相关性来验证**,
  本文件自带 `--validate` 正是做这件事(用已有对接结果做秩相关)。
  未经验证的代理不得用于筛选, 更不得替代 A/B 主指标。

用法
----
    # 1) 单个位姿打分(并报耗时)
    python src/scripts/fast_reward_proxy.py --pose X_out.pdbqt --target A2A

    # 2) 用已有对接结果验证相关性(必须做)
    python src/scripts/fast_reward_proxy.py --validate \
        --summary outputs/mass_dock/A2A/summary_merged.csv \
        --poses outputs/mass_dock/A2A/docking --target A2A --limit 2000
"""
import argparse
import csv
import json
import math
import os
import statistics as st
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, CONFIGS  # noqa: E402

# 接触判据(Å)。取自常见弱相互作用范围, 保守取值。
CONTACT = 4.5        # 非极性重原子接触
HBOND = 3.6          # 极性原子间
CLASH = 2.0          # 明显穿模
BURY_R = 5.0         # 判定"被埋藏"的最近蛋白原子距离
POLAR = {"N", "O"}


def _elem_from_pdbqt(tok):
    """PDBQT 的原子类型形如 C / A / OA / NA / HD, 归一到元素符号。"""
    t = tok.strip()
    if not t:
        return "C"
    if t in ("A", "C"):
        return "C"
    if t.startswith("O"):
        return "O"
    if t.startswith("N"):
        return "N"
    if t.startswith("S"):
        return "S"
    if t.startswith("H") or t == "HD":
        return "H"
    return t[0].upper()


def load_receptor(target, radius=18.0):
    """口袋中心 R 内的蛋白重原子 (坐标 + 元素)。只读一次, 供全部分子复用。"""
    reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    t = reg[target]
    pdb = os.path.join(ROOT, t["pdb"])
    c = np.asarray(t["center"], dtype=float)
    xyz, elems = [], []
    with open(pdb, encoding="utf-8", errors="ignore") as f:
        for line in f:
            # **只取 ATOM(蛋白), 不取 HETATM** —— 受体 PDB 里含共晶配体/水/离子,
            # 而 Vina 正是把分子对接到共晶配体所在的口袋 -> 若把 HETATM 也当口袋原子,
            # 每个对接位姿都会与共晶配体重叠, 实测 clash=1.0(所有原子 <2 Å),
            # 代理分完全失真。这个坑我踩过一次。
            if not line.startswith("ATOM"):
                continue
            try:
                p = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            e = (line[76:78].strip() or line[12:16].strip()[:1] or "C").upper()
            if e == "H":
                continue
            if np.linalg.norm(p - c) <= radius:
                xyz.append(p)
                elems.append(e)
    return np.asarray(xyz, dtype=np.float32), np.asarray(elems)


def load_pose(path, model=1):
    """读 PDBQT 位姿里**指定 MODEL** 的配体重原子坐标与元素。

    ⚠ 必须只取一个 MODEL: Vina 输出默认含 **9 个位姿**(MODEL 1..9)。若把全部
    ATOM 行当一个配体读, 会把 9 个互相叠合的构象混在一起 —— 实测"配体"变成
    207 个重原子、93% 原子与蛋白距离 <2 Å(物理上不可能), 相关性也必然无效。
    这个坑我踩过一次, 故在此写死并注明。
    """
    xyz, elems = [], []
    cur = 0
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("MODEL"):
                cur += 1
                if cur > model:
                    break
                continue
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if cur not in (0, model):      # 0 = 无 MODEL 标记的单构象文件
                if cur > model:
                    break
                continue
            try:
                p = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            tok = line[77:79].strip() or line[12:16].strip()[:1]
            e = _elem_from_pdbqt(tok)
            if e == "H":
                continue
            xyz.append(p)
            elems.append(e)
    if not xyz:
        return None, None
    return np.asarray(xyz, dtype=np.float32), np.asarray(elems)


def score(lig_xyz, lig_e, rec_xyz, rec_e):
    """接触分(越大越"贴得好")。全 numpy, 无搜索, 实测 <5 ms。

    组成(全部可解释):
      buried   被埋藏原子比例(最近蛋白原子 < BURY_R)
      contact  log1p(4.5 Å 内的重原子对数的归一化)
      hydro    疏水接触(C-C)占比
      polar    极性接触(N/O-N/O < 3.6 Å)计数
      clash    明显穿模原子比例(最近距离 < CLASH) —— 作为**惩罚**
    """
    d = np.linalg.norm(lig_xyz[:, None, :] - rec_xyz[None, :, :], axis=2)  # (M,N)
    dmin = d.min(axis=1)
    m = d.shape[0]
    buried = float((dmin < BURY_R).sum()) / m
    n_pair = int((d < CONTACT).sum())
    contact = math.log1p(n_pair) / math.log1p(max(1, m * 8))
    lp = np.isin(lig_e, list(POLAR))
    rp = np.isin(rec_e, list(POLAR))
    hydro_pairs = int((d < CONTACT)[~lp][:, ~rp].sum())
    hydro = hydro_pairs / max(1, n_pair)
    polar = int((d < HBOND)[lp][:, rp].sum())
    clash = float((dmin < CLASH).sum()) / m
    s = (1.0 * buried + 1.0 * contact + 0.5 * hydro + 0.10 * polar - 3.0 * clash)
    return dict(score=s, buried=buried, contact=contact, hydro=hydro,
                polar=polar, clash=clash, n_lig=m, n_pair=n_pair)


def _spearman(a, b):
    """秩相关(不依赖 scipy): 用于检验代理分与真实 Vina 分的单调一致性。"""
    def rank(x):
        idx = np.argsort(x, kind="mergesort")
        r = np.empty(len(x), dtype=float)
        r[idx] = np.arange(len(x), dtype=float)
        # 处理并列: 取平均秩
        xs = x[idx]
        i = 0
        while i < len(xs):
            j = i
            while j + 1 < len(xs) and xs[j + 1] == xs[i]:
                j += 1
            if j > i:
                r[idx[i:j + 1]] = (i + j) / 2.0
            i = j + 1
        return r
    ra, rb = rank(np.asarray(a, float)), rank(np.asarray(b, float))
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    den = math.sqrt(float((ra ** 2).sum()) * float((rb ** 2).sum()))
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(description="快速奖励代理(几何接触, 目标 <50ms/分子)")
    ap.add_argument("--pose", help="单个 PDBQT 位姿")
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--validate", action="store_true",
                    help="用已有对接结果检验代理分与 Vina 分的秩相关")
    ap.add_argument("--summary", help="对接 summary.csv(含 name/smiles/vina_score/status)")
    ap.add_argument("--poses", help="位姿目录")
    ap.add_argument("--limit", type=int, default=2000)
    args = ap.parse_args()

    rec_xyz, rec_e = load_receptor(args.target)
    print("口袋原子 %d 个 (target=%s)" % (len(rec_xyz), args.target))

    if args.validate:
        if not (args.summary and args.poses):
            raise SystemExit("--validate 需要 --summary 与 --poses")
        rows = [r for r in csv.DictReader(open(args.summary, encoding="utf-8-sig"))
                if r.get("status") == "ok"][:args.limit]
        print("读入 %d 条 ok 记录, 开始打分 ..." % len(rows))
        vs, ss, tms, miss = [], [], [], 0
        t0 = time.time()
        for r in rows:
            p = os.path.join(args.poses, r["name"] + "_out.pdbqt")
            if not os.path.exists(p):
                miss += 1
                continue
            lx, le = load_pose(p)
            if lx is None:
                miss += 1
                continue
            ta = time.perf_counter()
            out = score(lx, le, rec_xyz, rec_e)
            tms.append((time.perf_counter() - ta) * 1000.0)
            try:
                vs.append(float(r["vina_score"]))
            except ValueError:
                continue
            ss.append(out["score"])
        if len(vs) < 20:
            raise SystemExit("可用样本过少 (%d)" % len(vs))
        print("  缺位姿/解析失败 %d 个; 有效 %d 个; 总耗时 %.1f s" % (miss, len(vs), time.time() - t0))
        print("  单分子打分耗时: 中位 %.2f ms  最大 %.2f ms  (目标 <50 ms)"
              % (st.median(tms), max(tms)))
        rho = _spearman(ss, [-v for v in vs])   # Vina 越低越好 -> 取负使其与"分越高越好"同向
        print("  Spearman(代理分, -Vina) = %+.3f" % rho)
        print("  说明: 代理分与'更负的 Vina'应正相关; |rho| 越大越可用作筛选器。")
        print("  判读纪律: rho 明显偏低时**不得**用该代理做筛选; 它永远不能替代 A/B 主指标。")
        return 0

    if not args.pose:
        raise SystemExit("需要 --pose 或 --validate")
    lx, le = load_pose(args.pose)
    if lx is None:
        raise SystemExit("位姿无原子: %s" % args.pose)
    t = time.perf_counter()
    out = score(lx, le, rec_xyz, rec_e)
    dt = (time.perf_counter() - t) * 1000.0
    for k in ("score", "buried", "contact", "hydro", "polar", "clash", "n_lig", "n_pair"):
        print("  %-8s %s" % (k, out[k]))
    print("  耗时 %.2f ms" % dt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
