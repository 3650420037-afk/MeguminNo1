# -*- coding: utf-8 -*-
"""M7：奖励代理的**分布迁移**验证 —— 在【模型自生成位姿】上复核折外 Spearman。

为什么必须单独验证（记忆文件 §5.9）
------------------------------------
`learn_reward_proxy.py` 报出的 ρ=0.60 是在 **Vina 优化过的位姿**上算的；
而 ReST 循环要筛的是**模型自己生成的位姿**——两者分布不同，不验证就用 = 未经验证下结论。

关键设计：用 `vina --score_only`，**不要重新对接**
--------------------------------------------------
若把生成分子重新送进 `docking_pipeline.py`，Vina 会**重新搜索**并给出它自己的位姿，
那测到的仍然只是"Vina 优化位姿"上的代理质量，迁移问题根本没被触及。
`--score_only` 在**给定位姿上只算分、不搜索**，才回答得了"这个姿态好不好"。
（实测：同一生成位姿 score_only = −3.10 kcal/mol，而它在对接管线里的分数是 −9 量级 —— 差距正是迁移风险的来源。）

做法
----
1. 用对接数据训练 ridge 代理（口径与 `learn_reward_proxy.py` 完全一致），报**域内**折外 ρ；
2. 对每个生成位姿 SDF：`obabel -isdf ... -opdbqt`（**保坐标**，不加 --gen3d）
   → `vina --score_only`（带口袋 center/size）→ 同口径特征；
3. 报**域外（生成位姿）** ρ。判读：若 |ρ| 明显低于域内（经验阈值 <0.3），
   则代理**不可用于 ReST**，如实判失败；即便可用，也**不得替代 A/B 主指标**。

用法
----
    python src/scripts/verify_reward_proxy_migration.py \
        --gen-dirs outputs/ab_v9_gate/s2024 outputs/ab_combo/thrm05 \
        --target A2A --limit-gen 400 --workers 8 \
        --out outputs/m7_proxy/report_A2A.json
"""
import argparse
import csv
import glob
import json
import os
import random
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from paths import ROOT, CONFIGS, OUTPUTS, ensure_dir  # noqa: E402
from fast_reward_proxy import load_receptor, load_pose, _spearman  # noqa: E402
from learn_reward_proxy import features, ridge_fit  # noqa: E402
from docking_pipeline import prepare_receptor  # noqa: E402

OBABEL = os.path.join(ROOT, "tools", "openbabel", "Library", "bin", "obabel.exe")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
SIZE = 22.5
SCORE_RE = re.compile(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.?\d*)")


def target_ctx(target):
    reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    t = reg[target]
    return os.path.join(ROOT, t["pdb"]), [float(x) for x in t["center"]]


def lig_box(pdbqt):
    """从 PDBQT 读重原子坐标, 返回 (质心, 各轴最大跨度+8Å)。

    `--score_only` 只对**给定位姿**算分, 因此盒子只需包住该配体即可 ——
    用固定口袋盒会大量报 "The ligand is outside the grid box"（实测 10/12 失败）。
    盒子越小, Vina 计算网格越快, 所以按配体自适应。
    """
    xs, ys, zs = [], [], []
    for line in open(pdbqt, encoding="utf-8", errors="ignore"):
        if line.startswith(("ATOM", "HETATM")):
            try:
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
            except ValueError:
                continue
    if not xs:
        return None
    cx, cy, cz = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2
    sx = max(8.0, max(xs) - min(xs) + 8.0)
    sy = max(8.0, max(ys) - min(ys) + 8.0)
    sz = max(8.0, max(zs) - min(zs) + 8.0)
    return (cx, cy, cz), (sx, sy, sz)


def score_one(job):
    """SDF -> PDBQT（保坐标）-> vina --score_only。返回 (sdf, score|None, pdbqt|None)。"""
    sdf, work, center = job
    stem = os.path.splitext(os.path.basename(sdf))[0]
    parent = os.path.basename(os.path.dirname(os.path.dirname(sdf)))
    pdbqt = os.path.join(work, "ligs", "%s_%s.pdbqt" % (parent, stem))
    if not os.path.exists(pdbqt):
        r = subprocess.run([OBABEL, "-isdf", sdf, "-opdbqt", "-O", pdbqt],
                           capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(pdbqt) or os.path.getsize(pdbqt) == 0:
            return sdf, None, None
    box = lig_box(pdbqt)
    if box is None:
        return sdf, None, None
    c, sz = box
    cmd = [VINA, "--receptor", os.path.join(work, "receptor.pdbqt"), "--ligand", pdbqt,
           "--score_only",
           "--center_x=%.3f" % c[0], "--center_y=%.3f" % c[1], "--center_z=%.3f" % c[2],
           "--size_x=%.1f" % sz[0], "--size_y=%.1f" % sz[1], "--size_z=%.1f" % sz[2]]
    r = subprocess.run(cmd, capture_output=True, text=True)
    m = SCORE_RE.search(r.stdout or "")
    if not m:
        return sdf, None, None
    return sdf, float(m.group(1)), pdbqt


def main():
    ap = argparse.ArgumentParser(description="奖励代理的分布迁移验证（vina --score_only）")
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--train-summary", default=None, help="对接 summary_merged.csv（默认 outputs/mass_dock/<T>/summary_merged.csv）")
    ap.add_argument("--train-poses", default=None, help="对接位姿目录（默认 outputs/mass_dock/<T>/docking）")
    ap.add_argument("--gen-dirs", nargs="+", required=True, help="生成位姿 SDF 的搜索根目录（递归）")
    ap.add_argument("--limit-train", type=int, default=3000)
    ap.add_argument("--limit-gen", type=int, default=400)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(OUTPUTS, "m7_proxy", "report.json"))
    args = ap.parse_args()

    summary = args.train_summary or os.path.join(OUTPUTS, "mass_dock", args.target, "summary_merged.csv")
    poses_dir = args.train_poses or os.path.join(OUTPUTS, "mass_dock", args.target, "docking")
    work = ensure_dir(os.path.join(OUTPUTS, "m7_proxy", args.target))
    ensure_dir(os.path.join(work, "ligs"))
    ensure_dir(os.path.join(work, "logs"))   # prepare_receptor 会往 <work>/logs 写日志

    rec_xyz, rec_e = load_receptor(args.target)
    rows = [r for r in csv.DictReader(open(summary, encoding="utf-8-sig"))
            if r.get("status") == "ok"][:args.limit_train]
    X, y = [], []
    for r in rows:
        p = os.path.join(poses_dir, r["name"] + "_out.pdbqt")
        if not os.path.exists(p):
            continue
        lx, le = load_pose(p)
        if lx is None:
            continue
        try:
            v = float(r["vina_score"])
        except ValueError:
            continue
        X.append(features(lx, le, rec_xyz, rec_e))
        y.append(-v)
    X = np.asarray(X)
    y = np.asarray(y)
    print("域内训练样本 %d（对接优化位姿）, 特征维度 %d" % (len(X), X.shape[1]))

    # 域内折外 ρ（与 learn_reward_proxy 同口径）
    rng = np.random.default_rng(args.seed)
    idx = rng.permutation(len(X))
    folds = np.array_split(idx, args.folds)
    pred_in = np.zeros(len(X))
    for i, te in enumerate(folds):
        tr = np.concatenate([f for j, f in enumerate(folds) if j != i])
        mu, sd = X[tr].mean(axis=0), X[tr].std(axis=0) + 1e-9
        w = ridge_fit((X[tr] - mu) / sd, y[tr], args.lam)
        pred_in[te] = np.hstack([(X[te] - mu) / sd, np.ones((len(te), 1))]) @ w
    rho_in = _spearman(pred_in, y)

    # 全量拟合，用于迁移预测
    mu, sd = X.mean(axis=0), X.std(axis=0) + 1e-9
    w = ridge_fit((X - mu) / sd, y, args.lam)

    # 迁移：生成位姿
    pdb, center = target_ctx(args.target)
    rec_pdbqt, notes = prepare_receptor(pdb, work)
    print("受体: %s (%s)" % (rec_pdbqt, "; ".join(notes[:2])))
    sdfs = []
    for d in args.gen_dirs:
        for f in glob.glob(os.path.join(d, "**", "*.sdf"), recursive=True):
            parts = os.path.normpath(f).split(os.sep)
            # 只取**模型生成**的位姿: 必须位于 samples/ 下, 且不在 dock/ 里(那是重新对接过的位姿)
            if "samples" not in parts or "dock" in parts or "docking" in parts:
                continue
            sdfs.append(f)
    sdfs = sorted(set(sdfs))
    rng2 = random.Random(args.seed)
    rng2.shuffle(sdfs)
    sdfs = sdfs[:args.limit_gen]
    print("生成位姿候选 %d 个（抽样上限 %d）" % (len(sdfs), args.limit_gen))

    jobs = [(s, work, center) for s in sdfs]
    results = []
    with __import__("concurrent.futures").futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
        for sdf, sc, pdbqt in ex.map(score_one, jobs, chunksize=4):
            results.append((sdf, sc, pdbqt))
    ok = [(s, sc, p) for s, sc, p in results if sc is not None]
    print("score_only 成功 %d / %d" % (len(ok), len(results)))
    if not ok:
        raise SystemExit("没有成功的 score_only 结果")

    Xg, yg, used = [], [], []
    for sdf, sc, pdbqt in ok:
        lx, le = load_pose(pdbqt)
        if lx is None:
            continue
        Xg.append(features(lx, le, rec_xyz, rec_e))
        yg.append(-sc)
        used.append((sdf, sc))
    Xg = np.asarray(Xg)
    yg = np.asarray(yg)
    pred_out = np.hstack([(Xg - mu) / sd, np.ones((len(Xg), 1))]) @ w
    rho_out = _spearman(pred_out, yg)

    report = dict(target=args.target, n_train=len(X), n_gen=len(Xg),
                  rho_in_domain=round(float(rho_in), 4),
                  rho_transfer=round(float(rho_out), 4),
                  score_only_median=round(float(np.median([-v for v in yg])), 3),
                  verdict=("usable_as_filter" if abs(rho_out) >= 0.30 else "NOT_USABLE"),
                  gen_dirs=args.gen_dirs, limit_gen=args.limit_gen, seed=args.seed)
    ensure_dir(os.path.dirname(args.out))
    json.dump(report, open(args.out, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n" + "=" * 70)
    print("域内（Vina 优化位姿, %d 折折外）: Spearman = %+.3f" % (args.folds, rho_in))
    print("域外（模型自生成位姿, n=%d）  : Spearman = %+.3f" % (len(Xg), rho_out))
    print("生成位姿 score_only 中位 = %.3f kcal/mol" % report["score_only_median"])
    print("判读: %s" % ("|ρ|>=0.30 → 可作**筛选器**候选（仍不得替代 A/B 主指标）"
                        if abs(rho_out) >= 0.30 else
                        "|ρ|<0.30 → **代理不可用于 ReST**（按纪律如实判失败，不调参硬凑）"))
    print("报告: %s" % os.path.relpath(args.out, ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
