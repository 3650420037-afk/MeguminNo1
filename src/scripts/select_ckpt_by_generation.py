# -*- coding: utf-8 -*-
"""按「生成质量」而非「掩码补全损失」挑选检查点。

为什么需要它
------------
训练时的验证损失是**教师强制下的掩码补全**指标, 而评测看的是**无教师强制的
自回归生成**质量, 二者可以反向 —— 实测出现过"验证损失从 1.97 降到 0.77,
生成 QED 却从 0.597 掉到 0.408"。因此最终选点必须用生成指标, 不能只看 val loss。

做法
----
1. 从训练日志目录读取 checkpoints/val_history.jsonl, 按 val loss 选 top-K 候选
   (也可用 --ckpts 显式指定);
2. 对每个候选, 用**完全相同**的采样参数(靶点/种子/分子数/束宽/步数)跑一次生成;
3. 对生成结果计算:
   - 完成分子数 / 唯一规范化 SMILES 数
   - QED 中位、SA 中位
   - Murcko 骨架数
   - **口袋叠合度**: 生成分子与口袋蛋白的最近重原子距离中位(判据同
     diagnose_dataset.py: <0.5 A 记为穿模) —— 这是本项目新增的评估轴,
     官方权重从未被测过这一项;
4. 按 --rank-by 排序输出 CSV 与推荐检查点。

用法
----
    python src/scripts/select_ckpt_by_generation.py \
        --run logs/train_gpcr_ft_2026_xx_xx__xx_xx_xx \
        --target A2A --topk 3 --num-samples 10 --beam 50 --max-steps 30

    # 显式指定检查点
    python src/scripts/select_ckpt_by_generation.py \
        --ckpts logs/<run>/checkpoints/1200.pt logs/<run>/checkpoints/best.pt \
        --target A2A
"""
import argparse
import csv
import glob
import json
import os
import statistics as st
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import (ROOT, CONFIGS, DATA, OUTPUTS, TARGETS, ensure_dir)  # noqa: E402
from paths import src_on_path  # noqa: E402
src_on_path()

CLASH_MAX = 0.5      # 与 diagnose_dataset.py 的 LIGAND_CLASH_POCKET 同判据
DIVERSITY_W = 0.5    # 与 configs/sample_for_pdb_guided_l3.yml 及 design.py/predict.py
                     # 的 --diversity-w 默认值保持一致, 保证 A/B 与出货目标同源。


def pick_candidates(run_dir, topk):
    """按 val loss 从 val_history.jsonl 里取 top-K 检查点。"""
    hist = os.path.join(run_dir, "checkpoints", "val_history.jsonl")
    if not os.path.exists(hist):
        return []
    rows = []
    for line in open(hist, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    rows.sort(key=lambda r: r.get("val_total", 1e9))
    out = []
    for r in rows[:topk]:
        p = os.path.join(run_dir, "checkpoints", "%d.pt" % r["iter"])
        if os.path.exists(p):
            out.append((p, r.get("val_total"), r.get("iter")))
    return out


def evaluate_ckpt(ckpt, target, num_samples, beam, max_steps, seed, device, workdir,
                  frontier_threshold=0.0):
    """跑一次生成并对产物计算指标; 返回 dict。"""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import QED, AllChem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")

    sys.path.insert(0, os.path.join(ROOT, "src"))
    from scripts.diagnose_dataset import min_distances
    import numpy as np

    # 1) 用该检查点生成
    cfg = os.path.join(workdir, "cfg.yml")
    import yaml
    tpl = os.path.join(CONFIGS, "sample_for_pdb_guided_l3.yml")
    c = yaml.safe_load(open(tpl, encoding="utf-8-sig"))
    c["model"]["checkpoint"] = ckpt
    c["sample"]["seed"] = seed
    c["sample"]["num_samples"] = num_samples
    c["sample"]["beam_size"] = beam
    c["sample"]["max_steps"] = max_steps
    # 评测目标必须与出货目标一致: 显式固定引导项, 不依赖模板是否恰好写了这些键。
    # 历史缺陷: 本函数只覆盖 ckpt/seed/num_samples/beam/max_steps, diversity_w 沿用
    # 模板 -> 模板缺键时取代码默认 0.0, 于是 A/B 在"无多样性惩罚"下评测, 而药库在
    # 0.5 下生成, 评测目标 != 出货目标。
    g = c["sample"].setdefault("guided", {})
    g["enabled"] = True
    g.setdefault("qed_w", 1.0)
    g.setdefault("sa_w", 1.0)
    g.setdefault("lam", 3)
    g["diversity_w"] = float(DIVERSITY_W)
    # frontier 判定阈值: 决定"何时停", 直接控制分子大小。
    # 键名是 frontier_threshold(新接通的), **不是** focal_threshold ——
    # 后者是原版 sample.py 的 focal 概率阈值, 对 frontier 判定完全无效
    # (实测把 focal_threshold 从 0.5 改到 0.1 产出逐位相同)。
    c["sample"].setdefault("threshold", {})["frontier_threshold"] = float(frontier_threshold)
    yaml.safe_dump(c, open(cfg, "w", encoding="utf-8"), allow_unicode=True, sort_keys=False)

    reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    t = reg[target]
    pdb = os.path.join(ROOT, t["pdb"])
    outdir = ensure_dir(os.path.join(workdir, "samples"))
    cmd = [sys.executable, os.path.join(ROOT, "src", "sample_for_pdb.py"),
           "--pdb_path", pdb,
           "--center=" + ",".join("%.2f" % v for v in t["center"]),
           "--config", cfg, "--device", device, "--outdir", outdir]
    with open(os.path.join(workdir, "sample.log"), "w", encoding="utf-8") as lf:
        rc = subprocess.run(cmd, cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT).returncode

    sess = sorted(d for d in glob.glob(os.path.join(outdir, "*")) if os.path.isdir(d))
    res = {"ckpt": os.path.basename(ckpt), "returncode": rc, "n_finished": 0,
           "n_unique": 0, "qed_med": None, "sa_med": None, "n_scaffold": 0,
           "clash_rate": None, "mind_med": None, "mind_min": None}
    if not sess:
        return res
    sess = sess[-1]

    smis = []
    p = os.path.join(sess, "SMILES.txt")
    if os.path.exists(p):
        smis = sorted({s.strip() for s in open(p, encoding="utf-8", errors="ignore") if s.strip()})
    mols = [m for m in (Chem.MolFromSmiles(s) for s in smis) if m is not None]
    res["n_finished"] = len(smis)
    res["n_unique"] = len({Chem.MolToSmiles(m) for m in mols})
    if mols:
        try:
            from evaluation.sascorer import calculateScore as sa_calc
            res["qed_med"] = round(st.median([QED.qed(m) for m in mols]), 4)
            res["sa_med"] = round(st.median([sa_calc(m) for m in mols]), 4)
        except Exception:
            res["qed_med"] = round(st.median([QED.qed(m) for m in mols]), 4)
        res["n_scaffold"] = len({MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols})

    # 2) 口袋叠合度: 生成分子的 3D 坐标 vs 口袋蛋白原子
    from utils.protein_ligand import PDBProtein
    sdf_dir = os.path.join(sess, "SDF")
    sdfs = sorted(glob.glob(os.path.join(sdf_dir, "*.sdf"))) if os.path.isdir(sdf_dir) else []
    if sdfs:
        # 口袋: 与 sample_for_pdb 同样的半径口径(12 A, center_of_mass)
        prot = PDBProtein(pdb)
        residues = prot.query_residues_radius(t["center"], 12.0, criterion="center_of_mass")
        pk = os.path.join(workdir, "pocket.pdb")
        with open(pk, "w", encoding="ascii") as h:
            h.write(prot.residues_to_pdb_block(residues, name="POCKET"))
        pocket = PDBProtein(pk)
        ppos = np.asarray(pocket.pos, dtype=float)
        minds = []
        for f in sdfs[:200]:
            mol = Chem.MolFromMolFile(f, sanitize=False, removeHs=True)
            if mol is None or mol.GetNumConformers() == 0:
                continue
            c = mol.GetConformer()
            lp = np.asarray([[c.GetAtomPosition(i).x, c.GetAtomPosition(i).y,
                              c.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())],
                            dtype=float)
            _, mind = min_distances(lp, ppos)
            if mind is not None:
                minds.append(mind)
        if minds:
            a = np.asarray(minds)
            res["clash_rate"] = round(float((a < CLASH_MAX).mean()), 4)
            res["mind_med"] = round(float(np.median(a)), 4)
            res["mind_min"] = round(float(a.min()), 4)
    return res


def main():
    ap = argparse.ArgumentParser(description="按生成质量挑选检查点")
    ap.add_argument("--run", help="训练日志目录(内含 checkpoints/)")
    ap.add_argument("--ckpts", nargs="+", help="显式指定检查点(与 --run 二选一)")
    ap.add_argument("--topk", type=int, default=3, help="按 val loss 取前 K 个候选")
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--num-samples", type=int, default=10)
    ap.add_argument("--beam", type=int, default=50)
    ap.add_argument("--max-steps", type=int, default=30)
    ap.add_argument("--seed", type=int, default=2024, help="所有候选必须用同一种子")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--frontier-threshold", type=float, default=0.0,
                    help="frontier 判定阈值(决定何时停, 直接控制分子大小); 默认 0 与官方行为一致。"
                         "注意与 focal_threshold 不同 —— 后者是原版 sample.py 的 focal 概率阈值, "
                         "对 frontier 判定完全无效。")
    ap.add_argument("--rank-by", default="qed",
                    choices=["qed", "clash", "both"],
                    help="排序依据: qed=QED 中位降序; clash=穿模率升序; both=先穿模率后 QED")
    ap.add_argument("--keep-dir", default=None,
                    help="给定目录时保留每个检查点的采样产物到 <dir>/<检查点名>/, "
                         "便于再用 src/scripts/judge_promotion.py 复核(默认用临时目录不保留)")
    ap.add_argument("--out", default=os.path.join(OUTPUTS, "ckpt_selection.csv"))
    args = ap.parse_args()

    if args.ckpts:
        cands = [(p, None, None) for p in args.ckpts]
    elif args.run:
        cands = pick_candidates(args.run, args.topk)
    else:
        raise SystemExit("需要 --run 或 --ckpts")
    if not cands:
        raise SystemExit("没有可用检查点(检查 checkpoints/val_history.jsonl 是否存在)")

    print("=" * 88)
    print("按生成质量挑选检查点 | 靶点=%s 种子=%d 规模=%d/%d/%d 候选=%d"
          % (args.target, args.seed, args.num_samples, args.beam, args.max_steps, len(cands)))
    try:
        from paths import recommended_frontier_threshold as _rec
        for ck, _v, _it in cands:
            _val, _src = _rec(ck)
            if abs(_val - args.frontier_threshold) > 1e-9:
                print("[口径提示] %s 的推荐阈值是 %+.2f（来源 %s），本次用 %+.2f"
                      % (os.path.basename(ck), _val, _src, args.frontier_threshold), flush=True)
    except Exception:
        pass
    print("=" * 88, flush=True)

    rows = []
    for ckpt, vloss, it in cands:
        if not os.path.exists(ckpt):
            print("  跳过(不存在): %s" % ckpt)
            continue
        print("  [评估] %s (val_loss=%s)" % (os.path.basename(ckpt), vloss), flush=True)
        if args.keep_dir:
            # 保留产物: 每个检查点一个目录, 内含唯一一次采样的会话目录,
            # 可直接交给 src/scripts/judge_promotion.py 做晋级判定
            wd = ensure_dir(os.path.join(args.keep_dir, os.path.splitext(os.path.basename(ckpt))[0]))
            r = evaluate_ckpt(ckpt, args.target, args.num_samples, args.beam,
                              args.max_steps, args.seed, args.device, wd,
                              args.frontier_threshold)
        else:
            with tempfile.TemporaryDirectory(prefix="ckptsel_") as wd:
                r = evaluate_ckpt(ckpt, args.target, args.num_samples, args.beam,
                                  args.max_steps, args.seed, args.device, wd,
                                  args.frontier_threshold)
        r["val_loss"] = vloss
        r["iter"] = it
        r["run_dir"] = os.path.relpath(
            os.path.join(args.keep_dir, os.path.splitext(os.path.basename(ckpt))[0]), ROOT
        ) if args.keep_dir else ""
        rows.append(r)
        print("        完成 %d | 唯一 %d | QED中位 %s | SA中位 %s | 骨架 %d | "
              "穿模率 %s | 最近距离中位 %s"
              % (r["n_finished"], r["n_unique"], r["qed_med"], r["sa_med"],
                 r["n_scaffold"], r["clash_rate"], r["mind_med"]), flush=True)

    if not rows:
        raise SystemExit("没有成功评估任何检查点")

    def key(r):
        q = r["qed_med"] if r["qed_med"] is not None else -1
        c = r["clash_rate"] if r["clash_rate"] is not None else 9
        if args.rank_by == "qed":
            return (-q, c)
        if args.rank_by == "clash":
            return (c, -q)
        return (c, -q)

    rows.sort(key=key)
    ensure_dir(os.path.dirname(args.out))
    with open(args.out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    best = rows[0]
    print("=" * 88)
    print("推荐检查点: %s (QED中位 %s, 穿模率 %s, 骨架 %d)"
          % (best["ckpt"], best["qed_med"], best["clash_rate"], best["n_scaffold"]))
    print("明细: %s" % os.path.relpath(args.out, ROOT))
    if args.keep_dir:
        print()
        print("产物已保留, 可用项目原生判据复核(注意: 每个 run_dir 下只应有这一次采样):")
        print("  python src/scripts/judge_promotion.py <baseline_run_dir> <finetuned_run_dir> "
              "outputs/decision.txt")
        for r in rows:
            print("    %-14s %s" % (r["ckpt"], r["run_dir"]))
    print("=" * 88)
    return 0


if __name__ == "__main__":
    sys.exit(main())
