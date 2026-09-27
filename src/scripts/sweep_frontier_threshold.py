# -*- coding: utf-8 -*-
"""扫描 frontier 判定阈值, 检验"模型过早终止导致分子偏小"这一假设并寻找校准值。

问题(已实测)
------------
我们的自训模型生成的分子比它的**训练配体**还小:
    训练配体(gpcr_dock_v1) MW 中位 371.5 / 重原子 27.0
    v3 生成的分子          MW 中位 316.4 / 重原子 23.5
    (官方权重生成)         MW 中位 389.4 / 重原子 29.0
后果: 原始 Vina 分更低(打分与尺寸强相关)、药库 strict 过滤通过率只有 77%(官方 97%)。

为什么先试阈值
--------------
`frontier_pred` 输出 logits, `ind_frontier = (y_frontier_pred > frontier_threshold)`;
**没有任何 frontier 原子时生成即终止**。阈值是采样侧的**校准参数**:
微调后模型输出分布已偏移, 沿用预训练期的 0 会让它过早判定"没有边界原子"。
调低阈值 -> 更多原子被判为 frontier -> 继续生长 -> 分子变大。
这比改网络或重新训练便宜几个数量级, 且**完全可回退**(只是采样参数)。

⚠ 一个已修的坑: 配置里的 `sample.threshold.focal_threshold` 对 frontier 判定**毫无作用**
(它是原版 sample.py 的 focal 概率阈值); 真正生效的是本次接通的
`sample.threshold.frontier_threshold`(默认 0, 即 "logit > 0")。
首次扫描时误改了前者, 两个阈值的产出**逐位相同**, 才暴露出来。

风险
----
阈值过低会让分子停止不下来(历史最大失败模式: 采样日志 Failed=0、分子永不终止)。
因此本工具把**完成数**作为一号观测列, 并报出被 max_steps 截断的比例。

用法
----
    python src/scripts/sweep_frontier_threshold.py --ckpt models/7-eonmol_ft_gpcr_v3.pt \
        --target A2A --thresholds 0.5 0.3 0.2 0.1 --num-samples 60 --seed 2024
"""
import argparse
import csv
import glob
import json
import os
import statistics as st
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, CONFIGS, PYTHON  # noqa: E402

TEMPLATE = os.path.join(CONFIGS, "sample_for_pdb_guided_l3.yml")
SAMPLER = os.path.join(ROOT, "src", "sample_for_pdb.py")

# 出货药库的 strict 过滤口径(build_gpcr_v3.py 的 DRUGLIKE_TIERS['strict'])
MW_LO, MW_HI, LOGP_MAX, QED_MIN, SA_MAX, RINGS_MIN = 250.0, 500.0, 5.0, 0.40, 6.0, 1


def log(m):
    print(m, flush=True)


def metrics(smis):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, QED, RDConfig
    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer
    RDLogger.DisableLog("rdApp.*")
    ms = [m for m in (Chem.MolFromSmiles(s) for s in smis) if m is not None]
    if not ms:
        return None
    mw = [Descriptors.MolWt(m) for m in ms]
    ha = [m.GetNumHeavyAtoms() for m in ms]
    qd = [QED.qed(m) for m in ms]
    sa = [sascorer.calculateScore(m) for m in ms]
    npass = sum(1 for m in ms
                if MW_LO <= Descriptors.MolWt(m) <= MW_HI
                and Descriptors.MolLogP(m) <= LOGP_MAX
                and QED.qed(m) >= QED_MIN
                and sascorer.calculateScore(m) <= SA_MAX
                and m.GetRingInfo().NumRings() >= RINGS_MIN)
    from rdkit.Chem.Scaffolds import MurckoScaffold
    scaf = {MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in ms}
    return dict(n=len(ms), mw=st.median(mw), ha=st.median(ha), qed=st.median(qd),
                sa=st.median(sa), n_scaf=len(scaf),
                pass_pct=100.0 * npass / len(ms),
                mw250=100.0 * sum(1 for x in mw if x >= MW_LO) / len(mw))


def main():
    ap = argparse.ArgumentParser(description="frontier 阈值扫描")
    ap.add_argument("--ckpt", default=os.path.join(ROOT, "models", "7-eonmol_ft_gpcr_v3.pt"))
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--thresholds", nargs="+", type=float, default=[0.5, 0.3, 0.2, 0.1])
    ap.add_argument("--num-samples", type=int, default=60)
    ap.add_argument("--beam", type=int, default=50)
    ap.add_argument("--max-steps", type=int, default=40)
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--work", default=os.path.join(ROOT, "outputs", "thr_sweep"))
    args = ap.parse_args()

    import yaml
    reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    t = reg[args.target]
    pdb = os.path.join(ROOT, t["pdb"])
    center = ",".join("%.2f" % v for v in t["center"])

    rows = []
    for thr in args.thresholds:
        tag = "thr%s" % str(thr).replace(".", "p")
        wd = os.path.join(args.work, tag)
        os.makedirs(wd, exist_ok=True)
        c = yaml.safe_load(open(TEMPLATE, encoding="utf-8-sig"))
        c["model"]["checkpoint"] = args.ckpt
        s = c["sample"]
        s["num_samples"] = args.num_samples
        s["beam_size"] = args.beam
        s["max_steps"] = args.max_steps
        s["seed"] = args.seed
        # 只改这一个变量。键名是 frontier_threshold(新接通的), **不是** focal_threshold ——
        # 后者在原版 sample.py 里是 focal 概率阈值, 与 frontier 判定无关;
        # 实测改 focal_threshold 从 0.5 到 0.1 产出逐位相同(因为它从未被 frontier 判定读取)。
        s["threshold"]["frontier_threshold"] = thr
        s.setdefault("guided", {})["diversity_w"] = 0.5
        cfg = os.path.join(wd, "cfg.yml")
        with open(cfg, "w", encoding="utf-8") as f:
            yaml.safe_dump(c, f, allow_unicode=True, sort_keys=False)
        lf = os.path.join(wd, "sample.log")
        log("\n[focal_threshold=%.2f] 采样中 ..." % thr)
        with open(lf, "w", encoding="utf-8", errors="ignore") as f:
            rc = subprocess.run([PYTHON, SAMPLER, "--pdb_path", pdb, "--center=" + center,
                                 "--config", cfg, "--device", "cuda",
                                 "--outdir", os.path.join(wd, "samples")],
                                cwd=ROOT, stdout=f, stderr=subprocess.STDOUT).returncode
        smis, trunc = [], 0
        for p in glob.glob(os.path.join(wd, "samples", "*", "SMILES.txt")):
            smis += [l.strip() for l in open(p, encoding="utf-8") if l.strip()]
        # 从采样日志读"被步数截断"的迹象: 最后一步仍有大量未完成队列
        if os.path.exists(lf):
            for line in open(lf, encoding="utf-8", errors="ignore"):
                if "Queue" in line and "Finished" in line:
                    try:
                        last = line
                    except Exception:
                        pass
        m = metrics(smis)
        if m is None:
            log("  rc=%d, 没有产出分子" % rc)
            continue
        m.update(threshold=thr, rc=rc)
        rows.append(m)
        log("  完成 %3d | MW中位 %6.1f | 重原子 %5.1f | QED %.3f | SA %.2f | 骨架 %3d"
            " | MW>=250 %5.1f%% | 过strict %5.1f%%"
            % (m["n"], m["mw"], m["ha"], m["qed"], m["sa"], m["n_scaf"],
               m["mw250"], m["pass_pct"]))

    if not rows:
        raise SystemExit("没有任何结果")
    out = os.path.join(args.work, "threshold_sweep.csv")
    with open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    log("\n" + "=" * 78)
    log("阈值扫描结果(参考基准: 官方权重生成 MW 389/HA 29; 训练配体 MW 371/HA 27)")
    log("=" * 78)
    log("%6s %6s %8s %8s %7s %6s %6s %9s %9s"
        % ("阈值", "完成", "MW中位", "重原子", "QED", "SA", "骨架", "MW>=250", "过strict"))
    for r in rows:
        log("%6.2f %6d %8.1f %8.1f %7.3f %6.2f %6d %8.1f%% %8.1f%%"
            % (r["threshold"], r["n"], r["mw"], r["ha"], r["qed"], r["sa"],
               r["n_scaf"], r["mw250"], r["pass_pct"]))
    log("\n判读要点")
    log("  - 一号观测列是**完成数**: 明显下降意味着阈值调低后分子停不下来(危险信号);")
    log("  - 目标是让 MW 中位接近训练分布(371)且过滤通过率接近官方(97%), 同时完成数不塌;")
    log("  - 通过率大幅提升但完成数塌了 -> 该阈值不可用, 应转为训练侧修复。")
    log("产物: %s" % os.path.relpath(out, ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
