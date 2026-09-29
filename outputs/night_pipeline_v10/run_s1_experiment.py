# -*- coding: utf-8 -*-
"""S1（轨迹内尺寸调度）预登记实验：**同一批内跑"关/开"两臂**，协议完全一致。

为什么不用既有臂当对照：盘上 `outputs/ab_v10_s60/*` 是 **max_steps 60 / n=50** 的协议，
与预登记的 40/60 不同 → **不能当配对对照**（本仓库的血泪纪律：协议不一致的对比出过事故）。
故本脚本每种子各跑两臂，除 `--schedule-*` 外**所有参数逐字相同**。

预登记（记忆 §5.46）：
- 权重 v10@5500（当前默认）；靶点 A2A；n=60 / beam 50 / max_steps 40；阈值 0.0；种子 2024–2027。
- S1 参数：target_ha=28 / slack=4.0 / relax=0.35 / warmup=0.25。
- 判据：① 落区率 `340≤MW≤430` ≥3/4 种子；② 跨种子极差 ≤10%；③ 生成数 ≥n；
  ④ 后续对接级 2 种子不劣于官方（±0.20）且尺寸差 <15%。
- 成本 ≈ 8 次采样 × 60 分子 ≈ 20–30 分钟 GPU（**比 V2 训练便宜 ~3–6 倍**）。

用法：
    python outputs/night_pipeline_v10/run_s1_experiment.py            # 跑实验（GPU）
    python outputs/night_pipeline_v10/run_s1_experiment.py --dry-run  # 只打印命令
"""
import argparse
import csv
import glob
import io
import os
import subprocess
import sys

# 无人值守时 stdout 常是 GBK 管道：非 GBK 字符会直接崩脚本（A′ v1.0 真实踩过）
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PY = sys.executable
BASE = os.path.join(ROOT, "outputs", "s1_exp")
BAND = (340.0, 430.0)


def log(msg, path):
    line = "[%s] %s" % (__import__("time").strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with io.open(path, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def run_design(seed, arm, ckpt, n, beam, steps, thr, dry):
    tag = {"on": "s1on", "off": "s1off", "official": "official"}[arm]
    outdir = os.path.join(BASE, "%s_s%d" % (tag, seed))
    cmd = [PY, os.path.join(ROOT, "design.py"), "--target", "A2A",
           "--num-samples", str(n), "--beam", str(beam), "--max-steps", str(steps),
           "--seed", str(seed), "--ckpt", ckpt,
           "--frontier-threshold", str(thr), "--outdir", outdir]
    if arm == "on":
        cmd += ["--schedule-target-ha", "28", "--schedule-slack", "4",
                "--schedule-relax", "0.35", "--schedule-warmup", "0.25"]
    if dry:
        print(" ".join(cmd))
        return outdir, 0
    os.makedirs(BASE, exist_ok=True)
    with io.open(os.path.join(BASE, "run_%s_s%d.log" % (tag, seed)), "w",
                 encoding="utf-8") as fh:
        rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
    return outdir, rc


def find_smiles(d):
    c = glob.glob(os.path.join(d, "**", "SMILES.txt"), recursive=True)
    return sorted(c)[-1] if c else None


def mannwhitney_u(x, y):
    """小样本 U 检验（正态近似，双侧）—— 两臂跨种子尺寸是否有差。"""
    import math
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return None, None
    ranks = {}
    allv = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    i = 0
    order = {}
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            order.setdefault(allv[k][1], []).append(avg)
        i = j + 1
    r1 = sum(order.get(0, []))
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = (u1 - mu) / sd if sd else 0.0
    p = math.erfc(abs(z) / math.sqrt(2))
    return u1, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", nargs="+", type=int, default=[2024, 2025, 2026, 2027])
    ap.add_argument("--ckpt", default=os.path.join(ROOT, "models", "7-eonmol_ft_gpcr_v10.pt"))
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--beam", type=int, default=50)
    ap.add_argument("--max-steps", type=int, default=40)
    ap.add_argument("--thr", type=float, default=0.0)
    ap.add_argument("--official-seeds", nargs="+", type=int, default=[2024, 2025],
                    help="额外跑**同协议官方臂**的种子（判据④需要；默认前两个种子）")
    ap.add_argument("--official-ckpt",
                    default=os.path.join(ROOT, "models", "pretrained_Pocket2Mol.pt"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    logf = os.path.join(BASE, "s1_experiment.log")
    os.makedirs(BASE, exist_ok=True)
    log("S1 实验开始 | ckpt=%s | seeds=%s | n=%d beam=%d steps=%d thr=%.2f"
        % (os.path.basename(args.ckpt), args.seeds, args.n, args.beam, args.max_steps, args.thr), logf)

    arms = []
    plan = [(seed, arm) for seed in args.seeds for arm in ("off", "on")]
    plan += [(seed, "official") for seed in args.official_seeds]   # 同协议官方臂
    for seed, arm in plan:
        ck = args.official_ckpt if arm == "official" else args.ckpt
        d, rc = run_design(seed, arm, ck, args.n, args.beam, args.max_steps,
                           args.thr, args.dry_run)
        log("  seed %d | %-8s | rc=%d | %s" % (seed, arm, rc, os.path.relpath(d, ROOT)), logf)
        if not args.dry_run and rc == 0:
            smi = find_smiles(d)
            if smi:
                tag = {"on": "s1on", "off": "s1off", "official": "official"}[arm]
                arms.append(("%s_s%d" % (tag, seed), smi))
    if args.dry_run:
        return 0

    out_csv = os.path.join(BASE, "arm_gen_s1.csv")
    cmd = [PY, os.path.join(ROOT, "src", "scripts", "arm_gen_metrics.py"), "--arms"] + \
          ["%s=%s" % a for a in arms] + ["--baseline", "s1off_s%d" % args.seeds[0],
                                         "--out", out_csv]
    with io.open(os.path.join(BASE, "arm_gen_s1.log"), "w", encoding="utf-8") as fh:
        subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    if not os.path.exists(out_csv):
        log("arm_gen_metrics 未产出 CSV，终止", logf)
        return 1

    rows = {r["arm"]: r for r in csv.DictReader(io.open(out_csv, encoding="utf-8-sig"))}
    print("\n%-12s %8s %7s %8s %8s %8s" % ("arm", "MW中位", "HA中位", "strict%", "生成数", "落区"))
    mws = {"off": [], "on": []}
    inband = {"off": 0, "on": 0}
    for seed in args.seeds:
        for arm in ("off", "on"):
            k = "s1%s_s%d" % (arm, seed)
            r = rows.get(k)
            if not r:
                continue
            mw = float(r["mw_med"]); ha = float(r["ha_med"])
            sr = 100 * float(r["strict_rate"]); nf = int(float(r["n_finished"]))
            ok = BAND[0] <= mw <= BAND[1]
            inband[arm] += int(ok)
            mws[arm].append(mw)
            print("%-14s %8.1f %7.1f %8.1f %8d %8s" % (k, mw, ha, sr, nf, "IN" if ok else "--"))
    print()
    for arm in ("off", "on"):
        v = mws[arm]
        if v:
            spread = 100.0 * (max(v) - min(v)) / (sum(v) / len(v))
            print("S1 %-3s: 落区 %d/%d | 跨种子极差 %.1f%% | MW 均值 %.1f"
                  % (arm, inband[arm], len(v), spread, sum(v) / len(v)))
    if len(mws["off"]) == len(mws["on"]) and mws["off"]:
        u, p = mannwhitney_u(mws["off"], mws["on"])
        print("配对比较（种子级 MW）：U=%.1f, p≈%.3f（n=%d/臂，仅作参考——真正的判据是落区率与极差）"
              % (u, p, len(mws["off"])))
        verdict = ("PASS" if (inband["on"] >= max(3, len(mws["on"]) * 3 // 4)
                              and 100.0 * (max(mws["on"]) - min(mws["on"])) / (sum(mws["on"]) / len(mws["on"])) <= 10.0)
                   else "NOT_PASSED")
        print("\n预登记判据 ①②（生成级）：**%s** —— ③④（生成数/对接）见 arm CSV 与后续对接步骤" % verdict)
    return 0


if __name__ == "__main__":
    sys.exit(main())
