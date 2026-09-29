# -*- coding: utf-8 -*-
"""S3 调窗复验：把硬尺寸窗下调到 (floor 22, ceiling 30)，看能否让 MW 中位居中在 340–430。

为什么调（阶段 1 实测）：S3(24,32) 的 strict 最好（96.7%/89.0%），但 MW 中位 435–439 **略高于带上限 430**
—— 因为上限是硬截断，分子会在 32 个重原子处堆积。按"每步 +1 重原子、MW ≈ 13.5 Da/重原子"估算，
上限 30 个重原子 ≈ MW 405，落在带中央。

用法：python outputs/night_pipeline_v10/run_s3_tuned.py [--wait-pid <凹坑测试 PID>]
"""
import argparse
import glob
import os
import statistics as st
import subprocess
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PY = sys.executable
BASE = os.path.join(ROOT, "outputs", "s3_tuned")
SEEDS = (2024, 2025)
CKPT = os.path.join(ROOT, "models", "7-eonmol_ft_gpcr_v10.pt")
FLOOR, CEIL = 22, 30


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with open(os.path.join(BASE, "s3_tuned.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def pid_alive(pid):
    out = subprocess.run(["tasklist", "/FI", "PID eq %d" % pid], capture_output=True,
                         text=True, errors="ignore").stdout
    return str(pid) in out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wait-pid", type=int, default=0)
    ap.add_argument("--n", type=int, default=60)
    args = ap.parse_args()
    os.makedirs(BASE, exist_ok=True)
    if args.wait_pid:
        log("等待前置任务（PID %d）结束…" % args.wait_pid)
        t0 = time.time()
        while pid_alive(args.wait_pid):
            if time.time() - t0 > 2 * 3600:
                log("等待超时，放弃")
                return 1
            time.sleep(30)
        log("前置任务结束（等待 %.0f 分钟）" % ((time.time() - t0) / 60.0))
        time.sleep(10)

    for seed in SEEDS:
        outdir = os.path.join(BASE, "s3t_s%d" % seed)
        cmd = [PY, os.path.join(ROOT, "design.py"), "--target", "A2A",
               "--num-samples", str(args.n), "--beam", "50", "--max-steps", "40",
               "--seed", str(seed), "--ckpt", CKPT, "--frontier-threshold", "0.0",
               "--schedule-floor-steps", str(FLOOR), "--schedule-floor-thr", "-2.0",
               "--schedule-ceiling-ha", str(CEIL), "--outdir", outdir]
        t0 = time.time()
        with open(os.path.join(BASE, "run_s%d.log" % seed), "w", encoding="utf-8") as fh:
            rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
        log("  seed %d | rc=%d | %.1f 分钟" % (seed, rc, (time.time() - t0) / 60.0))

    # 汇总
    sys.path.insert(0, os.path.join(ROOT, "src"))
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    RDLogger.DisableLog("rdApp.*")
    print("%-14s %5s %9s %7s %8s %8s" % ("臂", "n", "MW中位", "HA", "落区率", "MW<340"))
    for seed in SEEDS:
        c = glob.glob(os.path.join(BASE, "s3t_s%d" % seed, "**", "SMILES.txt"), recursive=True)
        if not c:
            continue
        ms = [Chem.MolFromSmiles(l.strip()) for l in open(sorted(c)[-1], encoding="utf-8") if l.strip()]
        ms = [m for m in ms if m]
        if not ms:
            continue
        ha = sorted(m.GetNumHeavyAtoms() for m in ms)
        mw = sorted(round(Descriptors.MolWt(m), 1) for m in ms)
        print("S3(22,30) s%d  %5d %9.1f %7.1f %8.2f %8.2f"
              % (seed, len(ms), st.median(mw), st.median(ha),
                 sum(1 for x in mw if 340 <= x <= 430) / len(mw),
                 sum(1 for x in mw if x < 340) / len(mw)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
