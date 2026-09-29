# -*- coding: utf-8 -*-
"""**凹坑测试**：把尺寸控制机制用在"已知掉进小模式的检查点"上（v10@6500，今夜同协议实测 MW 302.5）。

为什么必须补这一步（§5.55）
--------------------------
阶段 1 把机制挂在了 **v10@5500（健康点，MW 397.5 已在带内）** 上，只能看到"过度变大"
（S1 恒温 447、S2 下限 493，落区率反而从 0.77 掉到 0.20–0.29）。
这些机制的**目标场景**是"检查点掉坑"——所以要在坑里测：v10@6500（今夜同协议 MW **302.5**、HA 22.0）。

对照臂**无需重跑**：`outputs/night_pipeline_v10/v10_matched_sel/6500`（今夜同协议同种子）。
处理臂：S1 恒温 / S2 硬下限 24 / S3 硬尺寸窗(24, 32)，各 2 种子（2024/2025）= **6 次采样 ≈ 30 分钟**。
判据沿用预登记（§5.41）：尺寸落区（MW 中位 ≥340、HA ≥25）+ strict ≥ 官方同协议臂 90% + 生成数 ≥n。

用法：python outputs/night_pipeline_v10/run_dip_test.py [--wait-pid <阶段1 PID>]
"""
import argparse
import os
import subprocess
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PY = sys.executable
BASE = os.path.join(ROOT, "outputs", "dip_test")
CKPT = os.path.join(ROOT, "logs", "train_gpcr_mass_v10_2026_09_28__02_29_15",
                    "checkpoints", "6500.pt")
SEEDS = (2024, 2025)
ARMS = {
    "off": [],                                                  # 对照（已有，可重跑做校验）
    "s1on": ["--schedule-target-ha", "28", "--schedule-slack", "4",
             "--schedule-relax", "0.35", "--schedule-warmup", "0.25"],
    "s2floor": ["--schedule-floor-steps", "24", "--schedule-floor-thr", "-2.0"],
    "s3window": ["--schedule-floor-steps", "24", "--schedule-floor-thr", "-2.0",
                 "--schedule-ceiling-ha", "32"],
}


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with open(os.path.join(BASE, "dip_test.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def pid_alive(pid):
    out = subprocess.run(["tasklist", "/FI", "PID eq %d" % pid],
                         capture_output=True, text=True, errors="ignore").stdout
    return str(pid) in out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wait-pid", type=int, default=0)
    ap.add_argument("--max-wait-h", type=float, default=2.0)
    ap.add_argument("--n", type=int, default=60)
    args = ap.parse_args()

    os.makedirs(BASE, exist_ok=True)
    if args.wait_pid:
        t0 = time.time()
        log("等待前置实验（PID %d）结束…" % args.wait_pid)
        while pid_alive(args.wait_pid):
            if time.time() - t0 > args.max_wait_h * 3600:
                log("等待超时，放弃")
                return 1
            time.sleep(30)
        log("前置实验已结束（等待 %.0f 分钟）" % ((time.time() - t0) / 60.0))
        time.sleep(15)

    log("凹坑测试开始：v10@6500（对照=今夜同协议 MW 302.5）| 臂=%s | 种子=%s"
        % (",".join(ARMS), SEEDS))
    for arm, extra in ARMS.items():
        for seed in SEEDS:
            outdir = os.path.join(BASE, "%s_s%d" % (arm, seed))
            cmd = [PY, os.path.join(ROOT, "design.py"), "--target", "A2A",
                   "--num-samples", str(args.n), "--beam", "50", "--max-steps", "40",
                   "--seed", str(seed), "--ckpt", CKPT,
                   "--frontier-threshold", "0.0", "--outdir", outdir] + extra
            with open(os.path.join(BASE, "run_%s_s%d.log" % (arm, seed)), "w",
                      encoding="utf-8") as fh:
                t0 = time.time()
                rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
            log("  %-9s seed %d | rc=%d | %.1f 分钟 | %s"
                % (arm, seed, rc, (time.time() - t0) / 60.0, os.path.relpath(outdir, ROOT)))
    log("凹坑测试采样结束（分析见 analyze_dip_test.py）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
