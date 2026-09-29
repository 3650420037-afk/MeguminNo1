# -*- coding: utf-8 -*-
"""v10（pos_weight=1.0）的**同协议配对参照臂**：等 v10term 扫描结束后自动跑。

为什么需要（否则又是一次协议不一致的对比）
------------------------------------------
§5.33 记录的 v10 尺寸（5500→MW 393–407 / 6500→302 / 7500→359）来自**别的协议**（n=50/steps60 等）。
要回答"`pos_weight=2.0` 是否消除了 6500 的尺寸凹坑"，必须让 v10 与 v10term
在**同靶点、同 n、同 beam、同 steps、同种子、同阈值**下采样 —— 即本脚本。

做法：等 PID（v10term 扫描）退出 → 用**完全相同的参数**扫 v10 的 5500/6500/7500/8000，
分子保留到 `outputs/night_pipeline_v10/v10_matched_sel/`，CSV 写 `sel_v10_matched.csv`。

用法：python outputs/night_pipeline_v10/run_v10_matched_reference.py --wait-pid 28504
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
NIGHT = os.path.join(ROOT, "outputs", "night_pipeline_v10")
RUN10 = os.path.join(ROOT, "logs", "train_gpcr_mass_v10_2026_09_28__02_29_15", "checkpoints")
ITERS = (5500, 6500, 7500, 8000)


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with open(os.path.join(NIGHT, "ref_v10_matched.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def pid_alive(pid):
    out = subprocess.run(["tasklist", "/FI", "PID eq %d" % pid],
                         capture_output=True, text=True, errors="ignore").stdout
    return str(pid) in out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wait-pid", type=int, default=0)
    ap.add_argument("--max-wait-h", type=float, default=3.0)
    args = ap.parse_args()

    if args.wait_pid:
        t0 = time.time()
        log("等待 v10term 扫描（PID %d）结束…" % args.wait_pid)
        while pid_alive(args.wait_pid):
            if time.time() - t0 > args.max_wait_h * 3600:
                log("等待超时，放弃")
                return 1
            time.sleep(30)
        log("前置扫描已结束（等待 %.0f 分钟），开始 v10 参照臂" % ((time.time() - t0) / 60.0))
        time.sleep(20)

    ckpts = [os.path.join(RUN10, "%d.pt" % it) for it in ITERS]
    miss = [p for p in ckpts if not os.path.exists(p)]
    if miss:
        log("缺少 v10 检查点: %s" % [os.path.basename(p) for p in miss])
        return 1
    cmd = [PY, os.path.join(ROOT, "src", "scripts", "select_ckpt_by_generation.py"),
           "--ckpts"] + ckpts + [
        "--target", "A2A", "--num-samples", "60", "--beam", "50", "--max-steps", "40",
        "--seed", "2024", "--frontier-threshold", "0.0", "--rank-by", "both",
        "--keep-dir", os.path.join(NIGHT, "v10_matched_sel"),
        "--out", os.path.join(NIGHT, "sel_v10_matched.csv")]
    log("v10 参照臂开始（%s，与 v10term 同协议）" % ", ".join(os.path.basename(p) for p in ckpts))
    t0 = time.time()
    with open(os.path.join(NIGHT, "ref_v10_matched_stdout.log"), "a", encoding="utf-8") as fh:
        fh.write("\n===== v10 参照臂 %s =====\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
        rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
    log("v10 参照臂结束：rc=%d，用时 %.1f 分钟" % (rc, (time.time() - t0) / 60.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
