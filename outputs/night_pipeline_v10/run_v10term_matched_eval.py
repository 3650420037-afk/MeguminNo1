# -*- coding: utf-8 -*-
"""v10term **同迭代匹配**扫描（替代等距铺开）：直接对 v10 的已知尺寸做配对比较。

为什么换设计（2026-09-29 19:55 决定）
------------------------------------
原计划"均匀铺开 8 个检查点"有两个问题：
1. **早期检查点极慢**：250.pt 这类未训练好的检查点分子**不终止**，会跑满 `max_steps=40`，
   实测 ~8–14 分钟/点（成熟点约 6 分钟）；9 个点要 ~1.5–2 小时，而信息量低。
2. **不能配对**：均匀铺开拿到的是 4750/5750/7000/8000，而 v10 的已知尺寸在
   **5500→MW 393–407 / 6500→302 / 7500→359**（§5.33）—— 迭代不对齐就只能"看趋势"，不能配对比较。

改为扫 **5500 / 6500 / 7500（与 v10 同迭代）+ 8000（本轮终点）+ best（val 最优，仅作对照）**：
5 个点 ≈ 40 分钟，且能直接回答"**pos_weight=2.0 是否消除了 6500 的尺寸凹坑**"。

判据沿用预登记（§5.41）：尺寸（MW≥340 且 HA≥25）+ strict ≥ 官方基线 90% + 生成数 ≥n。
"""
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
RUN = os.path.join(ROOT, "logs", "train_gpcr_mass_v10term_2026_09_29__18_30_02", "checkpoints")
NIGHT = os.path.join(ROOT, "outputs", "night_pipeline_v10")
CKPTS = [os.path.join(RUN, "%d.pt" % it) for it in (5500, 6500, 7500, 8000)] + \
        [os.path.join(RUN, "best.pt")]


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with open(os.path.join(NIGHT, "eval_v10term.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def main():
    for p in CKPTS:
        if not os.path.exists(p):
            print("缺少检查点: %s" % p)
            return 1
    log("v10term 同迭代匹配扫描开始（%d 个点，与 v10 的 5500/6500/7500 可直接配对）"
        % len(CKPTS))
    log("点: %s" % ", ".join(os.path.basename(p) for p in CKPTS))
    cmd = [PY, os.path.join(ROOT, "src", "scripts", "select_ckpt_by_generation.py"),
           "--ckpts"] + CKPTS + [
        "--target", "A2A", "--num-samples", "60", "--beam", "50", "--max-steps", "40",
        "--seed", "2024", "--frontier-threshold", "0.0", "--rank-by", "both",
        "--keep-dir", os.path.join(NIGHT, "v10term_sel"),
        "--out", os.path.join(NIGHT, "sel_v10term.csv")]
    t0 = time.time()
    with open(os.path.join(NIGHT, "eval_v10term_stdout.log"), "a", encoding="utf-8") as fh:
        fh.write("\n===== 同迭代匹配扫描 %s =====\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
        rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
    log("扫描结束：rc=%d，用时 %.1f 分钟" % (rc, (time.time() - t0) / 60.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
