# -*- coding: utf-8 -*-
"""v10term 训练结束后的自动下游：等训练进程退出 → 按检查点扫生成级指标。

纪律（预登记）：
- **不看 val loss 选点**（§5.33：val-loss-best 系统性偏早）；
- 用固定种子、同一阈值（v10 推荐 0.0），对**均匀铺开的检查点**扫生成指标；
- 分子保留（--keep-dir），便于后续算 strict 通过率与做对接 A/B。

用法：python outputs/night_pipeline_v10/run_v10term_eval.py [--wait-pid 28824]
"""
import argparse
import glob
import os
import re
import subprocess
import sys

# 无人值守运行时 stdout 常被重定向为 GBK 管道：非 GBK 字符（⚠/… 等）会直接崩脚本。
# 这是 A′ v1.0 采样器真实踩过的坑（当时两个臂全崩、CSV 0 行），此处统一加固。
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PY = sys.executable
OUTDIR = os.path.join(ROOT, "outputs", "night_pipeline_v10")


def log(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(os.path.join(OUTDIR, "eval_v10term.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def pid_alive(pid):
    out = subprocess.run(["tasklist", "/FI", "PID eq %d" % pid],
                         capture_output=True, text=True, errors="ignore").stdout
    return str(pid) in out


def find_run_dir():
    """取 checkpoints 最多的那个 v10term run 目录。"""
    best, best_n = None, -1
    for d in glob.glob(os.path.join(ROOT, "logs", "train_gpcr_mass_v10term_*")):
        ck = os.path.join(d, "checkpoints")
        n = len(glob.glob(os.path.join(ck, "[0-9]*.pt")))
        if n > best_n:
            best, best_n = d, n
    return best, best_n


def pick_ckpts(run_dir, k=8):
    ck_dir = os.path.join(run_dir, "checkpoints")
    iters = sorted(int(re.match(r"(\d+)\.pt$", os.path.basename(p)).group(1))
                   for p in glob.glob(os.path.join(ck_dir, "[0-9]*.pt")))
    if not iters:
        return []
    # 均匀铺开 k 个 + best.pt（best 是按 val loss 选的，只作为对照臂，不单独采信）
    if len(iters) <= k:
        chosen = iters
    else:
        step = (len(iters) - 1) / float(k - 1)
        chosen = sorted({iters[int(round(i * step))] for i in range(k)})
    paths = [os.path.join(ck_dir, "%d.pt" % i) for i in chosen]
    best = os.path.join(ck_dir, "best.pt")
    if os.path.exists(best):
        paths.append(best)
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wait-pid", type=int, default=0, help="先等这个 PID 退出(训练)")
    ap.add_argument("--max-wait-h", type=float, default=5.0)
    ap.add_argument("--num-samples", type=int, default=60)
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    if args.wait_pid:
        t0 = time.time()
        log("等待训练进程 PID %d 退出（上限 %.1f h）…" % (args.wait_pid, args.max_wait_h))
        while pid_alive(args.wait_pid):
            if time.time() - t0 > args.max_wait_h * 3600:
                log("等待超时，放弃本轮自动评估（训练可能仍在跑）")
                return 1
            time.sleep(60)
        log("训练进程已退出（等待 %.0f 分钟）" % ((time.time() - t0) / 60.0))
        time.sleep(30)  # 让最后一次落盘写完

    run_dir, n_ck = find_run_dir()
    log("run 目录 = %s（%d 个检查点）" % (run_dir, n_ck))
    ckpts = pick_ckpts(run_dir)
    if not ckpts:
        log("没有可用的检查点，退出")
        return 1
    log("将扫描 %d 个检查点：%s" % (len(ckpts), ", ".join(os.path.basename(p) for p in ckpts)))

    csv_out = os.path.join(OUTDIR, "sel_v10term.csv")
    keep = os.path.join(OUTDIR, "v10term_sel")
    cmd = [PY, os.path.join(ROOT, "src", "scripts", "select_ckpt_by_generation.py"),
           "--ckpts"] + ckpts + [
           "--target", "A2A", "--num-samples", str(args.num_samples),
           "--beam", "50", "--max-steps", "40", "--seed", "2024",
           "--frontier-threshold", "0.0",           # v10 的推荐阈值（随权重走）
           "--rank-by", "both", "--keep-dir", keep, "--out", csv_out]
    log("开始生成级扫描：%s" % " ".join(cmd[:6] + ["..."]))
    t0 = time.time()
    with open(os.path.join(OUTDIR, "eval_v10term_stdout.log"), "a", encoding="utf-8") as fh:
        fh.write("\n===== %s =====\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
        r = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    log("扫描结束：returncode=%d，用时 %.1f 分钟，CSV=%s" % (r.returncode, (time.time() - t0) / 60.0, csv_out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
