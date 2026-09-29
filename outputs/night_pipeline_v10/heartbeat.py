# -*- coding: utf-8 -*-
"""夜间心跳：每 30 分钟往 `heartbeat.log` 追加一条状态快照（用户要求的 30 分钟一次心跳）。

快照内容（全部只读，不干扰任何实验）：
- 时间戳；正在跑的 GPU 关键任务与其进度（凹坑测试 / 库生成 / 训练）
- 当前默认权重与 git 状态（本地领先 origin 多少提交，按用户要求**不推送**）
- 最近一次实测读数（若相关产物已落盘）
- 下一步计划（从固定的夜班路线图里取）

用法（后台常驻，不占 GPU）：
    python outputs/night_pipeline_v10/heartbeat.py --interval 1800
"""
import argparse
import glob
import io
import os
import re
import subprocess
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NIGHT = os.path.join(ROOT, "outputs", "night_pipeline_v10")
LOG = os.path.join(NIGHT, "heartbeat.log")
sys.path.insert(0, os.path.join(ROOT, "src"))


def sh(cmd):
    try:
        return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                              errors="ignore").stdout.strip()
    except Exception as e:
        return "(失败: %s)" % e


def procs():
    out = sh(["powershell", "-NoProfile", "-Command",
              "Get-Process python -ErrorAction SilentlyContinue | "
              "Select-Object -ExpandProperty Id"])
    return [x for x in out.split() if x.strip().isdigit()]


def tail(path, n=1):
    if not os.path.exists(path):
        return "(无)"
    lines = [l for l in io.open(path, encoding="utf-8", errors="ignore").read().splitlines() if l.strip()]
    return lines[-n] if lines else "(空)"


def dip_progress():
    d = os.path.join(ROOT, "outputs", "dip_test")
    arms = sorted(os.path.basename(x.rstrip("/")) for x in glob.glob(os.path.join(d, "*_s*")))
    return "%d 臂目录：%s" % (len(arms), ", ".join(arms) if arms else "（尚无）")


def lib_progress():
    d = os.path.join(ROOT, "outputs", "library_v10")
    runs = glob.glob(os.path.join(d, "**", "SMILES.txt"), recursive=True)
    n = 0
    for r in runs:
        n += sum(1 for _ in io.open(r, encoding="utf-8", errors="ignore"))
    return "%d 个采样会话，累计 %d 个分子" % (len(runs), n)


def snapshot():
    gpu_jobs = procs()
    lines = [
        "=" * 78,
        "[心跳] %s" % time.strftime("%Y-%m-%d %H:%M:%S"),
        "- 默认权重: models/7-eonmol_ft_gpcr_v10.pt（v10@5500；阈值侧车 v10→0.0、v3→−0.5）",
        "- python 进程: %s" % (", ".join(gpu_jobs) if gpu_jobs else "（无）"),
        "- 凹坑测试: %s" % dip_progress(),
        "- 凹坑测试日志末行: %s" % tail(os.path.join(ROOT, "outputs", "dip_test", "dip_test.log")),
        "- 新药库生成: %s" % lib_progress(),
        "- git: 本地领先 origin %s 个提交（**按用户要求不推送**）"
        % sh(["git", "rev-list", "--count", "origin/main..HEAD"]),
        "- 最近提交: %s" % sh(["git", "log", "--oneline", "-1"]),
        "- 夜班路线: ①凹坑测试判定 → ②S3 调窗复验(22/30) → ③(若都无效)V2 逐样本平衡训练 → "
        "④用 v10 重出四个药库（备份旧库后替换）→ ⑤更新 附件3.docx 区赛二–六部分 → ⑥同步文档",
        "- 待用户决定: 是否推送 GitHub；生成数硬规则是否修订；A′ 是否全量重验",
    ]
    with io.open(LOG, "a", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interval", type=int, default=1800, help="心跳间隔（秒，默认 1800 = 30 分钟）")
    ap.add_argument("--once", action="store_true", help="只打一次快照（调试用）")
    args = ap.parse_args()
    os.makedirs(NIGHT, exist_ok=True)
    snapshot()
    if args.once:
        return 0
    while True:
        time.sleep(args.interval)
        try:
            snapshot()
        except Exception as e:
            with io.open(LOG, "a", encoding="utf-8") as fh:
                fh.write("[心跳] %s 快照失败: %s\n" % (time.strftime("%H:%M:%S"), e))
    return 0


if __name__ == "__main__":
    sys.exit(main())
