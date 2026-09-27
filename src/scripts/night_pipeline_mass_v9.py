# -*- coding: utf-8 -*-
"""夜间流水线总控: 等全量对接完成 -> 建 gpcr_mass_v1 -> 校验 -> 训练 v9。

为什么需要它
------------
全量对接要跑数小时, 之后还 demands 依次执行"布局转换 -> 建数据集 -> 校验 -> 训练"。
人工盯着不现实, 而且中间任何一步失败若不检查就会把坏数据喂进训练。
本脚本把这条链固定下来, 每步都校验, 失败即停并留下明确原因; 且**可重复运行**:
   - 对接未完成 -> 直接退出(不误建半成品数据集);
   - 数据集已存在且校验通过 -> 跳过重建;
   - 训练日志目录已存在 -> 视为已启动, 不再重复启动。

用法
----
    python src/scripts/night_pipeline_mass_v9.py --stage all      # 全链
    python src/scripts/night_pipeline_mass_v9.py --stage build    # 只建数据集+校验
    python src/scripts/night_pipeline_mass_v9.py --stage train    # 只训练
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, DATA, LOGS, PYTHON  # noqa: E402

TARGETS = ["A2A", "B2AR", "D3", "5HT2B"]
MASS = os.path.join(ROOT, "outputs", "mass_dock")
BUILD = os.path.join(ROOT, "outputs", "mass_dock_build")
DS = os.path.join(DATA, "gpcr_mass_v1")
CFG = os.path.join(ROOT, "configs", "train_gpcr_mass_v9.yml")


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def run(cmd, tag, logdir):
    os.makedirs(logdir, exist_ok=True)
    lf = os.path.join(logdir, "%s.log" % tag)
    log("执行 %s -> %s" % (tag, os.path.relpath(lf, ROOT)))
    with open(lf, "w", encoding="utf-8", errors="ignore") as f:
        rc = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        log("!! %s 失败 (rc=%d), 详见 %s" % (tag, rc, os.path.relpath(lf, ROOT)))
        # 把日志尾部打出来, 便于无人值守时也能直接看到原因
        tail = open(lf, encoding="utf-8", errors="ignore").read().splitlines()[-15:]
        for l in tail:
            print("    | " + l)
        raise SystemExit(1)
    return rc


def docking_ready():
    """四个靶点的 summary_merged.csv 全部存在才算对接完成。"""
    miss = [T for T in TARGETS
            if not os.path.exists(os.path.join(MASS, T, "summary_merged.csv"))]
    if miss:
        log("对接未完成, 缺 %s 的 summary_merged.csv" % miss)
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description="夜间流水线: 建 mass_v1 数据集并训练 v9")
    ap.add_argument("--stage", choices=["all", "build", "train"], default="all")
    ap.add_argument("--wait-dock", action="store_true",
                    help="若对接未完成则轮询等待(默认直接退出, 避免误建半成品)")
    ap.add_argument("--poll", type=int, default=600, help="轮询间隔秒")
    args = ap.parse_args()

    logdir = os.path.join(ROOT, "outputs", "night_pipeline")

    if args.stage in ("all", "build"):
        while not docking_ready():
            if not args.wait_dock:
                log("加 --wait-dock 可轮询等待; 本次退出。")
                return 3
            time.sleep(args.poll)
        log("四个靶点对接均已完成")

        # 1) 布局转换 + 位姿完整性校验
        run([PYTHON, os.path.join(ROOT, "src", "scripts", "adapt_mass_dock.py"),
             "--src", MASS, "--dst", BUILD], "adapt_mass_dock", logdir)

        # 2) 建数据集(复用已完成的对接, 不重跑)
        if os.path.exists(os.path.join(DS, "split_by_name.pt")):
            log("gpcr_mass_v1 已存在, 跳过重建")
        else:
            run([PYTHON, os.path.join(ROOT, "src", "scripts", "build_docked_dataset.py"),
                 "--work", BUILD, "--reuse-dock", "--out", DS], "build_docked_dataset", logdir)

        # 3) 校验 —— 必须看**每口袋密度**(决定性指标)与一致性
        run([PYTHON, os.path.join(ROOT, "src", "scripts", "verify_gpcr_datasets.py"),
             "gpcr_mass_v1"], "verify_gpcr_datasets", logdir)
        # 把校验摘要打出来
        vf = os.path.join(logdir, "verify_gpcr_datasets.log")
        if os.path.exists(vf):
            for l in open(vf, encoding="utf-8", errors="ignore").read().splitlines():
                if any(k in l for k in ("index_pairs", "per_pocket_density", "consistency",
                                        "WARN", "pocket_files")):
                    log("  " + l.strip())

    if args.stage in ("all", "train"):
        if not os.path.exists(os.path.join(DS, "split_by_name.pt")):
            log("数据集不存在, 不能训练。先跑 --stage build")
            return 4
        existing = glob.glob(os.path.join(LOGS, "train_gpcr_mass_v9_*"))
        if existing:
            log("已有 v9 训练目录 %s, 不再重复启动" % os.path.basename(existing[-1]))
            return 0
        log("启动 v9 训练(长时间运行, 日志在 logs/train_gpcr_mass_v9_*)")
        run([PYTHON, os.path.join(ROOT, "train.py"), "--config", CFG, "--device", "cuda"],
            "train_gpcr_mass_v9", logdir)

    log("流水线结束")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
