# -*- coding: utf-8 -*-
"""补跑药库生成缺失的第 3 步：**对接 + 用对接分重新排序**。

背景（事故留痕）：`predict.py --library` 只做"打分排序"，**不做对接** —— 对接由
`src/scripts/dock_library.py` 单独完成。夜里的药库生成漏了这一步，导致新库 results.csv 的
"关键预测指标_对接分Vina(kcal/mol)" 列**全为空**、排序也只按类药性。本脚本补齐：

    对每个靶点：dock_library.py --library <T>/library --target <T>
              → predict.py --library <T>/library --results-dir <T>/results   （带对接分重排）

用法：python outputs/night_pipeline_v10/finish_library_docking.py [--targets ...] [--parallel 6]
"""
import argparse
import io
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
STAGING = os.path.join(ROOT, "outputs", "library_v10")


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with io.open(os.path.join(STAGING, "dock_finish.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", default=["A2A", "B2AR", "D3", "5HT2B"])
    ap.add_argument("--parallel", type=int, default=6)
    args = ap.parse_args()
    os.makedirs(STAGING, exist_ok=True)
    for t in args.targets:
        lib = os.path.join(STAGING, t, "library")
        if not os.path.exists(os.path.join(lib, "compounds.csv")):
            log("[%s] 缺少 library/compounds.csv，跳过" % t)
            continue
        log("[%s] 对接中（dock_library.py, parallel=%d）..." % (t, args.parallel))
        t0 = time.time()
        with io.open(os.path.join(STAGING, t, "04_dock.log"), "w", encoding="utf-8") as fh:
            rc = subprocess.run([PY, os.path.join(ROOT, "src", "scripts", "dock_library.py"),
                                 "--library", lib, "--target", t, "--parallel", str(args.parallel)],
                                cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
        log("[%s] 对接 rc=%d，用时 %.1f 分钟" % (t, rc, (time.time() - t0) / 60.0))
        if rc != 0:
            continue
        log("[%s] 用对接分重新排序（predict.py --library）..." % t)
        t0 = time.time()
        with io.open(os.path.join(STAGING, t, "05_predict.log"), "w", encoding="utf-8") as fh:
            rc2 = subprocess.run([PY, os.path.join(ROOT, "predict.py"), "--target", t,
                                  "--library", lib, "--top", "300", "--seed", "2024",
                                  "--results-dir", os.path.join(STAGING, t, "results"),
                                  "--outdir", os.path.join(STAGING, t, "predict")],
                                 cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
        log("[%s] 重排 rc=%d，用时 %.1f 分钟" % (t, rc2, (time.time() - t0) / 60.0))
    log("对接补齐完成")
    return 0


if __name__ == "__main__":
    sys.exit(main())
