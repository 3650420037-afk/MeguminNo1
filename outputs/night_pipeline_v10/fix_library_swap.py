# -*- coding: utf-8 -*-
"""补齐/校正药库替换：把 staging 的**库文件**（compounds.csv / sdf / docking_summary.csv / library.db /
rejected.csv / summary.txt / top_candidates.csv）一并搬到 results/ 对应位置。

背景（事故留痕）：`regenerate_libraries.py --swap` 只搬了 predict 的产物（results.csv / run_manifest.json /
structures），**漏了 library 目录** → 替换后 `results/compounds.csv` 仍是旧库的 521 个分子，
与新的 results.csv（300 候选）不一致。本脚本把两边对齐，并在搬动前把"当前 results/" 再备份一次。

用法：python outputs/night_pipeline_v10/fix_library_swap.py [--targets ...]
"""
import argparse
import io
import os
import shutil
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STAGE = os.path.join(ROOT, "outputs", "library_v10")
RESULTS = os.path.join(ROOT, "results")
LIB_FILES = ("compounds.csv", "docking_summary.csv", "library.db", "rejected.csv",
             "summary.txt", "top_candidates.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", default=["A2A", "B2AR", "D3", "5HT2B"])
    ap.add_argument("--no-backup", action="store_true")
    args = ap.parse_args()

    if not args.no_backup:
        dest = os.path.join(ROOT, "other", "library_backups",
                            "%s_after_swap" % time.strftime("%Y%m%d_%H%M%S"))
        for name in os.listdir(RESULTS):
            src = os.path.join(RESULTS, name)
            if name == "official_weight_library":
                continue
            dst = os.path.join(dest, name)
            if os.path.isdir(src):
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                os.makedirs(dest, exist_ok=True)
                shutil.copy2(src, dst)
        print("已把当前 results/ 再备份 → %s" % os.path.relpath(dest, ROOT))

    for t in args.targets:
        lib = os.path.join(STAGE, t, "library")
        if not os.path.exists(os.path.join(lib, "compounds.csv")):
            print("[%s] 缺 library/compounds.csv，跳过" % t)
            continue
        dst_dir = RESULTS if t == "A2A" else os.path.join(RESULTS, t)
        os.makedirs(dst_dir, exist_ok=True)
        moved = []
        for name in LIB_FILES:
            s = os.path.join(lib, name)
            if os.path.exists(s):
                shutil.copy2(s, os.path.join(dst_dir, name))
                moved.append(name)
        # sdf/ 目录
        s_sdf, d_sdf = os.path.join(lib, "sdf"), os.path.join(dst_dir, "sdf")
        if os.path.isdir(s_sdf):
            shutil.copytree(s_sdf, d_sdf, dirs_exist_ok=True)
            moved.append("sdf/")
        # predict 的阶段产物（results.csv / structures / run_manifest）
        pres = os.path.join(STAGE, t, "results")
        for name in ("results.csv", "run_manifest.json"):
            s = os.path.join(pres, name)
            if os.path.exists(s):
                shutil.copy2(s, os.path.join(dst_dir, name))
                moved.append(name)
        if t == "A2A":
            # A2A 的清单同时放根目录（历史结构如此）
            for name in ("results.csv", "run_manifest.json", "admet_screen.csv",
                         "top_candidates.csv"):
                s = os.path.join(pres, name) if "results" in name or "manifest" in name \
                    else os.path.join(lib, name)
                if os.path.exists(s):
                    shutil.copy2(s, os.path.join(RESULTS, name))
        print("[%s] 已同步库文件：%s → %s" % (t, ", ".join(moved), os.path.relpath(dst_dir, ROOT)))

    # 一致性校验
    print("\n一致性校验：")
    ok = True
    for t in args.targets:
        dst_dir = RESULTS if t == "A2A" else os.path.join(RESULTS, t)
        c = os.path.join(dst_dir, "compounds.csv")
        r = os.path.join(dst_dir, "results.csv")
        nc = sum(1 for _ in io.open(c, encoding="utf-8-sig")) - 1 if os.path.exists(c) else -1
        nr = sum(1 for _ in io.open(r, encoding="utf-8-sig")) - 1 if os.path.exists(r) else -1
        same = (nc >= nr >= 0)
        ok &= same
        print("  %-6s compounds=%d 候选=%d %s" % (t, nc, nr, "OK" if same else "**不一致**"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
