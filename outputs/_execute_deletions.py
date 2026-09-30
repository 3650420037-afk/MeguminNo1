# -*- coding: utf-8 -*-
"""按用户确认执行删除：组 1+2+4+5+8（**先干跑核对，再真删**）。

组 1: other/outputs_history/**/samples_<N>.pt        （逐步束状态快照，零读取者）
组 2: other/outputs_history/**/samples_all.pt        （历史归档的束状态）
组 4: data/*_processed.lmdb (+ -lock)                （数据集缓存，稀疏、占 0.25 GB、可自动重建）
组 5: B2AR 药库/中间产物/lmdb（**保留证据 CSV 与日志**）
组 8: outputs/mass_dock, outputs/mass_dock_build     （对接中间产物，可重跑）

**不删**：outputs/**/samples_all.pt（组 3）、results/official_weight_library 与 other/library_backups 其余（组 6）、
logs/*/checkpoints（组 7）。

用法：
    python outputs/_execute_deletions.py --dry-run     # 只列清单与体积
    python outputs/_execute_deletions.py --execute     # 真删（并写删除清单）
"""
import argparse
import glob
import io
import os
import shutil
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(ROOT, "other", "deletion_manifests",
                        "%s_delete_1_2_4_5_8.md" % time.strftime("%Y%m%d_%H%M%S"))
LOG = []


def is_file(p):
    return os.path.isfile(p)


def is_dir(p):
    return os.path.isdir(p)


def size_files(paths):
    tot = n = 0
    for p in paths:
        try:
            tot += os.path.getsize(p); n += 1
        except OSError:
            pass
    return tot, n


def size_any(paths):
    tot = n = 0
    for p in paths:
        if is_file(p):
            try:
                tot += os.path.getsize(p); n += 1
            except OSError:
                pass
        elif is_dir(p):
            for dp, _dn, fns in os.walk(p):
                for f in fns:
                    try:
                        tot += os.path.getsize(os.path.join(dp, f)); n += 1
                    except OSError:
                        pass
    return tot, n


def collect():
    groups = {}
    # 组 1/2
    groups["1. outputs_history 逐步快照 samples_<N>.pt"] = (
        "files", sorted(glob.glob(os.path.join(ROOT, "other", "outputs_history", "**",
                                              "samples_*.pt"), recursive=True)))
    # 上面 glob 会同时匹配 samples_all.pt，稍后按文件名分流
    g12 = groups.pop("1. outputs_history 逐步快照 samples_<N>.pt")
    step = [p for p in g12[1] if os.path.basename(p) != "samples_all.pt"]
    allp = [p for p in g12[1] if os.path.basename(p) == "samples_all.pt"]
    groups["1. outputs_history 逐步快照 samples_<N>.pt"] = ("files", step)
    groups["2. outputs_history 归档束状态 samples_all.pt"] = ("files", allp)

    # 组 4
    g4 = sorted(glob.glob(os.path.join(ROOT, "data", "*_processed.lmdb"))) + \
        sorted(glob.glob(os.path.join(ROOT, "data", "*_processed.lmdb-lock")))
    groups["4. data/*_processed.lmdb（数据集缓存，可重建）"] = ("files", g4)

    # 组 5（**只删库/中间产物/lmdb，保留证据 CSV/日志**）
    b2ar = [os.path.join(ROOT, "results", "B2AR"),
            os.path.join(ROOT, "data", "gpcr_b2ar"),
            os.path.join(ROOT, "data", "gpcr_b2ar_name2id.pt"),
            os.path.join(ROOT, "data", "gpcr_b2ar_processed.lmdb"),
            os.path.join(ROOT, "data", "gpcr_b2ar_processed.lmdb-lock"),
            os.path.join(ROOT, "outputs", "ab_dock_v3_B2AR")]
    b2ar += sorted(glob.glob(os.path.join(ROOT, "other", "library_backups", "*", "B2AR")))
    b2ar += sorted(glob.glob(os.path.join(ROOT, "other", "outputs_history", "*b2ar_library")))
    b2ar += sorted(glob.glob(os.path.join(ROOT, "other", "outputs_history", "*B2AR*")))
    b2ar += sorted(glob.glob(os.path.join(ROOT, "outputs", "*b2ar*")))
    b2ar += sorted(glob.glob(os.path.join(ROOT, "outputs", "*B2AR*")))
    # 证据保护：剔除 .csv/.log/.md（回归数字要留）
    keep_ext = (".csv", ".log", ".md", ".txt")
    b2ar = [p for p in set(b2ar) if not (is_file(p) and p.endswith(keep_ext))]
    groups["5. B2AR 药库/中间产物/lmdb（证据 CSV 保留）"] = ("any", sorted(b2ar))

    # 组 8
    groups["8. outputs/mass_dock + mass_dock_build（对接中间）"] = ("any", [
        os.path.join(ROOT, "outputs", "mass_dock"),
        os.path.join(ROOT, "outputs", "mass_dock_build")])
    return groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--execute", action="store_true")
    args = ap.parse_args()
    if not (args.dry_run or args.execute):
        print("请显式指定 --dry-run 或 --execute")
        return 2

    groups = collect()
    grand_b = grand_n = 0
    lines = ["# 删除清单（组 1/2/4/5/8，用户 2026-09-30 确认）", "",
             "> 生成时间：%s；模式：%s" % (time.strftime("%Y-%m-%d %H:%M:%S"),
                                        "EXECUTE" if args.execute else "DRY-RUN"), ""]
    for title, (kind, paths) in groups.items():
        paths = [p for p in paths if os.path.exists(p)]
        b, n = size_any(paths) if kind == "any" else size_files(paths)
        grand_b += b; grand_n += n
        lines += ["## %s" % title, "", "- 条目 **%d**，文件 **%d**，逻辑大小 **%.2f GB**" % (len(paths), n, b / 1e9)]
        for p in paths[:8]:
            if is_dir(p):
                sb, sn = size_any([p])
                lines.append("  - `%s/`（%d 文件，%.2f GB）" % (os.path.relpath(p, ROOT), sn, sb / 1e9))
            else:
                lines.append("  - `%s`" % os.path.relpath(p, ROOT))
        if len(paths) > 8:
            lines.append("  - …另有 %d 项" % (len(paths) - 8))
        lines.append("")
        print("%-52s 条目 %5d 文件 %7d %8.2f GB" % (title, len(paths), n, b / 1e9))
        LOG.append((title, paths))

    lines += ["## 合计", "", "- 文件 **%d**，逻辑大小 **%.2f GB**" % (grand_n, grand_b / 1e9), ""]
    print("\n合计：文件 %d，逻辑 %.2f GB" % (grand_n, grand_b / 1e9))

    if args.dry_run:
        print("\n[DRY-RUN] 未删除任何文件。")
    else:
        removed_b = removed_n = 0
        for title, paths in LOG:
            for p in paths:
                if not os.path.exists(p):
                    continue
                if is_dir(p):
                    b, n = size_any([p])
                    shutil.rmtree(p, ignore_errors=True)
                else:
                    try:
                        b, n = os.path.getsize(p), 1
                        os.remove(p)
                    except OSError:
                        b = n = 0
                removed_b += b; removed_n += n
            print("[已删] %s" % title)
        lines += ["## 执行结果", "",
                  "- 实际删除：文件 **%d**，逻辑 **%.2f GB**（%s）" % (removed_n, removed_b / 1e9,
                                                                    time.strftime("%H:%M:%S")), ""]
        print("\n实际删除：文件 %d，逻辑 %.2f GB" % (removed_n, removed_b / 1e9))

    os.makedirs(os.path.dirname(MANIFEST), exist_ok=True)
    io.open(MANIFEST, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("删除清单 → %s" % os.path.relpath(MANIFEST, ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
