# -*- coding: utf-8 -*-
"""全仓库磁盘盘点（**只读，绝不删除**）：单遍 os.scandir 遍历，产出可决策的清单。

为什么不用 du：Git Bash 的 du 在 Windows 上每个 stat 都要转换，且我上一版让它走了两遍
（`du -sh .` 之后又 `du -sh */`）。本脚本只走**一遍**，同时统计：
- 顶层目录体积/文件数；logs、outputs、data、other、results 的**二级目录**体积排名
- 最大文件 Top 40；按扩展名分组的体积（.pt / .sdf / .csv / .json / .tfevents / .log）
- 空目录、0 字节文件、__pycache__ 数量
输出：outputs/_disk_inventory.md（同时把进度写 outputs/_disk_inventory.progress）

用法（后台）：python outputs/_disk_inventory.py
"""
import io
import os
import sys
import time
from collections import defaultdict

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs", "_disk_inventory.md")
PROG = os.path.join(ROOT, "outputs", "_disk_inventory.progress")
SKIP_DIRS = {".git", ".venv", "env", "node_modules", "__pycache__", ".pytest_cache"}
MAJOR = ("logs", "outputs", "data", "other", "results", "models", "docs", "src",
         "configs", "notebooks", "attachments", "dist", "build", "tools", "tests")

GROUP_PATTERNS = {
    ".pt (权重/检查点)": (".pt",),
    ".sdf (分子结构)": (".sdf",),
    ".csv (数据/证据)": (".csv",),
    ".json (数据/配置)": (".json",),
    ".tfevents (训练曲线)": (".tfevents",),
    ".log/.txt (日志)": (".log", ".txt"),
    ".pdb/.pdbqt (结构/对接)": (".pdb", ".pdbqt"),
    "其他": (),
}


def group_of(name):
    low = name.lower()
    for tag, exts in GROUP_PATTERNS.items():
        for e in exts:
            if low.endswith(e):
                return tag
    return "其他"


def main():
    t0 = time.time()
    top_size = defaultdict(int)
    top_files = defaultdict(int)
    second_size = {d: defaultdict(int) for d in ("logs", "outputs", "data", "other", "results")}
    second_files = {d: defaultdict(int) for d in second_size}
    groups = defaultdict(int)
    groups_n = defaultdict(int)
    big = []
    empty_dirs = []
    zero_files = []
    pycache = 0
    total = 0
    nfiles = 0
    scanned = 0

    def log_progress(msg):
        with io.open(PROG, "a", encoding="utf-8") as fh:
            fh.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), msg))

    log_progress("开始单遍遍历 %s" % ROOT)
    for dirpath, dirnames, filenames in os.walk(ROOT, topdown=True):
        # 跳过明显无关的大目录
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        if "__pycache__" in os.listdir(dirpath) if False else False:
            pass
        rel = os.path.relpath(dirpath, ROOT)
        parts = [] if rel == "." else rel.split(os.sep)
        tl = parts[0] if parts else "(根)"
        sl = parts[1] if len(parts) > 1 else "(直接子文件)"
        if not filenames and not dirnames and rel != ".":
            empty_dirs.append(rel)
        for fn in filenames:
            fp = os.path.join(dirpath, fn)
            try:
                sz = os.path.getsize(fp)
            except OSError:
                continue
            total += sz
            nfiles += 1
            scanned += 1
            top_size[tl] += sz
            top_files[tl] += 1
            if tl in second_size:
                second_size[tl][sl] += sz
                second_files[tl][sl] += 1
            g = group_of(fn)
            groups[g] += sz
            groups_n[g] += 1
            if sz == 0:
                zero_files.append(fp if len(zero_files) < 200 else zero_files[0])
            if sz > 20 * 1024 * 1024:
                big.append((sz, os.path.relpath(fp, ROOT)))
            if fn == "__pycache__" or fn.endswith(".pyc"):
                pycache += 1
        if scanned // 20000 > (scanned - len(filenames)) // 20000:
            log_progress("已扫描 %d 个文件，累计 %.1f GB" % (nfiles, total / 1e9))

    big.sort(reverse=True)
    lines = ["# 磁盘盘点（只读，未删除任何文件）", "",
             "- 根目录：`%s`" % ROOT,
             "- **总计：%.1f GB / %d 个文件**（耗时 %.1f 分钟）" % (total / 1e9, nfiles, (time.time() - t0) / 60),
             "", "## 一、顶层目录体积排名", "",
             "| 目录 | 体积 | 文件数 |", "| --- | --- | --- |"]
    for k, v in sorted(top_size.items(), key=lambda kv: -kv[1]):
        lines.append("| %s | %.2f GB | %d |" % (k, v / 1e9, top_files[k]))

    for d in ("logs", "outputs", "data", "other", "results"):
        if not second_size.get(d):
            continue
        lines += ["", "## 二、`%s/` 下的二级目录 Top 25" % d, "",
                  "| 子目录 | 体积 | 文件数 |", "| --- | --- | --- |"]
        for k, v in sorted(second_size[d].items(), key=lambda kv: -kv[1])[:25]:
            lines.append("| %s | %.2f GB | %d |" % (k, v / 1e9, second_files[d][k]))

    lines += ["", "## 三、按类型分组", "", "| 类型 | 体积 | 文件数 |", "| --- | --- | --- |"]
    for k, v in sorted(groups.items(), key=lambda kv: -kv[1]):
        lines.append("| %s | %.2f GB | %d |" % (k, v / 1e9, groups_n[k]))

    lines += ["", "## 四、最大文件 Top 40（>20 MB）", "", "| 体积 | 文件 |", "| --- | --- |"]
    for sz, p in big[:40]:
        lines.append("| %.1f MB | `%s` |" % (sz / 1e6, p))

    lines += ["", "## 五、空目录 / 0 字节文件 / 缓存",
              "- 空目录：**%d** 个" % len(empty_dirs),
              "- 0 字节文件：**%d** 个（样本：%s）" % (len(zero_files), ", ".join(zero_files[:5])),
              "- `__pycache__`/`.pyc`：**%d** 个" % pycache, ""]
    if empty_dirs:
        lines += ["空目录清单（前 40）："] + ["- `%s`" % d for d in empty_dirs[:40]]

    io.open(OUT, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    log_progress("完成：%.1f GB / %d 文件 → %s" % (total / 1e9, nfiles, OUT))
    print("完成：%.1f GB / %d 文件；清单 → %s" % (total / 1e9, nfiles, OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
