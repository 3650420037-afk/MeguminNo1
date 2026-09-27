# -*- coding: utf-8 -*-
"""清理采样过程产生的巨型中间快照, 释放磁盘。

背景(实测, 2026-09-27)
----------------------
`src/sample_for_pdb.py` 原先把**整个束状态**(含蛋白特征)在每一步都序列化到
`samples_<N>.pt`, 单个 170-210 MB。一次只生成 60 个分子的采样就能写出数 GB;
仓库 `outputs/` 下累积了 **412 个此类快照 = 43.26 GB**, 加上 52 个
`samples_all.pt` 共 24.46 GB, 使可用空间从 77.7 GB 掉到 30.5 GB ——
**正在跑的 6 小时对接有爆盘风险**(此前已发生过一次 101 GB 被快照撑满的事故)。

代码侧已修(改为 `--save-snapshots` 显式开关, 默认关闭, 与官方 `sample.py` 行为对齐)。
本脚本负责清理**历史遗留**的快照。

安全措施
--------
- **只删** `samples_init.pt` 与 `samples_<数字>.pt` 两类; 任何下游脚本都不读它们
  (下游只读 `SMILES.txt` / `SDF/`; `samples_all.pt` 作为完成标记被
  `src/sample.py` / `e2e_full_run.py` 读取, **本脚本默认不动它**)。
- 默认跳过最近 `--min-age` 秒内被修改的文件 —— 正在写入的采样进程不能被删。
- `--dry-run` 只统计不删除。

用法
----
    python src/scripts/purge_disk_hogs.py --dry-run
    python src/scripts/purge_disk_hogs.py
    python src/scripts/purge_disk_hogs.py --include-pools --pool-prefix ab_ sel_
"""
import argparse
import glob
import os
import re
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, OUTPUTS  # noqa: E402

SNAP_RE = re.compile(r"^samples_(\d+|init)\.pt$")


def human(n):
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024 or u == "GB":
            return "%.2f %s" % (n, u)
        n /= 1024.0


def main():
    ap = argparse.ArgumentParser(description="清理采样巨型中间快照")
    ap.add_argument("--root", default=OUTPUTS)
    ap.add_argument("--min-age", type=float, default=600.0,
                    help="跳过最近 N 秒内被修改的文件(默认 600, 保护正在写入的采样)")
    ap.add_argument("--include-pools", action="store_true",
                    help="**额外**删除 samples_all.pt(完成标记). 默认不删 —— "
                         "只有确认这些目录不再需要作为完成标记/位姿来源时才用")
    ap.add_argument("--pool-prefix", nargs="*", default=None,
                    help="配合 --include-pools: 只删这些前缀目录下的 samples_all.pt")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = args.root if os.path.isabs(args.root) else os.path.join(ROOT, args.root)
    snaps, pools, fresh = [], [], []
    total = 0

    for dp, _dn, fn in os.walk(root):
        for f in fn:
            p = os.path.join(dp, f)
            try:
                sz = os.path.getsize(p)
                age = time.time() - os.path.getmtime(p)
            except OSError:
                continue
            if SNAP_RE.match(f):
                if age < args.min_age:
                    fresh.append((p, age))
                else:
                    snaps.append((p, sz))
                    total += sz
            elif f == "samples_all.pt" and args.include_pools:
                if args.pool_prefix:
                    rel = os.path.relpath(p, root)
                    if not any(rel.startswith(x) for x in args.pool_prefix):
                        continue
                if age < args.min_age:
                    fresh.append((p, age))
                else:
                    pools.append((p, sz))
                    total += sz

    print("=" * 72)
    print("巨型中间快照清理  根目录: %s" % os.path.relpath(root, ROOT))
    print("=" * 72)
    print("  逐步快照 samples_{init,<N>}.pt : %4d 个, %s" % (len(snaps), human(sum(s for _, s in snaps))))
    if args.include_pools:
        print("  束状态 samples_all.pt         : %4d 个, %s" % (len(pools), human(sum(s for _, s in pools))))
    else:
        print("  束状态 samples_all.pt         : 保留(完成标记, 需 --include-pools 才删)")
    print("  正在写入而跳过                : %4d 个" % len(fresh))
    print("-" * 72)
    print("  可释放合计: %s" % human(total))

    if args.dry_run:
        print("(dry-run, 未删除)")
        return 0

    n = 0
    for p, _s in snaps + pools:
        try:
            os.remove(p)
            n += 1
        except OSError as e:
            print("  [失败] %s: %s" % (p, e))
    print("已删除 %d 个文件, 释放 %s" % (n, human(total)))

    # 顺手清掉因此变空的目录, 让目录结构保持可读
    removed_dirs = 0
    for dp, dn, fn in os.walk(root, topdown=False):
        if dp == root:
            continue
        try:
            if not os.listdir(dp):
                os.rmdir(dp)
                removed_dirs += 1
        except OSError:
            pass
    print("清理空目录 %d 个" % removed_dirs)
    return 2 if fresh else 0


if __name__ == "__main__":
    raise SystemExit(main())
