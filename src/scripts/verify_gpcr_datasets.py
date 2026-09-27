# -*- coding: utf-8 -*-
"""Verify built GPCR datasets: pair counts, splits, reports, file integrity."""
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import DATA, src_on_path
src_on_path()

DATASETS = [
    "gpcr_a2a",
    "gpcr_b2ar",
    "gpcr_d3",
    "gpcr_5ht2b",
    "gpcr_multitarget_v2",
]


def verify(name):
    path = os.path.join(DATA, name)
    print("=== %s ===" % name)
    with open(os.path.join(path, "index.pkl"), "rb") as handle:
        entries = pickle.load(handle)
    with open(os.path.join(path, "split_by_name.pt"), "rb") as handle:
        import torch
        split = torch.load(handle)
    train = split.get("train", [])
    test = split.get("test", [])
    print("index_pairs=%d train=%d test=%d" % (len(entries), len(train), len(test)))
    missing = 0
    for pocket_name, ligand_name in train + test:
        if not os.path.isfile(os.path.join(path, pocket_name)):
            missing += 1
        if not os.path.isfile(os.path.join(path, ligand_name)):
            missing += 1
    print("missing_files=%d" % missing)
    # index <-> split 一致性校验 (旧版只查文件存在性, 会漏掉 split/index 不匹配、
    # 重复或遗漏条目 —— 实测曾出现 build_report 与 index 数字矛盾却"通过校验")
    idx_set = {(e[0], e[1]) for e in entries}
    sp_set = set(train + test)
    dup = len(train) + len(test) - len(sp_set)
    not_in_index = len(sp_set - idx_set)
    not_in_split = len(idx_set - sp_set)
    print("split_dupes=%d split_not_in_index=%d index_not_in_split=%d" % (dup, not_in_index, not_in_split))
    if dup or not_in_index or not_in_split:
        print("  WARN: index 与 split 不一致 (split 由 rebuild_split 重切后可能未同步, 检查 scripts 流程)")
    print("consistency=%s" % ("OK" if not (dup or not_in_index or not_in_split) else "MISMATCH"))
    pockets = sorted({item[0] for item in entries})
    print("pocket_files=%s" % pockets)
    # 每口袋样本密度 —— 这是微调成败的**决定性**指标(实测 1.04 样本/口袋会让停止策略
    # 不可学、分子永不终止; 40-127/口袋则恢复正常)。旧版只报总对数, 会掩盖
    # "总量很大但每口袋只有 1 条"这种致命分布; 因此必须显式报出来。
    per_pocket = {}
    for item in entries:
        per_pocket[item[0]] = per_pocket.get(item[0], 0) + 1
    if per_pocket:
        vals = sorted(per_pocket.values())
        print("per_pocket_density: n_pockets=%d min=%d median=%d max=%d"
              % (len(vals), vals[0], vals[len(vals) // 2], vals[-1]))
        for pk in sorted(per_pocket, key=lambda k: -per_pocket[k])[:12]:
            print("    %-40s %5d" % (pk, per_pocket[pk]))
        if len(vals) > 1 and vals[0] <= 5:
            print("  WARN: 存在每口袋 <=5 个样本的口袋 (%d 个) —— 该分布曾导致停止策略不可学"
                  % sum(1 for v in vals if v <= 5))
    report_path = os.path.join(path, "build_report.txt")
    if os.path.isfile(report_path):
        with open(report_path, "r", encoding="utf-8") as handle:
            content = handle.read()
        summary_lines = [
            line for line in content.splitlines()
            if not line.startswith("skipped_")
        ]
        print("build_report:")
        for line in summary_lines:
            print("  %s" % line)
    else:
        print("build_report: MISSING")
    return 0


def main():
    # 允许显式指定数据集, 这样新建的数据集(如 gpcr_mass_v1)不必改源码就能校验 ——
    # 旧版把清单写死在源码里, 新数据集只能靠改代码才能验证, 容易漏检。
    args = sys.argv[1:]
    names = args if args else DATASETS
    failed = 0
    for name in names:
        try:
            verify(name)
        except Exception as error:
            failed += 1
            print("VERIFY_FAILED %s: %r" % (name, error))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
