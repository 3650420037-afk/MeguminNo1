# -*- coding: utf-8 -*-
"""Verify built GPCR datasets: pair counts, splits, reports, file integrity."""
import os
import pickle
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

DATASETS = [
    "gpcr_a2a",
    "gpcr_b2ar",
    "gpcr_d3",
    "gpcr_5ht2b",
    "gpcr_multitarget_v2",
]


def verify(name):
    path = os.path.join(DATA_DIR, name)
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
    pockets = sorted({item[0] for item in entries})
    print("pocket_files=%s" % pockets)
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
    failed = 0
    for name in DATASETS:
        try:
            verify(name)
        except Exception as error:
            failed += 1
            print("VERIFY_FAILED %s: %r" % (name, error))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
