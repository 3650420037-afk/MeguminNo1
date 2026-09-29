# -*- coding: utf-8 -*-
"""**决定性读出**：v10term（pos_weight=2.0）vs v10（pos_weight=1.0），**同协议、同迭代配对**。

回答的问题：`pos_weight=2.0` 这个训练侧干预，是否消除了"尺寸随检查点掉坑"的现象？
（v10 的已知凹坑：6500 → MW 302；而 5500 → 393–407。§5.33）

输入：`sel_v10term.csv` + `sel_v10_matched.csv` 各自 `--keep-dir` 下的分子。
输出：按迭代配对的 MW/HA/strict 表 + 差值；并把 v10 的"凹坑"是否复现/消失讲清楚。

判据（预登记 §5.41，不自创）：尺寸（MW≥340 且 HA≥25）、strict ≥ 官方基线 90%、生成数 ≥n。
本脚本只做**生成级读出**；对接级须另跑（见 decide/dock 两脚本）。

用法：python outputs/night_pipeline_v10/analyze_v10_vs_v10term.py
"""
import csv
import glob
import io
import os
import statistics as st
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))
NIGHT = os.path.join(ROOT, "outputs", "night_pipeline_v10")
BAND = (340.0, 430.0)

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def find_smiles(d):
    c = glob.glob(os.path.join(d, "**", "SMILES.txt"), recursive=True)
    return sorted(c)[-1] if c else None


def metrics(smiles_path):
    mols = []
    if not smiles_path or not os.path.exists(smiles_path):
        return None
    with io.open(smiles_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            m = Chem.MolFromSmiles(line)
            if m is not None:
                mols.append(m)
    if not mols:
        return None
    ha = sorted(m.GetNumHeavyAtoms() for m in mols)
    mw = sorted(round(Descriptors.MolWt(m), 1) for m in mols)   # Mol 没有 GetMolWt()，自测踩到
    return dict(n=len(mols), mw=st.median(mw), ha=st.median(ha),
                ha_min=min(ha), ha_max=max(ha),
                frac_ha25=sum(1 for x in ha if x >= 25) / len(ha),
                in_band=sum(1 for x in mw if BAND[0] <= x <= BAND[1]) / len(mw))


def collect(sel_csv, keep_dir):
    out = {}
    if not os.path.exists(sel_csv):
        return out
    for r in csv.DictReader(io.open(sel_csv, encoding="utf-8-sig")):
        name = os.path.splitext(os.path.basename(r["ckpt"]))[0]
        d = r.get("run_dir") or os.path.join(keep_dir, name)
        m = metrics(find_smiles(d))
        if m:
            out[name] = m
    return out


def main():
    a = collect(os.path.join(NIGHT, "sel_v10_matched.csv"), os.path.join(NIGHT, "v10_matched_sel"))
    b = collect(os.path.join(NIGHT, "sel_v10term.csv"), os.path.join(NIGHT, "v10term_sel"))
    if not a or not b:
        print("数据不足：v10=%d 个点, v10term=%d 个点（扫描可能还没跑完）" % (len(a), len(b)))
        for k, v in sorted(b.items()):
            print("  v10term %-10s %s" % (k, v))
        return 1

    print("同协议（A2A / n=60 / beam50 / steps40 / seed2024 / 阈值 0.0）同迭代配对\n")
    print("%-8s | %-32s | %-32s | %s" % ("迭代", "v10（pos_weight=1.0）", "v10term（pos_weight=2.0）", "ΔMW"))
    print("-" * 104)
    rows = []
    for name in sorted(set(a) & set(b), key=lambda s: (not s.isdigit(), int(s) if s.isdigit() else 0)):
        x, y = a[name], b[name]
        d = y["mw"] - x["mw"]
        rows.append((name, x, y, d))
        print("%-8s | MW %6.1f HA %5.1f 落区 %.2f 生成 %3d | MW %6.1f HA %5.1f 落区 %.2f 生成 %3d | %+7.1f"
              % (name, x["mw"], x["ha"], x["in_band"], x["n"],
                 y["mw"], y["ha"], y["in_band"], y["n"], d))
    if rows:
        print("\nΔMW 均值 %+.1f（n=%d 对）；v10 的 MW 极差 %.1f，v10term 的 MW 极差 %.1f"
              % (st.mean([r[3] for r in rows]), len(rows),
                 max(r[1]["mw"] for r in rows) - min(r[1]["mw"] for r in rows),
                 max(r[2]["mw"] for r in rows) - min(r[2]["mw"] for r in rows)))
        dip_v10 = [r for r in rows if r[1]["mw"] < 340]
        dip_term = [r for r in rows if r[2]["mw"] < 340]
        print("掉出尺寸下限（MW<340）的点：v10 %d 个 %s | v10term %d 个 %s"
              % (len(dip_v10), [r[0] for r in dip_v10], len(dip_term), [r[0] for r in dip_term]))
        print("\n判读：若 v10term 的 MW 极差明显更小、且掉出下限的点更少 → **训练侧 pos_weight 确实抑制了尺寸掉坑**；")
        print("      若两者相当 → 该干预在生成级无效（只能靠 V2 逐样本平衡 / S1·S2 采样侧机制）。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
