# -*- coding: utf-8 -*-
"""凹坑测试分析：机制能不能把"掉进小模式"的检查点拉回带内（判据用预登记那套）。

对照臂 = 今夜同协议 v10@6500（`v10_matched_sel/6500`，MW 302.5 / HA 22.0）。
处理臂 = outputs/dip_test/{s1on,s2floor,s3window}_s{2024,2025}。

判据（预登记 §5.41，不自创）：① MW 中位 ≥340 且 HA 中位 ≥25；② strict ≥ 官方同协议臂的 90%；
③ 生成数 ≥ 请求值 n；④ 跨种子极差 ≤10%（§5.33 处方①）。
"""
import csv
import glob
import io
import os
import statistics as st
import subprocess
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))
BASE = os.path.join(ROOT, "outputs", "dip_test")
NIGHT = os.path.join(ROOT, "outputs", "night_pipeline_v10")

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")
BAND = (340.0, 430.0)
N_REQ = 60


def find_smiles(d):
    c = glob.glob(os.path.join(d, "**", "SMILES.txt"), recursive=True)
    return sorted(c)[-1] if c else None


def mw_ha(smi):
    ms = [Chem.MolFromSmiles(l.strip()) for l in io.open(smi, encoding="utf-8") if l.strip()]
    ms = [m for m in ms if m]
    if not ms:
        return None
    ha = sorted(m.GetNumHeavyAtoms() for m in ms)
    mw = sorted(round(Descriptors.MolWt(m), 1) for m in ms)
    return dict(n=len(ms), mw=st.median(mw), ha=st.median(ha),
                inband=sum(1 for x in mw if BAND[0] <= x <= BAND[1]) / len(mw),
                under=sum(1 for x in mw if x < BAND[0]) / len(mw),
                over=sum(1 for x in mw if x > BAND[1]) / len(mw))


def main():
    rows = []
    control = find_smiles(os.path.join(NIGHT, "v10_matched_sel", "6500"))
    if control:
        m = mw_ha(control)
        rows.append(("对照 v10@6500（无机制）", m))
    for arm in ("s1on", "s2floor", "s3window"):
        vals = []
        for seed in (2024, 2025):
            smi = find_smiles(os.path.join(BASE, "%s_s%d" % (arm, seed)))
            if smi:
                m = mw_ha(smi)
                if m:
                    vals.append(m)
                    rows.append(("%s seed %d" % (arm, seed), m))
        if len(vals) == 2:
            spread = 100.0 * abs(vals[0]["mw"] - vals[1]["mw"]) / st.mean([v["mw"] for v in vals])
            rows.append(("  %s 跨种子 MW 极差" % arm, dict(n=0, mw=spread, ha=0, inband=0, under=0, over=0)))

    print("%-26s %4s %9s %7s %8s %8s %8s" % ("臂", "n", "MW中位", "HA", "落区率", "MW<340", "MW>430"))
    for tag, m in rows:
        if m["n"] == 0:
            print("%-26s %4s %8.1f%%" % (tag, "-", m["mw"]))
            continue
        print("%-26s %4d %9.1f %7.1f %8.2f %8.2f %8.2f"
              % (tag, m["n"], m["mw"], m["ha"], m["inband"], m["under"], m["over"]))
    print("\n判据①（MW 中位 ≥340 且 HA ≥25）：看上面 MW中位/HA 两列")
    print("判据③（生成数 ≥%d）：看 n 列；判据④（跨种子极差 ≤10%%）：看极差行" % N_REQ)
    print("判据②（strict ≥ 官方同协议臂 90%）需跑 arm_gen_metrics（见下）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
