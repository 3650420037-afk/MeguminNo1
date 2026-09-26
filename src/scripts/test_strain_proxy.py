# -*- coding: utf-8 -*-
"""A4 可行性原型: 未完成(生长中)分子的应变能代理计算 —— 速度/成功率/区分度"""
import os, sys, time, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import OUTPUTS, RESULTS, src_on_path
src_on_path()
import numpy as np
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
from utils import misc  # torch.load 补丁
from utils.guidance import build_partial_mol

def strain_proxy(data):
    """部分分子构象紧张度: MMFF 单点能 / 重原子数 (kcal/mol/atom)"""
    mol = build_partial_mol(data)
    if mol is None:
        return None
    try:
        mh = Chem.AddHs(mol)
        props = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant="MMFF94s")
        if props is None:
            return None
        ff = AllChem.MMFFGetMoleculeForceField(mh, props)
        if ff is None:
            return None
        e = ff.CalcEnergy()
        n = max(1, mh.GetNumHeavyAtoms())
        return e / n
    except Exception:
        return None

# 取一个采样快照的 queue 候选 (生长中分子)
snaps = glob.glob(os.path.join(OUTPUTS, "reverify_baseline", "*", "samples_all.pt"))
if not snaps:
    snaps = glob.glob(os.path.join(OUTPUTS, "guided_5010050_l3", "*", "samples_all.pt"))
pool = torch.load(snaps[0], map_location="cpu", weights_only=False)
cands = pool["queue"][:60]
print("快照:", os.path.basename(os.path.dirname(snaps[0])), "| 候选:", len(cands))

t0 = time.time()
vals, fail = [], 0
for d in cands:
    v = strain_proxy(d)
    if v is None:
        fail += 1
    else:
        vals.append(v)
dt = time.time() - t0
print("部分分子: 成功 %d | 失败 %d | 总耗时 %.2fs | 平均 %.1f ms/候选" % (
    len(vals), fail, dt, 1000*dt/max(1, len(cands))))
if vals:
    a = np.array(vals)
    print("  能量/重原子 (kcal/mol/atom): 中位 %.2f | P10 %.2f | P90 %.2f | max %.2f" % (
        np.median(a), np.percentile(a, 10), np.percentile(a, 90), a.max()))

# 对照: 闭合完成的分子 (库 SDF 抽样 60 个)
files = sorted(glob.glob(os.path.join(RESULTS, "sdf", "*.sdf")))[:60]
cv = []
for p in files:
    m = Chem.MolFromMolBlock(open(p, encoding="utf-8", errors="ignore").read(), removeHs=False, sanitize=True)
    if m is None: continue
    mh = Chem.AddHs(m)
    props = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant="MMFF94s")
    if props is None: continue
    ff = AllChem.MMFFGetMoleculeForceField(mh, props)
    cv.append(ff.CalcEnergy() / max(1, mh.GetNumHeavyAtoms()))
c = np.array(cv)
print("闭合成品 (对照): n=%d 中位 %.2f | P10 %.2f | P90 %.2f" % (
    len(c), np.median(c), np.percentile(c, 10), np.percentile(c, 90)))
print()
print("判定:", "可接入引导层" if len(vals) > 0.8*len(cands) and dt/max(1,len(cands)) < 0.01 else "降级为库级过滤")
