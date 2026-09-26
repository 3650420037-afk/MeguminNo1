# -*- coding: utf-8 -*-
"""审计 GPCR 数据集 SDF 构象质量：sanitize失败/键长异常/原子重叠/共面塌缩"""
import os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import DATA, src_on_path
src_on_path()
from rdkit import Chem
import numpy as np

def audit_sdf(path):
    issues = []
    mol = next(iter(Chem.SDMolSupplier(path, removeHs=False, sanitize=False)), None)
    if mol is None:
        return ["parse_fail"]
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_PROPERTIES)
    except Exception as e:
        issues.append("sanitize:" + type(e).__name__)
    try:
        conf = mol.GetConformer()
        pos = np.array(conf.GetPositions())
    except Exception:
        return issues + ["no_conformer"]
    n = pos.shape[0]
    # 拓扑孤键: 无任何键连接的原子 (会导致训练 LigandBFSMask KeyError)
    bonded = set()
    for b in mol.GetBonds():
        bonded.add(b.GetBeginAtomIdx()); bonded.add(b.GetEndAtomIdx())
    orphans = [i for i in range(n) if i not in bonded]
    if orphans:
        issues.append("orphan_atoms=%s" % orphans[:5])
    # 键长
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if b.GetBeginAtom().GetAtomicNum() == 1 or b.GetEndAtom().GetAtomicNum() == 1:
            continue
        d = float(np.linalg.norm(pos[i] - pos[j]))
        if d < 0.9 or d > 2.2:
            issues.append("bondlen_%.2f_%d_%d" % (d, i, j))
    # 非键原子重叠（只查重原子，间隔采样防 O(n^2) 过大）
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    if len(heavy) <= 60:
        for a in range(len(heavy)):
            for b2 in range(a + 1, len(heavy)):
                i, j = heavy[a], heavy[b2]
                bonded = mol.GetBondBetweenAtoms(i, j) is not None
                d = float(np.linalg.norm(pos[i] - pos[j]))
                if not bonded and d < 0.7:
                    issues.append("overlap_%.2f_%d_%d" % (d, i, j))
    # 平面性检查 (PCA, 与坐标轴无关): 旧写法用 pos.var(axis=0).min()<1e-3 判共面,
    # 依赖分子是否恰好与坐标轴对齐 —— 实测 benzene 误报、naphthalene 漏检。
    # 改用 SVD 最小主轴占比; 且只有"含 sp3 中心却整体近平面"才算塌缩
    # (纯芳香平面分子本就该是平面, 不应误报)。
    if n >= 6:
        try:
            centered = pos - pos.mean(axis=0)
            sv = np.linalg.svd(centered, compute_uv=False)
            if sv[0] > 1e-9 and sv[2] / sv[0] < 0.02:
                has_sp3 = any(
                    a.GetHybridization() == Chem.HybridizationType.SP3 and a.GetAtomicNum() != 1
                    for a in mol.GetAtoms()
                )
                if has_sp3:
                    issues.append("planar_collapse")
        except Exception:
            pass
    return issues

def main(dataset):
    d = os.path.join(DATA, dataset)
    index = pickle.load(open(os.path.join(d, "index.pkl"), "rb"))
    bad = {}
    for pocket, lig, _, _ in index:
        p = os.path.join(d, lig)
        iss = audit_sdf(p)
        if iss:
            bad[lig] = iss
    print("%s: 总配对=%d 异常=%d (%.1f%%)" % (dataset, len(index), len(bad), 100.0*len(bad)/max(1,len(index))))
    for lig, iss in sorted(bad.items()):
        print("  %s: %s" % (lig, iss[:4]))
    return bad

if __name__ == "__main__":
    all_bad = {}
    for ds in sys.argv[1:] or ["gpcr_a2a", "gpcr_b2ar", "gpcr_multitarget"]:
        all_bad[ds] = main(ds)
