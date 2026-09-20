# -*- coding: utf-8 -*-
"""库级构象弛豫: 对入库 SDF 做 MMFF94s 优化, 记录弛豫能, 输出干净 3D 库
注意: 中文路径下 RDKit 文件 API 不可用 -> 全部走 Python open() + MolBlock
"""
import os, csv, glob, sys
import numpy as np
from multiprocessing import Pool
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

LIB = r"D:\MMModel\化合物库"
SRC = os.path.join(LIB, "sdf")
DST = os.path.join(LIB, "sdf_relaxed")
os.makedirs(DST, exist_ok=True)

def relax_one(p):
    name = os.path.basename(p)
    try:
        block = open(p, encoding="utf-8", errors="ignore").read()
        m = Chem.MolFromMolBlock(block, removeHs=False, sanitize=True)
        if m is None:
            return (name, None, None, "parse_fail")
        mh = Chem.AddHs(m)
        # MMFF94s 优先, 缺参数回退 UFF
        ff_kind = "MMFF"
        props = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant="MMFF94s")
        if props is None:
            ff_kind = "UFF"
            ff = AllChem.UFFGetMoleculeForceField(mh)
        else:
            ff = AllChem.MMFFGetMoleculeForceField(mh, props)
        if ff is None:
            return (name, None, None, "ff_fail")
        e0 = ff.CalcEnergy()
        ff.Minimize(maxIts=500)
        e1 = ff.CalcEnergy()
        m_opt = Chem.RemoveHs(mh)
        out = os.path.join(DST, name)
        open(out, "w", encoding="utf-8").write(Chem.MolToMolBlock(m_opt))
        return (name, round(e0 - e1, 2), round(e1, 2), ff_kind)
    except Exception as e:
        return (name, None, None, "err:" + type(e).__name__)

if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(SRC, "*.sdf")))
    print("待弛豫:", len(files))
    rows = []
    with Pool(6) as pool:
        for i, r in enumerate(pool.imap_unordered(relax_one, files, chunksize=8), 1):
            rows.append(r)
            if i % 500 == 0:
                print("  ...%d/%d" % (i, len(files)))
    ok = [(n, s, e, k) for n, s, e, k in rows if s is not None]
    bad = [(n, k) for n, s, e, k in rows if s is None]
    with open(os.path.join(LIB, "strain_energy.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["sdf", "relax_energy_kcal", "final_energy_kcal", "forcefield"])
        for r in sorted(ok):
            w.writerow(r)
    s = np.array([x[1] for x in ok])
    print()
    print("成功 %d | 失败 %d" % (len(ok), len(bad)))
    if len(s):
        print("弛豫能: 中位 %.2f | 均值 %.2f | P90 %.2f | max %.2f" % (np.median(s), s.mean(), np.percentile(s, 90), s.max()))
        print("分布: <5 %.0f%% | 5-20 %.0f%% | 20-50 %.0f%% | >50 %.0f%%" % (
            100*(s < 5).mean(), 100*((s >= 5) & (s < 20)).mean(), 100*((s >= 20) & (s < 50)).mean(), 100*(s >= 50).mean()))
    if bad:
        print("失败样例:", bad[:5])
    print("输出目录:", DST)
