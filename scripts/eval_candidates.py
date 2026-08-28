# -*- coding: utf-8 -*-
"""候选质量评估：对比 A/B 双方快照中全部未完成候选的 logp/生长进度/类药性（零GPU）"""
import os, sys, glob
sys.path.insert(0, r"D:\MMModel\Pocket2Mol")
import torch, numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, QED

def partial_mol(data):
    try:
        pos = data.ligand_context_pos.numpy()
        els = data.ligand_context_element.numpy().tolist()
        bi = data.ligand_context_bond_index.numpy().tolist()
        bt = data.ligand_context_bond_type.numpy().tolist()
    except Exception:
        return None
    m = Chem.RWMol()
    for z in els:
        m.AddAtom(Chem.Atom(int(z)))
    conf = Chem.Conformer(len(els))
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
    m.AddConformer(conf)
    btmap = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for k, t in enumerate(bt):
        i, j = bi[0][k], bi[1][k]
        if i < j and int(t) in btmap:
            try: m.AddBond(i, j, btmap[int(t)])
            except Exception: return None
    mol = m.GetMol()
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_SETAROMATICITY ^ Chem.SANITIZE_PROPERTIES, catchErrors=True)
    except Exception:
        return None
    return mol

def eval_pool(path):
    pool = torch.load(path, map_location="cpu", weights_only=False)
    rows = []
    for d in pool["queue"]:
        steps = len(d.logp_focal)
        logp_sum = float(np.sum(d.average_logp))
        natoms = int(d.ligand_context_pos.shape[0])
        mol = partial_mol(d)
        mw = logp_crippen = qed = float("nan")
        rings = 0
        if mol is not None:
            try:
                mw = Descriptors.MolWt(mol)
                logp_crippen = Crippen.MolLogP(mol)
                rings = Descriptors.RingCount(mol)
                try: qed = QED.qed(Chem.MolFromSmiles(Chem.MolToSmiles(mol)) or mol)
                except Exception: pass
            except Exception:
                pass
        rows.append((steps, natoms, logp_sum, mw, logp_crippen, qed, rings))
    return rows

def stats(vals):
    v = np.array([x for x in vals if not np.isnan(x)], dtype=float)
    if len(v) == 0: return "n/a"
    return "mean=%.2f med=%.2f n=%d" % (v.mean(), np.median(v), len(v))

if __name__ == "__main__":
    base = r"D:\MMModel\Pocket2Mol\outputs"
    runs = [("baseline", glob.glob(os.path.join(base, "ab_305020_baseline", "*", "samples_all.pt")))]
    runs += [("baseline_s%d" % s, glob.glob(os.path.join(base, "ab_305020_baseline_s%d" % s, "*", "samples_all.pt"))) for s in (2025, 2026)]
    runs += [("local", glob.glob(os.path.join(base, "ab_305020_local", "*", "samples_all.pt")))]
    runs += [("local_s%d" % s, glob.glob(os.path.join(base, "ab_305020_local_s%d" % s, "*", "samples_all.pt"))) for s in (2025, 2026)]
    groups = {}
    for name, files in runs:
        if not files: continue
        rows = []
        for f in files:
            rows += eval_pool(f)
        groups[name] = rows
    print("候选数/组:", {k: len(v) for k, v in groups.items()})
    print()
    for metric, idx, label in [("生长步数", 0, "steps"), ("原子数", 1, "natoms"), ("总logp", 2, "logp_sum"),
                               ("MW", 3, "mw"), ("CrippenLogP", 4, "clogp"), ("QED", 5, "qed"), ("环数", 6, "rings")]:
        line = "%-10s" % metric
        for side in ["baseline", "baseline_s2025", "baseline_s2026", "local", "local_s2025", "local_s2026"]:
            if side in groups:
                line += " | %s %s" % (side.replace("_s", "@"), stats([r[idx] for r in groups[side]]))
        print(line)
    print()
    # 类药通过率（MW<500, LogP<5, QED>0.4, 环>=1）
    for side in ["baseline", "baseline_s2025", "baseline_s2026", "local", "local_s2025", "local_s2026"]:
        rows = groups.get(side, [])
        q = [r for r in rows if not np.isnan(r[5])]
        ok = [r for r in q if r[3] < 500 and r[4] < 5 and r[5] > 0.4 and r[6] >= 1]
        print("%-16s 可评估QED=%3d/%3d 类药通过=%d (%.0f%%)" % (
            side, len(q), len(rows), len(ok), 100.0*len(ok)/max(1, len(q))))
