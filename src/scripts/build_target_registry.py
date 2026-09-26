# -*- coding: utf-8 -*-
"""从 data/targets/ 的目标结构自动生成靶点登记表 configs/targets.json。

对每个靶点提取(全部来自结构文件本身, 不手工录入, 避免写错):
  receptor      受体名称(取 _struct.title / PDB TITLE)
  pdb           受体坐标文件(相对仓库根)
  ligand_comp   共晶配体的化学组分编号
  ligand_smiles 该配体的 SMILES(取条目自带 CCD 的 _chem_comp_atom/_chem_comp_bond)
  center        口袋中心 = 共晶配体重原子的质心(实验坐标, 即真实结合位点)
  n_ligand_heavy 共晶配体重原子数

用法:
    python src/scripts/build_target_registry.py
输出:
    configs/targets.json
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, TARGETS, CONFIGS, ensure_dir  # noqa: E402

import warnings

warnings.simplefilter("ignore")

# 靶点名 -> (受体坐标文件, 含配体坐标/定义的文件, 共晶配体 comp_id)
TARGET_DEFS = [
    ("A2A", "4EIY_A2A受体.pdb", "4EIY_配体ZMA.cif", "ZMA"),
    ("B2AR", "2RH1_β2肾上腺素受体.pdb", "2RH1_配体CAU.cif", "CAU"),
    ("D3", "3PBL_D3多巴胺受体.pdb", "3PBL.cif", "ETQ"),
    ("5HT2B", "4IB4_5HT2B受体.pdb", "4IB4.cif", "ERM"),
    ("A1AR", "5UEN.pdb", "5UEN.cif", "DU1"),
]

BOND_ORDER = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 1.5, "DELO": 1}


def cif_title(path):
    """从 mmCIF 的 _struct.title 取受体名称。"""
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        d = MMCIF2Dict(path)
        t = d.get("_struct.title")
        if isinstance(t, list):
            t = " ".join(t)
        return (t or "").strip()
    except Exception:
        return ""


def pdb_title(path):
    """从 PDB 的 TITLE 记录取受体名称。"""
    out = []
    try:
        for line in open(path, encoding="utf-8", errors="ignore"):
            if line.startswith("TITLE"):
                out.append(line[10:].strip())
            elif out and not line.startswith("TITLE"):
                break
    except Exception:
        pass
    return " ".join(out).strip()


def smiles_from_ccd(path, comp_id):
    """从条目自带的 CCD 块(_chem_comp_atom/_chem_comp_bond)重建 SMILES。"""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    d = MMCIF2Dict(path)

    def g(k):
        v = d.get(k, [])
        return v if isinstance(v, list) else [v]

    ids, els = g("_chem_comp_atom.atom_id"), g("_chem_comp_atom.type_symbol")
    comps = g("_chem_comp_atom.comp_id")
    rw, idx = Chem.RWMol(), {}
    for i, (a, e) in enumerate(zip(ids, els)):
        if comps and comps[i] != comp_id:
            continue
        if str(e).upper() == "H":
            continue
        idx[a] = rw.AddAtom(Chem.Atom(str(e).capitalize()))
    a1, a2 = g("_chem_comp_bond.atom_id_1"), g("_chem_comp_bond.atom_id_2")
    vo, bc = g("_chem_comp_bond.value_order"), g("_chem_comp_bond.comp_id")
    for i in range(len(a1)):
        if bc and bc[i] != comp_id:
            continue
        if a1[i] not in idx or a2[i] not in idx:
            continue
        order = str(vo[i]).upper()
        bt = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
              3: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}.get(
                  BOND_ORDER.get(order, 1), Chem.BondType.SINGLE)
        rw.AddBond(idx[a1[i]], idx[a2[i]], bt)
        if bt == Chem.BondType.AROMATIC:
            for a in (a1[i], a2[i]):
                rw.GetAtomWithIdx(idx[a]).SetIsAromatic(True)
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
        return Chem.MolToSmiles(m)
    except Exception:
        return ""


def ligand_centroid_from_cif(path, comp_id):
    """共晶配体重原子质心(mmCIF atom_site)。

    注意: 只含 chem_comp 定义(无 _atom_site)的精简 cif 无法解析坐标,
    此时返回 (None, 0), 由调用方回退到 PDB 的 HETATM。
    """
    import numpy as np
    from Bio.PDB.MMCIFParser import MMCIFParser
    with open(path, encoding="utf-8", errors="ignore") as fh:
        if "_atom_site.id" not in fh.read():
            return None, 0
    s = MMCIFParser(QUIET=True).get_structure("x", path)
    for res in s.get_residues():
        if not str(res.id[0]).strip().startswith("H_"):
            continue
        if res.get_resname().strip() != comp_id:
            continue
        pts = [a.coord for a in res.get_atoms()
               if (getattr(a, "element", "") or "").strip().upper() not in ("", "H")]
        if len(pts) >= 5:
            return np.asarray(pts, dtype=float).mean(axis=0), len(pts)
    return None, 0


def ligand_centroid_from_pdb(path, comp_id):
    """共晶配体重原子质心(PDB HETATM, 按元素列/原子名判重原子)。"""
    import numpy as np
    pts = []
    for line in open(path, encoding="utf-8", errors="ignore"):
        if not (line.startswith("HETATM") and line[17:20].strip() == comp_id):
            continue
        elem = line[76:78].strip() or line[12:16].strip()[:1]
        if elem.upper() == "H":
            continue
        try:
            pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
    if len(pts) < 5:
        return None, 0
    return np.asarray(pts, dtype=float).mean(axis=0), len(pts)


def main():
    reg = {}
    for name, pdb_name, cif_name, comp in TARGET_DEFS:
        pdb = os.path.join(TARGETS, pdb_name)
        cif = os.path.join(TARGETS, cif_name)
        entry = {"pdb": os.path.relpath(pdb, ROOT).replace("\\", "/"),
                 "ligand_comp": comp}
        # 受体名
        title = cif_title(cif) if os.path.exists(cif) else ""
        if not title and os.path.exists(pdb):
            title = pdb_title(pdb)
        entry["receptor"] = title or name
        # 配体 SMILES
        entry["ligand_smiles"] = smiles_from_ccd(cif, comp) if os.path.exists(cif) else ""
        # 口袋中心
        center, n = (None, 0)
        if os.path.exists(cif):
            center, n = ligand_centroid_from_cif(cif, comp)
        if center is None and os.path.exists(pdb):
            center, n = ligand_centroid_from_pdb(pdb, comp)
        if center is None:
            print("  [跳过] %-10s 无法定位共晶配体 %s" % (name, comp))
            continue
        entry["center"] = [round(float(v), 2) for v in center]
        entry["n_ligand_heavy"] = int(n)
        # 一致性披露: 沉积坐标的重原子数 vs CCD 定义的重原子数。
        # 二者不等说明该条目的共晶配体未被完整建模(例如 5UEN/A1AR 的 DU172 漏建了
        # 磺酰氟上的 F), 建库时这类条目无法同时满足"用实验坐标"与"用完整分子"。
        try:
            from rdkit import Chem, RDLogger
            RDLogger.DisableLog("rdApp.*")
            m = Chem.MolFromSmiles(entry["ligand_smiles"]) if entry["ligand_smiles"] else None
            entry["n_smiles_heavy"] = m.GetNumAtoms() if m is not None else -1
        except Exception:
            entry["n_smiles_heavy"] = -1
        entry["deposited_ligand_complete"] = (entry["n_smiles_heavy"] == int(n))
        if not entry["deposited_ligand_complete"]:
            print("    ! %s 的沉积配体不完整: 沉积 %d 个重原子, CCD 定义 %d 个"
                  % (comp, n, entry["n_smiles_heavy"]))
        entry["source"] = os.path.relpath(cif if os.path.exists(cif) else pdb,
                                          ROOT).replace("\\", "/")
        reg[name] = entry
        print("  %-10s %-22s 配体=%-4s 重原子=%-3d 中心=%s"
              % (name, entry["receptor"][:22], comp, n, entry["center"]))

    ensure_dir(CONFIGS)
    out = os.path.join(CONFIGS, "targets.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(reg, f, ensure_ascii=False, indent=1)
    print("已写出 %s (%d 个靶点)" % (out, len(reg)))
    if not reg:
        raise SystemExit("未生成任何靶点, 请检查 data/targets/ 是否完整")


if __name__ == "__main__":
    main()
