import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit import Geometry


class MolReconsError(Exception):
    pass

def add_context(data):
    data.ligand_context_pos = data.ligand_pos
    data.ligand_context_element = data.ligand_element
    data.ligand_context_bond_index = data.ligand_bond_index
    data.ligand_context_bond_type = data.ligand_bond_type
    return data

def reconstruct_from_generated_with_edges(data, raise_error=True, sanitize=True):
    xyz = data.ligand_context_pos.clone().cpu().tolist()
    atomic_nums = data.ligand_context_element.clone().cpu().tolist()
    # indicators = data.ligand_context_feature_full[:, -len(ATOM_FAMILIES_ID):].clone().cpu().bool().tolist()
    bond_index = data.ligand_context_bond_index.clone().cpu().tolist()
    bond_type = data.ligand_context_bond_type.clone().cpu().tolist()
    n_atoms = len(atomic_nums)

    rd_mol = Chem.RWMol()
    rd_conf = Chem.Conformer(n_atoms)
    
    # add atoms and coordinates
    for i, atom in enumerate(atomic_nums):
        rd_atom = Chem.Atom(atom)
        rd_mol.AddAtom(rd_atom)
        rd_coords = Geometry.Point3D(*xyz[i])
        rd_conf.SetAtomPosition(i, rd_coords)
    rd_mol.AddConformer(rd_conf)
    
    # add bonds
    for i, type_this in enumerate(bond_type):
        node_i, node_j = bond_index[0][i], bond_index[1][i]
        if node_i < node_j:
            if type_this == 1:
                rd_mol.AddBond(node_i, node_j, Chem.BondType.SINGLE)
            elif type_this == 2:
                rd_mol.AddBond(node_i, node_j, Chem.BondType.DOUBLE)
            elif type_this == 3:
                rd_mol.AddBond(node_i, node_j, Chem.BondType.TRIPLE)
            elif type_this == 12:
                rd_mol.AddBond(node_i, node_j, Chem.BondType.AROMATIC)
            else:
                raise Exception('unknown bond order {}'.format(type_this))
    
    # modify
    try:
        rd_mol = modify_submol(rd_mol)
    except:
        if raise_error:
            raise MolReconsError()
        else:
            print('MolReconsError')
    # check valid
    rd_mol_check = Chem.MolFromSmiles(Chem.MolToSmiles(rd_mol))
    if rd_mol_check is None:
        if raise_error:
            raise MolReconsError()
        else:
            print('MolReconsError')
    
    rd_mol = rd_mol.GetMol()
    if 12 in bond_type:  # mol may directlu come from ture mols and contains aromatic bonds
        Chem.Kekulize(rd_mol, clearAromaticFlags=True)
    if sanitize:
        try:
            Chem.SanitizeMol(rd_mol, Chem.SANITIZE_ALL^Chem.SANITIZE_KEKULIZE^Chem.SANITIZE_SETAROMATICITY)
        except:
            Chem.SanitizeMol(rd_mol)
           
    return rd_mol


def modify_submol(mol):  # modify mols containing C=N(C)O
    submol = Chem.MolFromSmiles('C=N(C)O', sanitize=False)
    sub_fragments = mol.GetSubstructMatches(submol)
    for fragment in sub_fragments:
        atomic_nums = np.array([mol.GetAtomWithIdx(atom).GetAtomicNum() for atom in fragment])
        idx_atom_N = fragment[np.where(atomic_nums == 7)[0][0]]
        idx_atom_O = fragment[np.where(atomic_nums == 8)[0][0]]
        mol.GetAtomWithIdx(idx_atom_N).SetFormalCharge(1)  # set N to N+
        mol.GetAtomWithIdx(idx_atom_O).SetFormalCharge(-1)  # set O to O-
    return mol


def relax_mol_geometry(mol, max_iters=500, keep_pose=True, n_confs=4):
    """修正模型生成构象的局部几何畸变, 同时保留其在口袋中的姿态。

    背景: 模型逐原子生成只保证键长级局部几何, 部分分子的键角畸变会导致
    RDKit AddHs 位置计算退化(出现 H-H 重合, MMFF 能量爆炸到 1e9), 此时
    直接对模型坐标做力场优化无法恢复。本函数改为:
      1) 用 ETKDG 重新生成内部几何合理的构象(带 H, 多构象取最优)
      2) 刚体叠合(AlignMol)到模型坐标 -> 保留模型预测的口袋姿态
      3) MMFF94s 优化(从合理起点, 稳定收敛)
    返回 (优化后的分子, 弛豫能 kcal/mol 或 None); 任何失败均回退原分子不抛异常。
    """
    try:
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdMolAlign
        base = Chem.Mol(mol)
        if base.GetNumConformers() == 0:
            return mol, None
        ref_pos = base.GetConformer().GetPositions().copy()
        heavy_idx = [a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() != 1]

        # 1) ETKDG 生成合理内部几何 (在加 H 的副本上)
        mh = Chem.AddHs(base, addCoords=False)
        ok = AllChem.EmbedMultipleConfs(mh, numConfs=n_confs, randomSeed=42,
                                        useRandomCoords=False, maxAttempts=50)
        if len(ok) == 0:
            ok = AllChem.EmbedMultipleConfs(mh, numConfs=n_confs, randomSeed=42, useRandomCoords=True)
        if len(ok) == 0:
            return mol, None

        # 2) 叠合到模型姿态: 选叠合后重原子 RMSD 最小的构象
        best, best_rms = None, None
        for cid in ok:
            probe = Chem.Mol(mh)
            amap = [(i, i) for i in heavy_idx]
            try:
                rms = rdMolAlign.AlignMol(probe, base, prbCid=cid, atomMap=amap)
            except Exception:
                continue
            if best_rms is None or rms < best_rms:
                best, best_rms = probe, rms
        if best is None:
            return mol, None

        # 3) MMFF94s 优化 (从合理起点)
        props = AllChem.MMFFGetMoleculeProperties(best, mmffVariant="MMFF94s")
        ff = AllChem.MMFFGetMoleculeForceField(best, props) if props is not None else None
        if ff is None:
            try:
                ff = AllChem.UFFGetMoleculeForceField(best)
            except Exception:
                return mol, None
        if ff is None:
            return mol, None
        e0 = ff.CalcEnergy()
        if not np.isfinite(e0) or e0 > 1e6:      # 起点仍异常 -> 放弃
            return mol, None
        ff.Minimize(maxIts=max_iters)
        e1 = ff.CalcEnergy()
        # 保留氢直接交付: 无氢 SDF 被下游 AddHs 时会重新推算氢位置, 在模型姿态上可能退化出 H 重合
        m_opt = best
        strain = float(e0 - e1)
        if not np.isfinite(strain):
            return mol, None
        return m_opt, (float(best_rms) if best_rms is not None else None, strain)
    except Exception:
        return mol, None

