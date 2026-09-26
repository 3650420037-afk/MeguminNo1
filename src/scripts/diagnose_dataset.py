# -*- coding: utf-8 -*-
"""gpcr_v3_merged 数据集只读诊断 (Windows / 中文路径安全 / 可重复运行)

用途
----
逐个样本体检 data/gpcr_v3_merged, 定位"非有限损失(NaN/Inf)"的成因, 并给出可执行的过滤规则。

检查项
------
1. 口袋 PDB: ATOM/HETATM 数, 残基数, 是否有非标准残基名(会导致 PDBProtein 抛 KeyError 而被静默丢弃)
2. 配体 SDF: RDKit 可读性, 声明原子数 vs 实际, 重原子数, 3D 构象数, 碎片数, 孤立原子数,
   芳香键标志(4)条数 (会污染 valence/tri_edge 特征), 重复坐标
3. 坐标合法性: NaN / 非有限 / |坐标| > 1e4
4. 配体是否在口袋内: 重原子落在"口袋包围盒外扩 2 A"内的比例, 配体到最近口袋原子的距离,
   配体质心到口袋质心距离, 重原子中"4 A 内有口袋原子"的比例 (对应 FocalBuilder 的 r=4.0 初筛)
5. mask 覆盖度: 复现 utils/transforms.py 的 LigandMixedMask, 统计"配体上下文为空"的概率
   —— 这是训练中 loss_frontier/loss_edge 在空张量上归约 -> NaN 的直接来源
6. 训练可用性: 与 <dataset>_name2id.pt 比对, 判断样本是否真的进入了 lmdb

只读保证: 仅读数据文件; 仅写 --out 指定的 CSV 与 --summary-out 指定的文本报告。
不训练, 不采样, 不修改数据集。

用法
----
python scripts/diagnose_dataset.py
python scripts/diagnose_dataset.py --data data/gpcr_v3_merged --mask-draws 400
"""
import argparse
import csv
import math
import os
import pickle
import random
import sys
import warnings
from collections import Counter, OrderedDict

# ---- 中文路径安全: 全部用绝对路径, 不做 chdir ----
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import CONFIGS, DATA, src_on_path
src_on_path()

warnings.simplefilter("ignore")

import numpy as np  # noqa: E402

from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

AA_STD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

CSV_COLUMNS = [
    "idx", "split", "source", "pocket_file", "ligand_file",
    "pocket_exists", "ligand_exists",
    "pocket_atom_count", "pocket_hetatm_count", "pocket_hydrogen_count",
    "pocket_residue_count", "pocket_chain_count", "pocket_unknown_resnames",
    "pocket_parse_ok", "pocket_parse_error",
    "pocket_nonfinite_coords", "pocket_max_abs_coord",
    "pocket_bbox_dx", "pocket_bbox_dy", "pocket_bbox_dz",
    "ligand_readable", "ligand_parse_ok", "ligand_parse_error",
    "ligand_declared_atoms", "ligand_declared_bonds",
    "ligand_atom_count", "ligand_heavy_atom_count", "ligand_bond_count",
    "ligand_num_conformers", "ligand_is3d", "ligand_num_fragments", "ligand_orphan_atoms",
    "ligand_aromatic_flag_bonds", "ligand_duplicate_coords",
    "ligand_nonfinite_coords", "ligand_max_abs_coord",
    "ligand_centroid_x", "ligand_centroid_y", "ligand_centroid_z",
    "pocket_centroid_x", "pocket_centroid_y", "pocket_centroid_z",
    "pocket_ligand_centroid_distance", "ligand_min_dist_to_pocket",
    "frac_heavy_inside_pocket_bbox_pad2", "frac_heavy_inside_pocket_bbox_pad0",
    "frac_heavy_within_4A_of_pocket",
    "in_name2id",
    "mask_p_empty_context_theory", "mask_p_empty_context_sim",
    "flags",
]


# --------------------------------------------------------------------------- #
# 基础工具
# --------------------------------------------------------------------------- #
def fnum(value, ndigits=4):
    """把 float/np 数值格式化成 CSV 友好的字符串, 非有限或缺失返回空串。"""
    if value is None:
        return ""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(value):
        return ""
    return round(value, ndigits)


def is_source_legacy(pocket_name):
    """老数据 (gpcr_multitarget_v2 并入) 的口袋名形如 gpcr_a2a_pocket.pdb; 新抓取为 pocket_<PDBID>_<LIG>.pdb"""
    base = os.path.basename(pocket_name)
    return not base.startswith("pocket_")


def scan_pdb_lines(path):
    """轻量文本扫描: HETATM 数, 氢原子数, 非标准残基名。只读。"""
    hetatm = 0
    hydro = 0
    unknown = set()
    n_atom_lines = 0
    try:
        with open(path, "r", encoding="ascii", errors="ignore") as fh:
            for line in fh:
                rec = line[0:6].strip()
                if rec == "HETATM":
                    hetatm += 1
                elif rec == "ATOM":
                    n_atom_lines += 1
                    resname = line[17:20].strip()
                    if resname and resname not in AA_STD:
                        unknown.add(resname)
                    elem = line[76:78].strip().upper()
                    if not elem:
                        elem = line[12:16].strip()[:1].upper()
                    if elem == "H":
                        hydro += 1
    except OSError:
        return None
    return {"hetatm": hetatm, "hydro": hydro, "unknown": sorted(unknown),
            "atom_lines": n_atom_lines}


def read_sdf_raw_text(path):
    """按 utils/protein_ligand.parse_sdf_file 的口径读 SDF 文本块, 返回统计量。"""
    out = {
        "declared_atoms": None, "declared_bonds": None,
        "aromatic_flag_bonds": None, "raw_error": None,
    }
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.read().splitlines()
        if len(lines) < 5:
            out["raw_error"] = "fewer than 5 lines"
            return out
        na, nb = int(lines[3][0:3]), int(lines[3][3:6])
        out["declared_atoms"], out["declared_bonds"] = na, nb
        arom = 0
        block = lines[4 + na:4 + na + nb]
        if len(block) != nb:
            out["raw_error"] = "bond block truncated (%d != %d)" % (len(block), nb)
        for bl in block:
            try:
                if int(bl[6:9]) == 4:
                    arom += 1
            except ValueError:
                pass
        out["aromatic_flag_bonds"] = arom
    except Exception as exc:  # noqa: BLE001 - 诊断脚本需要报告任意解析异常
        out["raw_error"] = "%s: %s" % (type(exc).__name__, exc)
    return out


def ligand_geometry(path):
    """RDKit 视角的配体信息。sanitize=True 失败时退回 sanitize=False, 保持可读性结论真实。"""
    info = {
        "readable": False, "atom_count": None, "heavy_atom_count": None,
        "bond_count": None, "num_conformers": None, "is3d": None,
        "num_fragments": None, "orphan_atoms": None, "duplicate_coords": None,
        "nonfinite_coords": None, "max_abs_coord": None, "centroid": None,
        "coords": None, "error": None,
    }
    mol = None
    try:
        supp = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
        mol = next(iter(supp), None)
    except Exception as exc:  # noqa: BLE001
        info["error"] = "sanitize=True: %s" % exc
    if mol is None:
        try:
            supp = Chem.SDMolSupplier(path, removeHs=False, sanitize=False)
            mol = next(iter(supp), None)
            if mol is not None:
                info["error"] = "needs sanitize=False"
        except Exception as exc:  # noqa: BLE001
            info["error"] = "sanitize=False: %s" % exc
    if mol is None:
        return info

    info["readable"] = True
    n = mol.GetNumAtoms()
    info["atom_count"] = n
    info["bond_count"] = mol.GetNumBonds()
    info["heavy_atom_count"] = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 1)
    info["num_fragments"] = len(Chem.GetMolFrags(mol)) if n else 0
    info["num_conformers"] = mol.GetNumConformers()
    conf = mol.GetConformer() if mol.GetNumConformers() > 0 else None
    info["is3d"] = bool(conf.Is3D()) if conf is not None else None

    bonded = set()
    for b in mol.GetBonds():
        bonded.add(b.GetBeginAtomIdx())
        bonded.add(b.GetEndAtomIdx())
    info["orphan_atoms"] = sum(1 for k in range(n) if k not in bonded)

    if conf is not None:
        coords = np.array(conf.GetPositions(), dtype=np.float64)
        heavy = np.array([a.GetAtomicNum() != 1 for a in mol.GetAtoms()])
        info["nonfinite_coords"] = int((~np.isfinite(coords)).sum())
        fin = coords[np.isfinite(coords)]
        info["max_abs_coord"] = float(np.abs(fin).max()) if fin.size else None
        info["coords"] = coords[heavy]
        if heavy.sum() > 0:
            info["centroid"] = coords[heavy].mean(axis=0)
        seen = Counter(map(tuple, np.round(coords[heavy], 4).tolist()))
        info["duplicate_coords"] = int(sum(c - 1 for c in seen.values() if c > 1))
    return info


def pocket_geometry(path, use_repo_parser=True):
    """口袋几何: 优先用训练同款 utils.protein_ligand.PDBProtein, 失败退回手工解析。"""
    info = {
        "parse_ok": False, "parse_error": None, "element": None, "pos": None,
        "residue_count": None, "chain_count": None,
    }
    if use_repo_parser:
        try:
            from utils.protein_ligand import PDBProtein
            protein = PDBProtein(path)
            info["element"] = np.asarray(protein.element, dtype=int)
            info["pos"] = np.asarray(protein.pos, dtype=np.float64)
            info["residue_count"] = len(protein.residues)
            info["chain_count"] = len({r.get("chain", "") for r in protein.residues})
            info["parse_ok"] = True
            return info
        except Exception as exc:  # noqa: BLE001
            info["parse_error"] = "%s: %s" % (type(exc).__name__, exc)
    # 退回: 手工解析 ATOM 行 (与 PDBProtein 口径一致, 忽略 HETATM)
    try:
        pos, elem, residues, chains = [], [], set(), set()
        with open(path, "r", encoding="ascii", errors="ignore") as fh:
            for line in fh:
                if line[0:6].strip() != "ATOM":
                    continue
                pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                symb = line[76:78].strip()
                if not symb:
                    symb = line[13:14]
                elem.append(Chem.GetPeriodicTable().GetAtomicNumber(symb.capitalize()))
                chains.add(line[21:22].strip())
                residues.add((line[21:22], line[22:26], line[26:27], line[17:20]))
        info["element"] = np.asarray(elem, dtype=int)
        info["pos"] = np.asarray(pos, dtype=np.float64).reshape(-1, 3)
        info["residue_count"] = len(residues)
        info["chain_count"] = len(chains)
        info["parse_ok"] = info["pos"].shape[0] > 0
        if info["parse_error"] is None:
            info["parse_error"] = "fallback parser used"
    except Exception as exc:  # noqa: BLE001
        if info["parse_error"] is None:
            info["parse_error"] = "%s: %s" % (type(exc).__name__, exc)
    return info


def min_distances(lig_pos, pocket_pos, chunk=2000):
    """返回 (每个配体原子到口袋的最近距离, 最近的距离)。分块避免大口袋爆内存。"""
    if lig_pos is None or pocket_pos is None or len(lig_pos) == 0 or len(pocket_pos) == 0:
        return None, None
    best = np.full(len(lig_pos), np.inf)
    for start in range(0, len(pocket_pos), chunk):
        block = pocket_pos[start:start + chunk]
        d = np.linalg.norm(lig_pos[:, None, :] - block[None, :, :], axis=-1)
        best = np.minimum(best, d.min(axis=1))
    return best, float(best.min())


# --------------------------------------------------------------------------- #
# mask 覆盖度: 复现 utils/transforms.py 的 LigandMixedMask 逻辑
# --------------------------------------------------------------------------- #
def mask_num_masked(num_atoms, ratio, min_num_masked, min_num_unmasked):
    """与 LigandRandomMask / LigandBFSMask 中完全一致的 num_masked 计算。"""
    num_masked = int(num_atoms * ratio)
    if num_masked < min_num_masked:
        num_masked = min_num_masked
    if (num_atoms - num_masked) < min_num_unmasked:
        num_masked = num_atoms - min_num_unmasked
    return num_masked


def p_empty_context_theory(num_atoms, bfs_len, mask_cfg, n_draws=20000, rng=None):
    """蒙特卡洛估计"配体上下文为空"的概率。

    关键: LigandRandomMask 取 idx[:num_masked] / idx[num_masked:]；
          LigandBFSMask 取 perm[-num_masked:] / perm[:-num_masked]。
    因此 context 为空 <=> num_masked >= (N 或 bfs_len)。
    """
    rng = rng or np.random.default_rng(12345)
    mn, mx = mask_cfg["min_ratio"], mask_cfg["max_ratio"]
    mn_masked, mn_unmasked = mask_cfg["min_num_masked"], mask_cfg["min_num_unmasked"]
    ratios = np.clip(rng.uniform(mn, mx, size=n_draws), 0.0, 1.0)
    rand_empty = 0
    bfs_empty = 0
    for r in ratios:
        nm = mask_num_masked(num_atoms, r, mn_masked, mn_unmasked)
        if nm >= num_atoms:
            rand_empty += 1
        if nm >= bfs_len:
            bfs_empty += 1
    p_rand = rand_empty / float(n_draws)
    p_bfs = bfs_empty / float(n_draws)
    p_rand_w = mask_cfg.get("p_random", 0.15)
    p_bfs_w = mask_cfg.get("p_bfs", 0.6) + mask_cfg.get("p_invbfs", 0.25)
    total_w = p_rand_w + p_bfs_w
    if total_w <= 0:
        return 0.0
    return (p_rand_w * p_rand + p_bfs_w * p_bfs) / total_w


class _MiniLigandData(object):
    """给真实 LigandMixedMask 用的最小数据容器 (不需要 torch_geometric)。"""

    def __init__(self, element, pos, bond_index, bond_type):
        import torch

        n = int(element.shape[0])
        self.ligand_element = torch.as_tensor(element, dtype=torch.long)
        self.ligand_pos = torch.as_tensor(pos, dtype=torch.float32)
        self.ligand_bond_index = torch.as_tensor(bond_index, dtype=torch.long)
        self.ligand_bond_type = torch.as_tensor(bond_type, dtype=torch.long)
        self.ligand_atom_feature_full = torch.zeros(n, 13, dtype=torch.long)
        deg = torch.zeros(n, dtype=torch.long)
        for i in self.ligand_bond_index[0].tolist():
            deg[i] += 1
        self.ligand_num_neighbors = deg
        nb = {}
        for k, j in enumerate(self.ligand_bond_index[1].tolist()):
            nb.setdefault(self.ligand_bond_index[0, k].item(), []).append(j)
        self.ligand_nbh_list = nb


def simulate_mask_reality(samples, mask_cfg, draws, seed=2021):
    """用仓库真实的 LigandMixedMask 类验证上面的理论概率。返回 {idx: measured_p}"""
    try:
        from utils.transforms import LigandMixedMask, LigandBFSMask
    except Exception as exc:  # noqa: BLE001
        return {}, "import utils.transforms failed: %s" % exc
    masker = LigandMixedMask(
        mask_cfg["min_ratio"], mask_cfg["max_ratio"],
        mask_cfg["min_num_masked"], mask_cfg["min_num_unmasked"],
        mask_cfg["p_random"], mask_cfg["p_bfs"], mask_cfg["p_invbfs"],
    )
    random.seed(seed)
    np.random.seed(seed)
    measured = {}
    for idx, element, pos, bond_index, bond_type in samples:
        empty = 0
        for _ in range(draws):
            data = _MiniLigandData(element, pos, bond_index, bond_type)
            try:
                data = masker(data)
            except Exception:  # noqa: BLE001
                empty += 1
                continue
            if data.context_idx.nelement() == 0:
                empty += 1
        measured[idx] = empty / float(draws)
    return measured, None


def bfs_perm_length(element, pos, bond_index, bond_type):
    """用真实 LigandBFSMask.get_bfs_perm 求 BFS 访问序列长度 (连通分子应等于原子数)。"""
    try:
        from utils.transforms import LigandBFSMask
    except Exception:  # noqa: BLE001
        return None
    data = _MiniLigandData(element, pos, bond_index, bond_type)
    random.seed(0)
    try:
        perm, _, _ = LigandBFSMask.get_bfs_perm(data.ligand_nbh_list)
        return int(len(perm))
    except Exception:  # noqa: BLE001
        return None


# --------------------------------------------------------------------------- #
# 过滤规则
# --------------------------------------------------------------------------- #
def evaluate_flags(row):
    """返回该样本命中的异常标签列表。阈值集中在此处, 便于复核。"""
    flags = []
    if not row["pocket_exists"]:
        flags.append("POCKET_MISSING")
    if not row["ligand_exists"]:
        flags.append("LIGAND_MISSING")
    if row["pocket_exists"] and not row["pocket_parse_ok"]:
        flags.append("POCKET_PARSE_ERROR")
    if row["pocket_unknown_resnames"]:
        flags.append("POCKET_NONSTD_RESIDUE")
    if row["pocket_atom_count"] is not None:
        if row["pocket_atom_count"] == 0:
            flags.append("POCKET_EMPTY")
        elif row["pocket_atom_count"] < 40:
            flags.append("POCKET_TOO_FEW_ATOMS")
    if row["pocket_residue_count"] is not None and row["pocket_residue_count"] < 8:
        flags.append("POCKET_TOO_FEW_RESIDUES")
    if row["pocket_nonfinite_coords"]:
        flags.append("POCKET_NONFINITE_COORDS")
    if row["pocket_max_abs_coord"] is not None and row["pocket_max_abs_coord"] > 1e4:
        flags.append("POCKET_COORD_TOO_LARGE")

    if row["ligand_exists"]:
        if not row["ligand_readable"]:
            flags.append("LIGAND_UNREADABLE")
        elif not row["ligand_parse_ok"]:
            flags.append("LIGAND_PARSE_ERROR")
        if row["ligand_num_conformers"] is not None and row["ligand_num_conformers"] < 1:
            flags.append("LIGAND_NO_CONFORMER")
        if row["ligand_is3d"] is False:
            flags.append("LIGAND_NOT_3D")
        if row["ligand_nonfinite_coords"]:
            flags.append("LIGAND_NONFINITE_COORDS")
        if row["ligand_max_abs_coord"] is not None and row["ligand_max_abs_coord"] > 1e4:
            flags.append("LIGAND_COORD_TOO_LARGE")
        if row["ligand_heavy_atom_count"] is not None and row["ligand_heavy_atom_count"] < 5:
            flags.append("LIGAND_TOO_FEW_HEAVY_ATOMS")
        if row["ligand_num_fragments"] is not None and row["ligand_num_fragments"] > 1:
            flags.append("LIGAND_FRAGMENTED")
        if row["ligand_orphan_atoms"]:
            flags.append("LIGAND_ORPHAN_ATOMS")
        if row["ligand_aromatic_flag_bonds"]:
            flags.append("LIGAND_AROMATIC_FLAG_BONDS")
        if row["ligand_duplicate_coords"]:
            flags.append("LIGAND_DUPLICATE_COORDS")
        if row["in_name2id"] is False:
            flags.append("NOT_IN_LMDB_NAME2ID")

        frac2 = row["frac_heavy_inside_pocket_bbox_pad2"]
        if frac2 is not None:
            if frac2 < 0.5:
                flags.append("LIGAND_OUTSIDE_POCKET")
            elif frac2 < 0.9:
                flags.append("LIGAND_PARTLY_OUTSIDE_POCKET")
        mind = row["ligand_min_dist_to_pocket"]
        if mind is not None:
            if mind > 4.0:
                flags.append("LIGAND_FAR_FROM_POCKET")
            if mind < 0.5:
                # 与真实实验位姿 (最小重原子距离 2.5~3.3 A) 相比, <0.5 A 是原子级穿模
                flags.append("LIGAND_CLASH_POCKET")
            elif mind < 1.5:
                flags.append("LIGAND_CLOSE_CONTACT")
        f4 = row["frac_heavy_within_4A_of_pocket"]
        if f4 is not None:
            if f4 <= 0.0:
                flags.append("LIGAND_NO_ATOM_WITHIN_4A")
            if f4 < 0.5:
                flags.append("LIGAND_LOW_POCKET_CONTACT")
        cdist = row["pocket_ligand_centroid_distance"]
        if cdist is not None and cdist > 15.0:
            flags.append("CENTROID_MISMATCH")
    return flags


# 建议丢弃的标签 (DROP_* 组成推荐过滤规则)
DROP_FLAGS = [
    "POCKET_MISSING", "LIGAND_MISSING", "POCKET_PARSE_ERROR", "POCKET_NONSTD_RESIDUE",
    "POCKET_EMPTY", "POCKET_TOO_FEW_ATOMS", "POCKET_TOO_FEW_RESIDUES",
    "POCKET_NONFINITE_COORDS", "POCKET_COORD_TOO_LARGE",
    "LIGAND_UNREADABLE", "LIGAND_PARSE_ERROR", "LIGAND_NO_CONFORMER", "LIGAND_NOT_3D",
    "LIGAND_NONFINITE_COORDS", "LIGAND_COORD_TOO_LARGE",
    "LIGAND_TOO_FEW_HEAVY_ATOMS", "LIGAND_FRAGMENTED", "LIGAND_ORPHAN_ATOMS",
    "NOT_IN_LMDB_NAME2ID",
    "LIGAND_CLASH_POCKET", "LIGAND_LOW_POCKET_CONTACT",
    "LIGAND_OUTSIDE_POCKET", "LIGAND_FAR_FROM_POCKET", "LIGAND_NO_ATOM_WITHIN_4A",
    "CENTROID_MISMATCH",
]
REVIEW_FLAGS = ["LIGAND_CLOSE_CONTACT", "LIGAND_PARTLY_OUTSIDE_POCKET",
                "LIGAND_DUPLICATE_COORDS"]
REMEDIATE_FLAGS = ["LIGAND_AROMATIC_FLAG_BONDS"]

# 三档可执行规则组合 (名称 -> 标签集合)
RULE_COMBOS = [
    ("conservative (R2+R3+R4: outside/low-contact/far)",
     ["LIGAND_OUTSIDE_POCKET", "LIGAND_LOW_POCKET_CONTACT",
      "LIGAND_FAR_FROM_POCKET", "LIGAND_NO_ATOM_WITHIN_4A"]),
    ("recommended  (R1+R2+R3+R4: + gross clash <0.5A)",
     ["LIGAND_CLASH_POCKET", "LIGAND_OUTSIDE_POCKET", "LIGAND_LOW_POCKET_CONTACT",
      "LIGAND_FAR_FROM_POCKET", "LIGAND_NO_ATOM_WITHIN_4A"]),
    ("aggressive   (recommended + close contact 0.5~1.5A)",
     ["LIGAND_CLASH_POCKET", "LIGAND_CLOSE_CONTACT", "LIGAND_OUTSIDE_POCKET",
      "LIGAND_LOW_POCKET_CONTACT", "LIGAND_FAR_FROM_POCKET",
      "LIGAND_NO_ATOM_WITHIN_4A"]),
]


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #
def load_split(path):
    import torch
    try:
        obj = torch.load(path, weights_only=False)
    except TypeError:
        obj = torch.load(path)
    out = {}
    for key, items in obj.items():
        for pocket, ligand in items:
            out[(pocket, ligand)] = key
    return out


def load_name2id(data_dir, dataset_name):
    import torch
    path = os.path.join(os.path.dirname(os.path.abspath(data_dir)),
                        dataset_name + "_name2id.pt")
    if not os.path.exists(path):
        return None, path
    try:
        obj = torch.load(path, weights_only=False)
    except TypeError:
        obj = torch.load(path)
    return set(obj.keys()), path


def default_mask_cfg():
    """默认取 configs/train_gpcr_v3.yml 的 transform.mask, 保证诊断与真实训练同口径。"""
    cfg = {"min_ratio": 0.0, "max_ratio": 1.1, "min_num_masked": 1,
           "min_num_unmasked": 0, "p_random": 0.15, "p_bfs": 0.6, "p_invbfs": 0.25}
    path = os.path.join(CONFIGS, "train_gpcr_v3.yml")
    try:
        import yaml
        with open(path, "r", encoding="utf-8") as fh:
            loaded = yaml.safe_load(fh)
        cfg.update(loaded["train"]["transform"]["mask"])
        return cfg, path
    except Exception:  # noqa: BLE001
        return cfg, None


def main():
    ap = argparse.ArgumentParser(description="gpcr_v3_merged read-only diagnosis")
    ap.add_argument("--data", default=os.path.join(DATA, "gpcr_v3_merged"))
    ap.add_argument("--out", default=None, help="CSV 路径, 默认 <data>/diagnosis.csv")
    ap.add_argument("--summary-out", default=None,
                    help="文本汇总路径, 默认 <data>/diagnosis_summary.txt")
    ap.add_argument("--mask-draws", type=int, default=400,
                    help="用真实 LigandMixedMask 验证的抽样次数/样本")
    ap.add_argument("--mask-samples", type=int, default=150,
                    help="参与真实 mask 模拟的样本数 (0 表示跳过)")
    ap.add_argument("--mask-draws-theory", type=int, default=20000)
    ap.add_argument("--top", type=int, default=20, help="打印最严重样本条数")
    ap.add_argument("--emit-drop-list", default=None,
                    help="可选: 把推荐规则的 DROP/KEEP 清单写到该路径 (默认不写, 保持只读)")
    args = ap.parse_args()

    data_dir = os.path.abspath(args.data)
    out_csv = args.out or os.path.join(data_dir, "diagnosis.csv")
    out_summary = args.summary_out or os.path.join(data_dir, "diagnosis_summary.txt")
    dataset_name = os.path.basename(data_dir)

    index_path = os.path.join(data_dir, "index.pkl")
    split_path = os.path.join(data_dir, "split_by_name.pt")
    if not os.path.exists(index_path):
        raise SystemExit("index.pkl not found: %s" % index_path)
    with open(index_path, "rb") as fh:
        index = pickle.load(fh)
    split_map = load_split(split_path) if os.path.exists(split_path) else {}
    name2id, name2id_path = load_name2id(data_dir, dataset_name)
    mask_cfg, mask_cfg_path = default_mask_cfg()

    report = []

    def emit(text=""):
        print(text)
        report.append(text)

    emit("=" * 78)
    emit("Pocket2Mol dataset diagnosis (read-only)")
    emit("=" * 78)
    emit("data dir        : %s" % data_dir)
    emit("index entries   : %d" % len(index))
    emit("split names     : %d (train=%d, test=%d)"
         % (len(split_map),
            sum(1 for v in split_map.values() if v == "train"),
            sum(1 for v in split_map.values() if v == "test")))
    emit("name2id file    : %s" % (name2id_path if name2id is not None
                                   else "NOT FOUND"))
    if name2id is not None:
        emit("name2id entries : %d  -> %d of %d index entries are usable"
             % (len(name2id), sum(1 for e in index if (e[0], e[1]) in name2id), len(index)))
    emit("mask config     : %s" % mask_cfg)
    emit("mask config from: %s" % (mask_cfg_path or "built-in default"))
    emit("")

    rng_mask = np.random.default_rng(12345)
    rows = []
    mask_sim_pool = []
    flag_counter = Counter()
    flag_by_source = {"legacy": Counter(), "rcsb": Counter()}
    per_source = Counter()
    per_split_source = Counter()
    theory_pool = []

    for idx, entry in enumerate(index):
        pocket_fn, ligand_fn = entry[0], entry[1]
        pocket_path = os.path.join(data_dir, pocket_fn)
        ligand_path = os.path.join(data_dir, ligand_fn)
        source = "legacy" if is_source_legacy(pocket_fn) else "rcsb"
        split = split_map.get((pocket_fn, ligand_fn), "not_in_split")
        per_source[source] += 1
        per_split_source[(split, source)] += 1

        row = OrderedDict((k, None) for k in CSV_COLUMNS)
        row["idx"] = idx
        row["split"] = split
        row["source"] = source
        row["pocket_file"] = pocket_fn
        row["ligand_file"] = ligand_fn
        row["pocket_exists"] = os.path.exists(pocket_path)
        row["ligand_exists"] = os.path.exists(ligand_path)
        row["in_name2id"] = (None if name2id is None
                             else ((pocket_fn, ligand_fn) in name2id))
        row["mask_p_empty_context_theory"] = None
        row["mask_p_empty_context_sim"] = None

        # ---- 口袋 ----
        pocket_pos = None
        if row["pocket_exists"]:
            scan = scan_pdb_lines(pocket_path) or {"hetatm": 0, "hydro": 0,
                                                   "unknown": [], "atom_lines": 0}
            pg = pocket_geometry(pocket_path)
            row["pocket_parse_ok"] = pg["parse_ok"]
            row["pocket_parse_error"] = pg["parse_error"] or ""
            row["pocket_hetatm_count"] = scan["hetatm"]
            row["pocket_hydrogen_count"] = scan["hydro"]
            row["pocket_unknown_resnames"] = ";".join(scan["unknown"])
            pocket_pos = pg["pos"]
            row["pocket_residue_count"] = pg["residue_count"]
            row["pocket_chain_count"] = pg["chain_count"]
            if pocket_pos is not None and len(pocket_pos):
                row["pocket_atom_count"] = int(len(pocket_pos))
                p_heavy = pocket_pos[pg["element"] != 1] if pg["element"] is not None else pocket_pos
                if len(p_heavy) == 0:
                    p_heavy = pocket_pos
                row["pocket_nonfinite_coords"] = int((~np.isfinite(pocket_pos)).sum())
                fin = pocket_pos[np.isfinite(pocket_pos)]
                row["pocket_max_abs_coord"] = float(np.abs(fin).max()) if fin.size else None
                span = p_heavy.max(axis=0) - p_heavy.min(axis=0)
                row["pocket_bbox_dx"], row["pocket_bbox_dy"], row["pocket_bbox_dz"] = (
                    float(span[0]), float(span[1]), float(span[2]))
                row["pocket_centroid_x"], row["pocket_centroid_y"], row["pocket_centroid_z"] = (
                    float(p_heavy.mean(axis=0)[0]), float(p_heavy.mean(axis=0)[1]),
                    float(p_heavy.mean(axis=0)[2]))
                row["_pocket_heavy"] = p_heavy
            else:
                row["pocket_atom_count"] = 0

        # ---- 配体 ----
        lig = {}
        if row["ligand_exists"]:
            raw = read_sdf_raw_text(ligand_path)
            row["ligand_declared_atoms"] = raw["declared_atoms"]
            row["ligand_declared_bonds"] = raw["declared_bonds"]
            row["ligand_aromatic_flag_bonds"] = raw["aromatic_flag_bonds"]
            lig = ligand_geometry(ligand_path)
            row["ligand_readable"] = lig["readable"]
            row["ligand_atom_count"] = lig["atom_count"]
            row["ligand_heavy_atom_count"] = lig["heavy_atom_count"]
            row["ligand_bond_count"] = lig["bond_count"]
            row["ligand_num_conformers"] = lig["num_conformers"]
            row["ligand_is3d"] = lig["is3d"]
            row["ligand_num_fragments"] = lig["num_fragments"]
            row["ligand_orphan_atoms"] = lig["orphan_atoms"]
            row["ligand_duplicate_coords"] = lig["duplicate_coords"]
            row["ligand_nonfinite_coords"] = lig["nonfinite_coords"]
            row["ligand_max_abs_coord"] = lig["max_abs_coord"]
            if lig["centroid"] is not None:
                row["ligand_centroid_x"] = float(lig["centroid"][0])
                row["ligand_centroid_y"] = float(lig["centroid"][1])
                row["ligand_centroid_z"] = float(lig["centroid"][2])
            # 训练同款解析 (元素/坐标/键), 用于复现 mask 行为
            try:
                from utils.protein_ligand import parse_sdf_file
                pd = parse_sdf_file(ligand_path)
                row["ligand_parse_ok"] = True
                row["_train_element"] = np.asarray(pd["element"], dtype=int)
                row["_train_pos"] = np.asarray(pd["pos"], dtype=np.float64)
                row["_train_bond_index"] = np.asarray(pd["bond_index"], dtype=int)
                row["_train_bond_type"] = np.asarray(pd["bond_type"], dtype=int)
            except Exception as exc:  # noqa: BLE001
                row["ligand_parse_ok"] = False
                row["ligand_parse_error"] = "%s: %s" % (type(exc).__name__, exc)

        # ---- 几何关系 ----
        if lig.get("coords") is not None and row.get("_pocket_heavy") is not None:
            lig_heavy = np.asarray(lig["coords"], dtype=np.float64)
            pock = row["_pocket_heavy"]
            if len(lig_heavy) and len(pock):
                lo, hi = pock.min(axis=0) - 2.0, pock.max(axis=0) + 2.0
                inside2 = np.all((lig_heavy >= lo) & (lig_heavy <= hi), axis=1)
                row["frac_heavy_inside_pocket_bbox_pad2"] = float(inside2.mean())
                lo0, hi0 = pock.min(axis=0), pock.max(axis=0)
                inside0 = np.all((lig_heavy >= lo0) & (lig_heavy <= hi0), axis=1)
                row["frac_heavy_inside_pocket_bbox_pad0"] = float(inside0.mean())
                best, mind = min_distances(lig_heavy, pock)
                row["ligand_min_dist_to_pocket"] = mind
                if best is not None:
                    row["frac_heavy_within_4A_of_pocket"] = float((best <= 4.0).mean())
                if lig.get("centroid") is not None:
                    pc = pock.mean(axis=0)
                    row["pocket_ligand_centroid_distance"] = float(
                        np.linalg.norm(np.asarray(lig["centroid"]) - pc))

        # ---- mask 覆盖度理论值 ----
        if row.get("_train_element") is not None:
            n_atoms = int(row["_train_element"].shape[0])
            blen = bfs_perm_length(row["_train_element"], row["_train_pos"],
                                   row["_train_bond_index"], row["_train_bond_type"])
            if blen is None:
                blen = n_atoms
            p_empty = p_empty_context_theory(n_atoms, blen, mask_cfg,
                                             n_draws=args.mask_draws_theory,
                                             rng=rng_mask)
            row["mask_p_empty_context_theory"] = float(p_empty)
            theory_pool.append((idx, n_atoms, blen, float(p_empty)))
            if len(mask_sim_pool) < args.mask_samples:
                mask_sim_pool.append((idx, row["_train_element"], row["_train_pos"],
                                      row["_train_bond_index"], row["_train_bond_type"]))

        flags = evaluate_flags(row)
        row["flags"] = ";".join(flags)
        for f in flags:
            flag_counter[f] += 1
            flag_by_source[source][f] += 1
        rows.append(row)

    # ---- 真实 mask 类验证 ----
    sim_note = None
    if args.mask_samples > 0 and mask_sim_pool:
        measured, sim_note = simulate_mask_reality(mask_sim_pool, mask_cfg,
                                                   args.mask_draws)
        for row in rows:
            if row["idx"] in measured:
                row["mask_p_empty_context_sim"] = measured[row["idx"]]

    # ---- 写 CSV (中文路径安全: 绝对路径 + utf-8-sig) ----
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.writer(fh)
        writer.writerow(CSV_COLUMNS)
        for row in rows:
            out = []
            for col in CSV_COLUMNS:
                val = row.get(col)
                if col in ("idx",):
                    out.append(int(val))
                elif isinstance(val, bool):
                    out.append("True" if val else "False")
                elif isinstance(val, float):
                    out.append(fnum(val))
                elif val is None:
                    out.append("")
                else:
                    out.append(val)
            writer.writerow(out)

    # ---- 汇总 ----
    n = len(rows)
    anomalous_rows = [r for r in rows if r["flags"]]
    drop_rows = [r for r in rows if any(f in DROP_FLAGS for f in r["flags"].split(";") if f)]
    emit("-" * 78)
    emit("1) 异常样本总览")
    emit("-" * 78)
    emit("样本总数                : %d" % n)
    emit("至少命中一个异常标签    : %d (%.2f%%)" % (len(anomalous_rows),
                                                 100.0 * len(anomalous_rows) / max(n, 1)))
    emit("命中推荐丢弃规则        : %d (%.2f%%)" % (len(drop_rows),
                                                 100.0 * len(drop_rows) / max(n, 1)))
    emit("")
    emit("按标签统计 (标签: 命中数 | legacy | rcsb):")
    for flag, cnt in flag_counter.most_common():
        emit("  %-32s %5d | %5d | %5d"
             % (flag, cnt, flag_by_source["legacy"][flag], flag_by_source["rcsb"][flag]))
    emit("")

    emit("-" * 78)
    emit("2) 来源分布")
    emit("-" * 78)
    for src in ("legacy", "rcsb"):
        sub = [r for r in rows if r["source"] == src]
        bad = [r for r in sub if r["flags"]]
        emit("来源 %-6s : %4d 样本, 异常 %3d (%.2f%%), 训练集 %d / 验证集 %d"
             % (src, len(sub), len(bad),
                100.0 * len(bad) / max(len(sub), 1),
                sum(1 for r in sub if r["split"] == "train"),
                sum(1 for r in sub if r["split"] == "test")))
    for flag, _ in flag_counter.most_common(8):
        sub = {"legacy": flag_by_source["legacy"][flag], "rcsb": flag_by_source["rcsb"][flag]}
        tot = sub["legacy"] + sub["rcsb"]
        if tot:
            emit("  标签 %-30s legacy %5d (%.1f%%) | rcsb %5d (%.1f%%)"
                 % (flag, sub["legacy"], 100.0 * sub["legacy"] / max(per_source["legacy"], 1),
                    sub["rcsb"], 100.0 * sub["rcsb"] / max(per_source["rcsb"], 1)))
    emit("")

    emit("-" * 78)
    emit("3) mask 覆盖度 -> 非有限损失的直接来源")
    emit("-" * 78)
    if theory_pool:
        ps = [t[3] for t in theory_pool]
        emit("参与计算的样本          : %d" % len(ps))
        emit("预测 P(配体上下文为空)  : min %.4f | mean %.4f | max %.4f"
             % (min(ps), sum(ps) / len(ps), max(ps)))
        uniq = len({round(p, 3) for p in ps})
        emit("不同取值的概率个数      : %d (为 1 表示与样本无关)" % uniq)
        if sim_note:
            emit("真实 LigandMixedMask 验证: 跳过 (%s)" % sim_note)
        else:
            sims = [r["mask_p_empty_context_sim"] for r in anomalous_rows + rows
                    if r["mask_p_empty_context_sim"] is not None]
            sims = [s for s in sims if s is not None]
            if sims:
                emit("真实 LigandMixedMask 实测: min %.4f | mean %.4f | max %.4f (n=%d, draws=%d)"
                     % (min(sims), sum(sims) / len(sims), max(sims), len(sims), args.mask_draws))
        emit("理论闭式 (max_ratio>1 且 min_num_unmasked=0): P = (max_ratio-1)/(max_ratio-min_ratio) = %.4f"
             % (((mask_cfg["max_ratio"] - 1.0) / (mask_cfg["max_ratio"] - mask_cfg["min_ratio"]))
                if (mask_cfg["max_ratio"] > 1.0 and mask_cfg["min_num_unmasked"] == 0
                    and mask_cfg["max_ratio"] > mask_cfg["min_ratio"]) else 0.0))
    emit("")

    emit("-" * 78)
    emit("4) 建议过滤规则 (DROP) 与影响面")
    emit("-" * 78)
    def impact(tag):
        hit = [r for r in rows if tag in (r["flags"].split(";") if r["flags"] else [])]
        return (len(hit),
                sum(1 for r in hit if r["split"] == "train"),
                sum(1 for r in hit if r["split"] == "test"))
    for tag in DROP_FLAGS:
        c, tr, te = impact(tag)
        if c:
            emit("  %-32s 命中 %4d | train %4d | test %3d  -> 丢弃" % (tag, c, tr, te))
    for tag in REVIEW_FLAGS:
        c, tr, te = impact(tag)
        if c:
            emit("  %-32s 命中 %4d | train %4d | test %3d  -> 人工复核" % (tag, c, tr, te))
    for tag in REMEDIATE_FLAGS:
        c, tr, te = impact(tag)
        if c:
            emit("  %-32s 命中 %4d | train %4d | test %3d  -> 先 kekulize 修复, 不必丢弃"
                 % (tag, c, tr, te))
    emit("")
    emit("推荐规则并集 = 上面所有 DROP 标签。合计丢弃 %d / %d (%.2f%%), "
         "其中 train %d, test %d。"
         % (len(drop_rows), n, 100.0 * len(drop_rows) / max(n, 1),
            sum(1 for r in drop_rows if r["split"] == "train"),
            sum(1 for r in drop_rows if r["split"] == "test")))
    emit("")
    emit("三档规则组合的量化影响:")
    combo_sets = {}
    for label, tags in RULE_COMBOS:
        sel = [r for r in rows
               if any(t in (r["flags"].split(";") if r["flags"] else []) for t in tags)]
        combo_sets[label] = set(r["idx"] for r in sel)
        keep = [r for r in rows if r["idx"] not in combo_sets[label]]
        emit("  %-52s drop %4d (%.2f%%) | train -%4d | test -%3d | keep train %4d/test %3d"
             % (label, len(sel), 100.0 * len(sel) / max(n, 1),
                sum(1 for r in sel if r["split"] == "train"),
                sum(1 for r in sel if r["split"] == "test"),
                sum(1 for r in keep if r["split"] == "train"),
                sum(1 for r in keep if r["split"] == "test")))
    rec_label = RULE_COMBOS[1][0]
    rec_ids = combo_sets.get(rec_label, set())
    if rec_ids:
        comp = Counter(r["ligand_file"].split("_")[-1].replace(".sdf", "")
                       for r in rows if r["idx"] in rec_ids)
        emit("")
        emit("推荐丢弃集合的配体类型 top10: %s" % (comp.most_common(10),))
        emit("推荐丢弃集合涉及不同口袋文件数: %d" % len({r["pocket_file"] for r in rows
                                                       if r["idx"] in rec_ids}))
        all_comp = Counter(r["ligand_file"].split("_")[-1].replace(".sdf", "") for r in rows)
        emit("全数据集该配体类型占比 top5 : %s" % (all_comp.most_common(5),))
        emit("(提示: 若某类配体在丢弃集合中占比远高于全集占比, 就是系统性污染源)")
    if args.emit_drop_list:
        with open(args.emit_drop_list, "w", encoding="utf-8") as fh:
            for r in rows:
                mark = "DROP" if r["idx"] in rec_ids else "KEEP"
                fh.write("%s\t%d\t%s\t%s\t%s\n"
                         % (mark, r["idx"], r["split"], r["source"], r["flags"] or "-"))
        emit("")
        emit("过滤清单: %s (DROP=%d)" % (args.emit_drop_list, len(rec_ids)))

    emit("")
    emit("-" * 78)
    emit("5) 最可疑样本 (按异常标签数排序)")
    emit("-" * 78)
    ranked = sorted(rows, key=lambda r: -len([f for f in r["flags"].split(";") if f]))
    shown = 0
    for r in ranked:
        if not r["flags"]:
            break
        emit("  #%-5d %-28s %-34s %s"
             % (r["idx"], r["pocket_file"], r["ligand_file"], r["flags"]))
        shown += 1
        if shown >= args.top:
            break
    if shown == 0:
        emit("  (无)")

    emit("")
    emit("CSV   : %s" % out_csv)
    emit("SUMMARY: %s" % out_summary)
    with open(out_summary, "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
