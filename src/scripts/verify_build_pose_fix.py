# -*- coding: utf-8 -*-
"""离线端到端验证: 配体是否还用"重建构象 + 质心平移"(会穿模), 还是改用实验坐标。

不依赖网络: 使用本地已有的真实 RCSB mmCIF + 从 CCD 块抽出的权威 SMILES
(scripts/ligand_smiles_ccd_offline.json, 取自各条目自带的 _chem_comp_atom /
_chem_comp_bond, 与 fetch_gpcr_structures.py 用的 rcsb_chem_comp_descriptor 同源)。
指标与 scripts/diagnose_dataset.py 完全一致(同一 PDBProtein 解析、同一 min_distances、
同一判据 LIGAND_CLASH_POCKET = 最近重原子距离 < 0.5 A)。

对每个 (结构, 配体) 计算三条:
  EXP  实验位姿(参考真值)
  NEW  build_from_experimental(): 实验坐标 + SMILES 校正键级  <- 本次修复
  OLD  build_conformer(): SMILES 生成构象 + 质心平移          <- 修复前行为

用法(默认 cif 目录为本地靶点结构目录, 可按需 --src):
    python scripts/verify_build_pose_fix.py

配套的流水线级冒烟(验证 --on-mismatch 两种策略与产出报告):
    python scripts/build_gpcr_v3.py --src <含 cif + pairs.csv 的目录> \
        --out data/_v3_smoke_out --min-samples 2
    python scripts/diagnose_dataset.py --data data/_v3_smoke_out --mask-samples 0
"""
import os, sys, json, warnings, argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import DATA, TARGETS, src_on_path
src_on_path()

warnings.simplefilter("ignore")

import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

from scripts.build_gpcr_v3 import (parse_cif_once, build_from_experimental,
                                   base_index, heavy_atom_count)
from scripts.build_gpcr_dataset import build_conformer
from scripts.diagnose_dataset import pocket_geometry, ligand_geometry, min_distances
from utils.protein_ligand import PDBProtein

RADIUS = 12.0
SEED = 2021

# 判据 (与 diagnose_dataset.py 完全一致)
CLASH_MAX = 0.5
CLOSE_MAX = 1.5
FAR_MIN = 4.0
BBOX_MIN = 0.5
WITHIN4_MIN = 0.5


def flags_for(mind, frac_bbox2, frac4):
    fl = []
    if frac_bbox2 is not None and frac_bbox2 < BBOX_MIN:
        fl.append("OUTSIDE_POCKET")
    if mind is not None:
        if mind > FAR_MIN:
            fl.append("FAR_FROM_POCKET")
        if mind < CLASH_MAX:
            fl.append("CLASH_POCKET")
        elif mind < CLOSE_MAX:
            fl.append("CLOSE_CONTACT")
    if frac4 is not None and frac4 < WITHIN4_MIN:
        fl.append("LOW_POCKET_CONTACT")
    return fl


def measure(lig_pos, pocket):
    """返回 (mind, frac_bbox_pad2, frac_within_4A)"""
    if lig_pos is None or pocket is None or len(lig_pos) == 0 or len(pocket) == 0:
        return None, None, None
    lo, hi = pocket.min(axis=0) - 2.0, pocket.max(axis=0) + 2.0
    inside = np.all((lig_pos >= lo) & (lig_pos <= hi), axis=1)
    best, mind = min_distances(lig_pos, pocket)
    f4 = float((best <= 4.0).mean()) if best is not None else None
    return mind, float(inside.mean()), f4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=TARGETS)
    ap.add_argument("--smiles-json",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "ligand_smiles_ccd_offline.json"))
    ap.add_argument("--work", default=os.path.join(DATA, "_pose_verify_tmp"))
    ap.add_argument("--out", default=os.path.join(DATA, "_pose_verify.json"))
    ap.add_argument("--radius", type=float, default=RADIUS)
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    smap = json.load(open(args.smiles_json, encoding="utf-8"))

    cases = [(os.path.join(args.src, "3PBL.cif"), "3PBL", "ETQ"),
             (os.path.join(args.src, "4IB4.cif"), "4IB4", "ERM"),
             (os.path.join(args.src, "5UEN.cif"), "5UEN", "DU1")]

    rows = []
    for cif, pid, cid in cases:
        print("=" * 74)
        print("%s / %s" % (pid, cid), flush=True)
        if not os.path.exists(cif):
            print("  缺 cif, 跳过"); continue
        smiles = smap.get(cid)
        if not smiles:
            print("  缺 SMILES, 跳过"); continue

        pdb = os.path.join(args.work, "%s.pdb" % pid)
        cents, lig_atoms = parse_cif_once(cif, pdb)
        if cid not in cents:
            print("  cif 里没有配体 %s, 跳过" % cid); continue
        protein = PDBProtein(pdb)
        residues = protein.query_residues_radius(
            cents[cid], args.radius, criterion="center_of_mass")
        print("  口袋残基=%d  (半径 %.1f A)" % (len(residues), args.radius))
        if len(residues) < 8:
            print("  口袋残基过少, 跳过"); continue

        pk = os.path.join(args.work, "pocket_%s_%s.pdb" % (pid, cid))
        with open(pk, "w", encoding="ascii") as h:
            h.write(protein.residues_to_pdb_block(residues, name="GPCR_POCKET"))
        pg = pocket_geometry(pk)
        pocket = pg["pos"] if pg["parse_ok"] else None
        print("  口袋重原子=%s  parse_ok=%s" % (
            None if pocket is None else len(pocket), pg["parse_ok"]))

        atoms = lig_atoms[cid]
        exp_pos = np.array([a[2] for a in atoms], dtype=np.float64)

        # ---- EXP ----
        m_exp, fb_exp, f4_exp = measure(exp_pos, pocket)

        # ---- NEW ----
        mol_new, new_status = None, "?"
        try:
            mol_new, new_status = build_from_experimental(atoms, smiles)
        except Exception as e:
            print("  NEW 异常: %s" % e)
            new_status = "exception"
        new_pos, new_dev, new_sdf = None, None, None
        if mol_new is not None:
            new_sdf = os.path.join(args.work, "lig_new_%s_%s.sdf" % (pid, cid))
            w = Chem.SDWriter(new_sdf); w.write(mol_new); w.close()
            lg = ligand_geometry(new_sdf)
            new_pos = lg["coords"]
            if new_pos is not None and len(new_pos) == len(exp_pos):
                new_dev = float(np.abs(new_pos - exp_pos).max())
        m_new, fb_new, f4_new = measure(new_pos, pocket)
        print("  NEW 状态=%s (沉积重原子=%d, SMILES 重原子=%d)"
              % (new_status, len(atoms), heavy_atom_count(smiles)))

        # ---- OLD ----
        old_pos, old_dev, old_sdf = None, None, None
        try:
            mol_old = build_conformer(smiles, cents[cid], SEED + base_index(pid, 0))
            old_sdf = os.path.join(args.work, "lig_old_%s_%s.sdf" % (pid, cid))
            w = Chem.SDWriter(old_sdf); w.write(mol_old); w.close()
            lg2 = ligand_geometry(old_sdf)
            old_pos = lg2["coords"]
            if old_pos is not None and len(old_pos) == len(exp_pos):
                old_dev = float(np.abs(old_pos - exp_pos).max())
        except Exception as e:
            print("  OLD 异常: %s" % e)
        m_old, fb_old, f4_old = measure(old_pos, pocket)

        for tag, mind, fb, f4, dev, nsdf in (
                ("EXP", m_exp, fb_exp, f4_exp, 0.0, None),
                ("NEW", m_new, fb_new, f4_new, new_dev, new_sdf),
                ("OLD", m_old, fb_old, f4_old, old_dev, old_sdf)):
            print("  %-3s  最近重原子距离=%-7s bbox内占比=%-6s 4A内占比=%-6s "
                  "坐标偏差=%-8s %s"
                  % (tag,
                     "None" if mind is None else "%.3f" % mind,
                     "None" if fb is None else "%.2f" % fb,
                     "None" if f4 is None else "%.2f" % f4,
                     "None" if dev is None else "%.4f" % dev,
                     ",".join(flags_for(mind, fb, f4)) or "-"))
            rows.append({"pdb": pid, "comp": cid, "kind": tag,
                         "min_dist": mind, "frac_bbox": fb, "frac_within4": f4,
                         "dev_from_exp": dev, "status": new_status if tag == "NEW" else
                         ("ok" if mind is not None else "unusable"),
                         "flags": flags_for(mind, fb, f4), "sdf": nsdf})

    json.dump(rows, open(args.out, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)

    print("=" * 74)
    print("汇总 (判据与 diagnose_dataset.py 一致: CLASH=最近重原子距离<%.1f, "
          "CLOSE=%.1f~%.1f)" % (CLASH_MAX, CLASH_MAX, CLOSE_MAX))
    print("%-6s %-4s %-4s %-9s %-9s %-18s %s"
          % ("结构", "配体", "路径", "最近距离", "坐标偏差", "状态", "标签"))
    for r in rows:
        print("%-6s %-4s %-4s %-9s %-9s %-18s %s" % (
            r["pdb"], r["comp"], r["kind"],
            "None" if r["min_dist"] is None else "%.3f" % r["min_dist"],
            "None" if r["dev_from_exp"] is None else "%.4f" % r["dev_from_exp"],
            r.get("status") or "-",
            ",".join(r["flags"]) or "-"))

    print("-" * 74)
    for kind in ("EXP", "NEW", "OLD"):
        sub = [r for r in rows if r["kind"] == kind]
        if not sub:
            continue
        usable = [r for r in sub if r["min_dist"] is not None]
        n_clash = sum(1 for r in usable if "CLASH_POCKET" in r["flags"])
        n_close = sum(1 for r in usable if "CLOSE_CONTACT" in r["flags"])
        print("%-4s 样本=%d 可用=%d 弃用=%d | 可用中: 穿模=%d (%.1f%%) 紧接触=%d (%.1f%%)"
              % (kind, len(sub), len(usable), len(sub) - len(usable),
                 n_clash, 100.0 * n_clash / max(1, len(usable)),
                 n_close, 100.0 * n_close / max(1, len(usable))))
    # 修复的核心主张: NEW 的坐标 == 实验坐标
    new_ok = [r for r in rows if r["kind"] == "NEW" and r["dev_from_exp"] is not None]
    if new_ok:
        worst = max(r["dev_from_exp"] for r in new_ok)
        print("NEW 与实验位姿的最大坐标偏差: %.4f A (n=%d) -> %s"
              % (worst, len(new_ok),
                 "逐原子复现实验坐标" if worst < 0.01 else "存在偏差, 需复查"))
    print("明细: %s" % args.out)


if __name__ == "__main__":
    main()
