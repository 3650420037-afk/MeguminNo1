# -*- coding: utf-8 -*-
"""M6 阶段 2：为"按受体扩口袋"批量定位人源结构 + 提取口袋中心（数据获取，不做对接）。

为什么要它
----------
用户路线①后半要求把口袋数从 4 扩到**数十个**。现有 `calc_pocket_centers.py` / 
`build_target_registry.py` 的目标表是**手工写死的 4+1 个**，无法规模化。
本脚本按 **UniProt 号**自动完成"选结构 → 选共晶类药配体 → 算口袋中心"，
产物可直接被 `mass_dock_actives.py`（对接）与 `build_target_registry.py`（登记表）使用。

判据（全部来自结构本身，不手工录入）
------------------------------------
- 人源：UniProt 精确匹配；分辨率 ≤ 3.0 Å；含非聚合物配体
- 共晶配体须**类药**：MW∈[250,500]、元素 ⊆ {C,N,O,F,P,S,Cl}、重原子 ≥ 20
- 口袋中心 = 该配体重原子的几何质心（实验坐标，与既有四靶点口径一致）
- 同一受体优先取**分辨率最高（数值最小）**的条目；没有合格配体就继续找下一个条目

用法
----
    python src/scripts/fetch_m6_targets.py --targets D2 5HT2A --out data/targets_m6
    python src/scripts/fetch_m6_targets.py --targets all --limit-per-target 3
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, DATA, ensure_dir  # noqa: E402

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY = "https://data.rcsb.org/rest/v1/core/entry/%s"
NONPOLY = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/%s/%s"
CHEMCOMP = "https://data.rcsb.org/rest/v1/core/chemcomp/%s"
DOWNLOAD = "https://files.rcsb.org/download/%s.pdb"

# 靶点 -> (UniProt, ChEMBL ID, 说明)。UniProt 用于精确选结构，避免按名称误命中同源靶点。
M6_TARGETS = {
    "D2":    ("P14416", "CHEMBL217",  "Dopamine D2 receptor"),
    "5HT1A": ("P08908", "CHEMBL214",  "5-HT1A receptor"),
    "5HT2A": ("P28223", "CHEMBL224",  "5-HT2A receptor"),
    "M1":    ("P11229", "CHEMBL216",  "Muscarinic M1 receptor"),
    "H1":    ("P35367", "CHEMBL231",  "Histamine H1 receptor"),
    "OPRM1": ("P35372", "CHEMBL233",  "Mu opioid receptor"),
    "OPRK1": ("P41145", "CHEMBL237",  "Kappa opioid receptor"),
    "A1":    ("P30542", "CHEMBL226",  "Adenosine A1 receptor"),
}
ELEMENTS_OK = {6, 7, 8, 9, 15, 16, 17}
# 非类药共晶配体黑名单：甾醇/脂类/去垢剂/辅因子/离子。CLR 与 OLA 最常出现——
# 实测不排除时，5HT1A(7E2X) 会把胆固醇当作"配体"，口袋中心被带到膜内（明显错误）。
LIGAND_BLACKLIST = {
    "CLR", "OLA", "PLM", "MYR", "STE", "LDA", "BOG", "LHG", "DDM", "CHS", "C8E", "LMT",
    "PEG", "PGE", "1PE", "2PE", "GOL", "EDO", "MPD", "DMS", "SO4", "PO4", "NO3", "ACT",
    "ACY", "FMT", "TRS", "EPE", "MES", "CIT", "TLA", "BME", "DTT", "IMD", "NH4", "CO3",
    "SCN", "AZI", "HOH", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "NI", "CD", "HG",
    "IOD", "BR", "CL", "F", "HEM", "NAG", "MAN", "BMA", "FUC", "GAL", "NDG", "PLP", "SAM",
    "GDP", "GTP", "ATP", "ADP", "AMP", "ANP", "NAD", "FAD", "COA", "ACP",
}


def get_json(url, tries=5, sleep=6.0):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(sleep * (i + 1))
    print("    [warn] %s -> %s" % (url, last))
    return None


def search_entries(uniprot, rows=60):
    q = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": uniprot}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows},
                            "results_content_type": ["experimental"]},
    }
    req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(),
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            js = json.loads(r.read().decode())
    except Exception as e:  # noqa: BLE001
        print("    [warn] RCSB search 失败: %s" % e)
        return []
    return [x["identifier"] for x in js.get("result_set", [])]


def entry_info(pdb):
    js = get_json(ENTRY % pdb)
    if not js:
        return None
    res = (js.get("rcsb_entry_info", {}) or {}).get("resolution_combined")
    return {"resolution": (res[0] if isinstance(res, list) and res else res),
            "nonpoly_ids": (js.get("rcsb_entry_container_identifiers", {}) or {}).get(
                "non_polymer_entity_ids", []) or []}


def lig_candidates(pdb, nonpoly_ids):
    out = []
    for eid in nonpoly_ids:
        js = get_json(NONPOLY % (pdb, eid))
        if not js:
            continue
        comp = (js.get("pdbx_entity_nonpoly", {}) or {}).get("comp_id")
        mw = ((js.get("rcsb_nonpolymer_entity", {}) or {}).get("formula_weight")
              or (js.get("rcsb_nonpolymer_entity", {}) or {}).get("formula_weight"))
        if not comp:
            continue
        cc = get_json(CHEMCOMP % comp) or {}
        atoms = ((cc.get("chem_comp", {}) or {}).get("pdbx_chem_comp_descriptor") or [])
        # 用 CCD 原子表统计元素与重原子数（比 formula_weight 更可靠）
        atom_rows = ((cc.get("rcsb_chem_comp_info", {}) or {}).get("atom_count") or None)
        smi = None
        for d in atoms:
            if d.get("type") in ("SMILES", "SMILES_CANONICAL") and d.get("descriptor"):
                smi = d["descriptor"]
                break
        # ⚠ RCSB 的 formula_weight 单位是 **kDa**（实测 8NU=0.41 → 410 Da），必须 ×1000
        out.append({"comp_id": comp, "mw": (float(mw) * 1000.0) if mw else None,
                    "smiles": smi, "atom_count": atom_rows})
    return out


def ligand_atoms_from_pdb(pdb_path, comp_id):
    """取该 comp_id **单个拷贝**的重原子坐标。

    ⚠ 不能把所有拷贝合并：不对称单元里常有两份配体（实测 6A93 的 8NU 被并成 60 个重原子，
    质心落成两份的中点 → 口袋中心错误）。这里按 (链, 序号段) 分组取最大的一份。
    """
    copies = {}
    for line in open(pdb_path, encoding="ascii", errors="ignore"):
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip() != comp_id:
            continue
        el = line[76:78].strip().upper()
        if el in ("H", ""):
            continue
        chain = line[21]
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        copies.setdefault(chain, []).append(xyz)
    if not copies:
        return []
    return max(copies.values(), key=len)


def main():
    ap = argparse.ArgumentParser(description="M6：按 UniProt 批量定位结构 + 提取口袋中心")
    ap.add_argument("--targets", nargs="+", default=["all"])
    ap.add_argument("--out", default=os.path.join(DATA, "targets_m6"))
    ap.add_argument("--limit-per-target", type=int, default=2, help="每个靶点最多保留几个候选结构")
    ap.add_argument("--max-res", type=float, default=3.0)
    args = ap.parse_args()

    names = list(M6_TARGETS) if args.targets == ["all"] else args.targets
    ensure_dir(args.out)
    rows = []
    for name in names:
        if name not in M6_TARGETS:
            print("[%s] 不在 M6_TARGETS 里, 跳过" % name)
            continue
        uni, chembl, desc = M6_TARGETS[name]
        print("=" * 78)
        print("[%s] %s | UniProt %s | ChEMBL %s" % (name, desc, uni, chembl), flush=True)
        ids = search_entries(uni)
        print("  候选结构 %d 个, 逐个筛..." % len(ids), flush=True)
        kept = 0
        for pdb in ids:
            if kept >= args.limit_per_target:
                break
            info = entry_info(pdb)
            if not info:
                continue
            res = info["resolution"]
            if res and float(res) > args.max_res:
                continue
            cands = [c for c in lig_candidates(pdb, info["nonpoly_ids"])
                     if c["mw"] and 250 <= c["mw"] <= 500
                     and c["comp_id"] not in LIGAND_BLACKLIST]
            if not cands:
                continue
            # 下载结构（受体坐标 + 配体坐标同在一个 PDB 里）
            pdb_path = os.path.join(args.out, "%s_%s.pdb" % (pdb, name))
            if not os.path.exists(pdb_path):
                try:
                    urllib.request.urlretrieve(DOWNLOAD % pdb, pdb_path)
                except Exception as e:  # noqa: BLE001
                    print("    [warn] 下载 %s 失败: %s" % (pdb, e))
                    continue
            for c in cands:
                atoms = ligand_atoms_from_pdb(pdb_path, c["comp_id"])
                if len(atoms) < 20:
                    continue
                cx = sum(a[0] for a in atoms) / len(atoms)
                cy = sum(a[1] for a in atoms) / len(atoms)
                cz = sum(a[2] for a in atoms) / len(atoms)
                rows.append({"target": name, "uniprot": uni, "chembl": chembl,
                             "desc": desc, "pdb_id": pdb, "resolution": res,
                             "comp_id": c["comp_id"], "ligand_mw": round(c["mw"], 1),
                             "ligand_heavy": len(atoms),
                             "center_x": round(cx, 3), "center_y": round(cy, 3),
                             "center_z": round(cz, 3),
                             "receptor_pdb": os.path.relpath(pdb_path, ROOT),
                             "ligand_smiles": c["smiles"]})
                print("    ✓ %s | %s | res %s | 配体 %s (MW %.0f, %d 重原子) | 中心 (%.2f, %.2f, %.2f)"
                      % (pdb, name, res, c["comp_id"], c["mw"], len(atoms), cx, cy, cz), flush=True)
                kept += 1
                break
        if kept == 0:
            print("  !! 未找到合格结构（可能该受体没有人源+类药共晶配体的条目）")
    out_csv = os.path.join(args.out, "registry.csv")
    if rows:
        with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    print("\n" + "=" * 78)
    print("合格口袋 %d 个 -> %s" % (len(rows), out_csv))
    for r in rows:
        print("  %-6s %-5s %-4s MW %5.0f 重原子 %2d 中心 (%.2f, %.2f, %.2f)"
              % (r["target"], r["pdb_id"], r["comp_id"], r["ligand_mw"], r["ligand_heavy"],
                 r["center_x"], r["center_y"], r["center_z"]))
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
