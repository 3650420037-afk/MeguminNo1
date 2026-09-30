# -*- coding: utf-8 -*-
"""v10 库新颖性复核：对 737 已知活性集的最大 Tanimoto + 骨架重合统计。

用于 PPT/交付件的新颖性数字取证（可复现）。
输出: outputs/night_pipeline_v10/novelty_v10_summary.csv / .txt
"""
import csv, os, statistics, sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")

KNOWN = {
    "A2A":   os.path.join(ROOT, "data", "known_drugs", "A2A_known_set.txt"),
    "D3":    os.path.join(ROOT, "data", "known_drugs", "dopamine_D3_receptor_known_set.txt"),
    "5HT2B": os.path.join(ROOT, "data", "known_drugs", "5-HT2B_receptor_known_set.txt"),
}
LIB = {
    "A2A":   [os.path.join(ROOT, "results", "compounds_A2A.csv")],
    "D3":    [os.path.join(ROOT, "results", "D3", "compounds.csv")],
    "5HT2B": [os.path.join(ROOT, "results", "5HT2B", "compounds.csv")],
}
TOP = {
    "A2A":   [os.path.join(ROOT, "results", "top_candidates.csv")],
    "D3":    [os.path.join(ROOT, "results", "D3", "top_candidates.csv")],
    "5HT2B": [os.path.join(ROOT, "results", "5HT2B", "top_candidates.csv")],
}


def read_smiles(path, col=None):
    out = []
    with open(path, encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f)
        col = col or ("smiles" if "smiles" in (rd.fieldnames or []) else rd.fieldnames[0])
        for r in rd:
            s = (r.get(col) or "").strip()
            if s:
                out.append(s)
    return out


def first_existing(paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return None


def scaf_of(m):
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        return ""


rows_out = []
for tgt, kpath in KNOWN.items():
    if not os.path.exists(kpath):
        print("[%s] 已知集缺失: %s" % (tgt, kpath)); continue
    ksmi = read_smiles(kpath, col=None)
    kfps, kscafs = [], set()
    for s in ksmi:
        m = Chem.MolFromSmiles(s)
        if not m:
            continue
        kfps.append(AllChem.GetMorganFingerprint(m, 2))
        sc = scaf_of(m)
        if sc:
            kscafs.add(sc)
    lpath = first_existing(LIB[tgt])
    lsmi = read_smiles(lpath)
    tpath = first_existing(TOP[tgt])
    tsmi = set(read_smiles(tpath)) if tpath else set()

    sims, hits_scaf, top_sims = [], 0, []
    for s in lsmi:
        m = Chem.MolFromSmiles(s)
        if not m:
            continue
        fp = AllChem.GetMorganFingerprint(m, 2)
        mx = max(DataStructs.BulkTanimotoSimilarity(fp, kfps)) if kfps else 0.0
        sims.append(mx)
        if scaf_of(m) in kscafs:
            hits_scaf += 1
        if s in tsmi:
            top_sims.append(mx)
    sims.sort()
    n = len(sims)
    line = ("[%s] 已知集 %d 分子/%d 骨架 | 库 %d 分子 | 最大相似 中位 %.3f 最大 %.3f | "
            "<0.4 %.1f%% | >=0.85 %.1f%% | 骨架重合 %d (%.1f%%)") % (
        tgt, len(kfps), len(kscafs), n, statistics.median(sims), sims[-1],
        100.0 * sum(1 for x in sims if x < 0.4) / n,
        100.0 * sum(1 for x in sims if x >= 0.85) / n,
        hits_scaf, 100.0 * hits_scaf / n)
    if top_sims:
        line += " | Top%d 最大 %.3f 全部<0.56:%s" % (
            len(top_sims), max(top_sims), all(x < 0.56 for x in top_sims))
    print(line)
    rows_out.append({
        "target": tgt, "known_n": len(kfps), "known_scaffolds": len(kscafs),
        "lib_n": n, "sim_median": round(statistics.median(sims), 3), "sim_max": round(sims[-1], 3),
        "pct_lt_0.4": round(100.0 * sum(1 for x in sims if x < 0.4) / n, 1),
        "pct_ge_0.85": round(100.0 * sum(1 for x in sims if x >= 0.85) / n, 1),
        "scaffold_overlap": hits_scaf,
        "top_n": len(top_sims), "top_sim_max": round(max(top_sims), 3) if top_sims else "",
        "top_all_lt_0.56": all(x < 0.56 for x in top_sims) if top_sims else "",
        "lib_path": os.path.relpath(lpath, ROOT).replace("\\", "/"),
    })

out = os.path.join(ROOT, "outputs", "night_pipeline_v10", "novelty_v10_summary.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
    w.writeheader()
    w.writerows(rows_out)
print("WROTE", out)
