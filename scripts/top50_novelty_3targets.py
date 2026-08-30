# -*- coding: utf-8 -*-
"""三靶点各自: 拉 ChEMBL 已知活性集 -> Top50 综合排序 -> 新颖性验证."""
import csv, json, os, sys, time, urllib.request
sys.path.insert(0, r"D:\MMModel\Pocket2Mol")
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, DataStructs

LIB = r"D:\MMModel\化合物库"
OUT = r"D:\MMModel\Pocket2Mol\outputs\docking_%s_library\summary.csv"

TARGETS = {
    "b2ar":  ("CHEMBL279",  "beta-2 adrenergic receptor"),
    "d3":    ("CHEMBL234",  "dopamine D3 receptor"),
    "5ht2b": ("CHEMBL1837", "5-HT2B receptor"),
}

def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2)

def known_set(chembl_id, name):
    ids = set()
    for off in range(0, 600, 200):
        acts = get("https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=%s&pchembl_value__gte=8&limit=200&offset=%d" % (chembl_id, off))["activities"]
        for a in acts:
            ids.add(a["molecule_chembl_id"])
        if len(acts) < 200: break
    smis = []
    idl = sorted(ids)
    for i in range(0, len(idl), 80):
        ms = get("https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id__in=%s&limit=100" % ",".join(idl[i:i+80]))["molecules"]
        for m in ms:
            sm = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if sm: smis.append(sm)
        time.sleep(0.3)
    with open(r"D:\MMModel\已知药物库\ChEMBL_已知药物_100.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["canonical_smiles"].strip(): smis.append(row["canonical_smiles"].strip())
    uniq = sorted(set(smis))
    p = os.path.join(r"D:\MMModel\已知药物库", "%s_known_set.txt" % name.replace(" ", "_"))
    open(p, "w", encoding="utf-8").write("\n".join(uniq))
    return uniq, p

def fps_and_scafs(smis):
    fps, scafs = [], set()
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if not m: continue
        fps.append(AllChem.GetMorganFingerprint(m, 2))
        try: scafs.add(MurckoScaffold.MurckoScaffoldSmiles(mol=m))
        except Exception: pass
    return fps, scafs

for key, (cid, name) in TARGETS.items():
    print("=" * 70)
    print("[%s] %s (%s)" % (key, name, cid))
    known, kpath = known_set(cid, name)
    kfps, kscafs = fps_and_scafs(known)
    print("  已知集: %d 分子, %d 骨架 -> %s" % (len(kfps), len(kscafs), kpath))

    # Top50 综合排序
    rows = [r for r in csv.DictReader(open(OUT % key, encoding="utf-8")) if r["status"] == "ok"]
    lib = {r["smiles"]: r for r in csv.DictReader(open(os.path.join(LIB, key, "compounds.csv"), encoding="utf-8"))}
    vmin, vmax = min(float(r["vina_score"]) for r in rows), max(float(r["vina_score"]) for r in rows)
    scored = []
    for r in rows:
        lib_r = lib.get(r["smiles"])
        if not lib_r: continue
        v, qed = float(r["vina_score"]), float(lib_r["qed"])
        vn = (v - vmax) / (vmin - vmax)
        scored.append((0.6 * vn + 0.4 * qed, v, qed, float(lib_r["sa"]), float(lib_r["mw"]), lib_r["murcko_scaffold"], r["smiles"]))
    scored.sort(reverse=True)
    top = scored[:50]
    with open(os.path.join(LIB, key, "top_candidates.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["rank", "综合分", "vina_kcal", "qed", "sa", "mw", "murcko", "smiles"])
        for i, (t, v, q, s, m, sc, smi) in enumerate(top, 1):
            w.writerow([i, round(t, 4), v, q, s, round(m, 1), sc, smi])
    print("  Top50: vina %.2f 至 %.2f" % (min(x[1] for x in top), max(x[1] for x in top)))

    # 新颖性
    cats = {"KNOWN_SCAFFOLD": 0, "NEAR_KNOWN": 0, "mid_similarity": 0, "HIGH_NOVEL": 0}
    det = []
    for t, v, q, s, m, sc, smi in top:
        mm = Chem.MolFromSmiles(smi)
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mm)
        mx = max(DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprint(mm, 2), kfps)) if kfps else 0.0
        cat = "KNOWN_SCAFFOLD" if scaf in kscafs else ("NEAR_KNOWN" if mx >= 0.85 else ("mid_similarity" if mx >= 0.4 else "HIGH_NOVEL"))
        cats[cat] += 1
        det.append((cat, mx, v, q, smi))
    print("  新颖性: %s" % cats)
    with open(os.path.join(LIB, key, "top50_novelty.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["rank", "category", "max_tanimoto_vs_known", "vina_kcal", "qed", "smiles"])
        for i, (cat, mx, v, q, smi) in enumerate(det, 1):
            w.writerow([i, cat, round(mx, 3), v, q, smi])
print("ALL DONE")
