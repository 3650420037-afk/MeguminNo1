# -*- coding: utf-8 -*-
"""v10 库 Top50（按库内综合分排序）新颖性复核：最大 Tanimoto、是否全部 <0.4 / <0.56。"""
import csv, os, statistics, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")

KNOWN = {"A2A": "A2A_known_set.txt", "D3": "dopamine_D3_receptor_known_set.txt", "5HT2B": "5-HT2B_receptor_known_set.txt"}
TOP = {"A2A": "results/top_candidates.csv", "D3": "results/D3/top_candidates.csv", "5HT2B": "results/5HT2B/top_candidates.csv"}

for tgt, fn in KNOWN.items():
    kfps, kscaf = [], set()
    for s in open(os.path.join(ROOT, "data", "known_drugs", fn), encoding="utf-8"):
        s = s.strip()
        if not s:
            continue
        m = Chem.MolFromSmiles(s)
        if not m:
            continue
        kfps.append(AllChem.GetMorganFingerprint(m, 2))
        try:
            kscaf.add(MurckoScaffold.MurckoScaffoldSmiles(mol=m))
        except Exception:
            pass
    rows = list(csv.DictReader(open(os.path.join(ROOT, TOP[tgt]), encoding="utf-8-sig", newline="")))
    res = []
    for r in rows[:50]:
        m = Chem.MolFromSmiles(r["smiles"])
        mx = max(DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprint(m, 2), kfps))
        try:
            sc = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
        except Exception:
            sc = ""
        res.append((int(r["rank"]), float(r["vina_kcal"]), float(r["qed"]), round(mx, 3), sc in kscaf))
    sims = sorted(x[3] for x in res)
    print("[%s] Top50: 最大相似 %.3f 中位 %.3f | <0.4 数 %d/50 | <0.56 数 %d/50 | 骨架重合 %d | vina %.2f~%.2f" % (
        tgt, sims[-1], statistics.median(sims),
        sum(1 for x in sims if x < 0.4), sum(1 for x in sims if x < 0.56),
        sum(1 for x in res if x[4]), min(x[1] for x in res), max(x[1] for x in res)))
    worst = max(res, key=lambda x: x[3])
    print("      相似度最高者 rank=%d vina=%.2f qed=%.2f sim=%.3f 骨架重合=%s" % worst)
