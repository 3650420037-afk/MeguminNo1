# -*- coding: utf-8 -*-
"""重生成 Top50 候选榜与新颖性验证 (A2A + 三靶点)

背景: 此前 top_candidates.csv 由内联脚本生成, 且受 docking summary.csv 的 BOM 丢样
污染(每靶点少 1 个候选参与排序)。summary 已用 norm_docking_summary.py 修复, 本脚本
按修复后的数据重建榜单, 并把生成逻辑固化为可复跑脚本。

输入: outputs/docking_<target>_library/summary.csv + results/<target>/compounds.csv
      + data/known_drugs/<target> 的 known_set (新颖性对照)
输出: results/<target>/top_candidates.csv, top50_novelty.csv
"""
import os, csv, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import KNOWN_DRUGS, OUTPUTS, RESULTS, src_on_path
src_on_path()
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, DataStructs

LIB = RESULTS
OUTROOT = OUTPUTS
KNOWN = KNOWN_DRUGS

# target -> (docking 目录, 化合物库目录, 已知集文件)
TARGETS = {
    "A2A":   ("docking_a2a_library",  LIB,                     "A2A_known_set.txt"),
    "b2ar":  ("docking_b2ar_library", os.path.join(LIB, "b2ar"),  "beta-2_adrenergic_receptor_known_set.txt"),
    "d3":    ("docking_d3_library",   os.path.join(LIB, "d3"),    "dopamine_D3_receptor_known_set.txt"),
    "5ht2b": ("docking_5ht2b_library", os.path.join(LIB, "5ht2b"), "5-HT2B_receptor_known_set.txt"),
}


def known_fps(scaf_file):
    p = os.path.join(KNOWN, scaf_file)
    smis = [s for s in open(p, encoding="utf-8").read().splitlines() if s.strip()] if os.path.exists(p) else []
    scafs, fps = set(), []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        fps.append(AllChem.GetMorganFingerprint(m, 2))
        try:
            scafs.add(MurckoScaffold.MurckoScaffoldSmiles(mol=m))
        except Exception:
            pass
    return scafs, fps


def run(target):
    dock_dir, lib_dir, known_file = TARGETS[target]
    summ = os.path.join(OUTROOT, dock_dir, "summary.csv")
    comp = os.path.join(lib_dir, "compounds.csv")
    if not (os.path.exists(summ) and os.path.exists(comp)):
        print("  [%s] 缺输入, 跳过" % target); return
    rows = [r for r in csv.DictReader(open(summ, encoding="utf-8-sig"))
            if r["status"] == "ok" and r["vina_score"]]
    lib = {r["smiles"]: r for r in csv.DictReader(open(comp, encoding="utf-8-sig"))}
    matched = [r for r in rows if r["smiles"] in lib]
    print("  [%s] 对接成功 %d | 与库匹配 %d (%.1f%%)" % (target, len(rows), len(matched), 100.0*len(matched)/max(1,len(rows))))
    v = np.array([float(r["vina_score"]) for r in matched])
    vmin, vmax = v.min(), v.max()
    scored = []
    for r in matched:
        q = float(lib[r["smiles"]]["qed"])
        vn = (float(r["vina_score"]) - vmax) / (vmin - vmax) if vmax > vmin else 0.0
        scored.append((0.6 * vn + 0.4 * q, float(r["vina_score"]), q,
                       float(lib[r["smiles"]]["sa"]), float(lib[r["smiles"]]["mw"]),
                       lib[r["smiles"]]["murcko_scaffold"], r["smiles"]))
    scored.sort(reverse=True)
    top = scored[:50]
    with open(os.path.join(lib_dir, "top_candidates.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["rank", "综合分", "vina_kcal", "qed", "sa", "mw", "murcko", "smiles"])
        for i, (t, vv, q, s, m, sc, smi) in enumerate(top, 1):
            w.writerow([i, round(t, 4), vv, q, s, round(m, 1), sc, smi])
    print("    Top50: vina %.2f ~ %.2f | Top1 综合 %.4f" % (min(x[1] for x in top), max(x[1] for x in top), top[0][0]))

    # 新颖性
    kscaf, kfps = known_fps(known_file)
    cats = {"KNOWN_SCAFFOLD": 0, "NEAR_KNOWN": 0, "mid_similarity": 0, "HIGH_NOVEL": 0}
    det = []
    for t, vv, q, s, m, sc, smi in top:
        mm = Chem.MolFromSmiles(smi)
        if mm is None:
            continue
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mm)
        except Exception:
            scaf = ""
        mx = max(DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprint(mm, 2), kfps)) if kfps else 0.0
        cat = "KNOWN_SCAFFOLD" if scaf in kscaf else ("NEAR_KNOWN" if mx >= 0.85 else ("mid_similarity" if mx >= 0.4 else "HIGH_NOVEL"))
        cats[cat] += 1
        det.append((cat, mx, vv, q, smi))
    with open(os.path.join(lib_dir, "top50_novelty.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["rank", "category", "max_tanimoto_vs_known", "vina_kcal", "qed", "smiles"])
        for i, (cat, mx, vv, q, smi) in enumerate(det, 1):
            w.writerow([i, cat, round(mx, 3), vv, q, smi])
    print("    新颖性: %s" % cats)


if __name__ == "__main__":
    targets = sys.argv[1:] or list(TARGETS.keys())
    for t in targets:
        run(t)
