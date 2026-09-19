# -*- coding: utf-8 -*-
"""任务9 模型自洽性: kNN 域外检测 (候选 vs 已知化学空间的距离)
输出: compounds_tagged.csv 追加 extrapolation 列 (inner/transition/outer) + 交叉表报告
"""
import os, sys, csv, pickle
sys.path.insert(0, r"D:\MMModel\Pocket2Mol")
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import AllChem, DataStructs

K = 5
KNOWN = r"D:\MMModel\已知药物库"
MTV2 = r"D:\MMModel\Pocket2Mol\data\gpcr_multitarget_v2"
LIB = r"D:\MMModel\化合物库"

def fps_from_smis(smis):
    out = []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m:
            out.append((s, AllChem.GetMorganFingerprint(m, 2)))
    return out

# ---------- 参照库: 训练集配体 + 各靶点已知集 ----------
ref = []
idx = pickle.load(open(os.path.join(MTV2, "index.pkl"), "rb"))
for pocket, lig, _, _ in idx:
    m = Chem.MolFromMolFile(os.path.join(MTV2, lig), removeHs=True, sanitize=True)
    if m is not None:
        ref.append(Chem.MolToSmiles(m))
print("训练集配体:", len(ref))
for fn in ["A2A_known_set.txt", "beta-2_adrenergic_receptor_known_set.txt",
           "dopamine_D3_receptor_known_set.txt", "5-HT2B_receptor_known_set.txt"]:
    p = os.path.join(KNOWN, fn)
    if os.path.exists(p):
        smis = [s for s in open(p, encoding="utf-8").read().splitlines() if s.strip()]
        ref += smis
        print("  +", fn, len(smis))
ref = sorted(set(ref))
ref_fps = fps_from_smis(ref)
print("参照库去重后:", len(ref_fps))

# ---------- 全库计算 ----------
rows = list(csv.DictReader(open(os.path.join(LIB, "compounds_tagged.csv"), encoding="utf-8")))
dists, labels = [], []
for r in rows:
    m = Chem.MolFromSmiles(r["smiles"])
    if m is None:
        dists.append(np.nan); labels.append("unknown"); continue
    fp = AllChem.GetMorganFingerprint(m, 2)
    sims = DataStructs.BulkTanimotoSimilarity(fp, [f for _, f in ref_fps])
    sims.sort(reverse=True)
    sim_k = float(np.mean(sims[:K]))
    dists.append(1.0 - sim_k)
a = np.array([d for d in dists if not np.isnan(d)])
p50, p75 = np.percentile(a, 50), np.percentile(a, 75)
print("kNN 距离分布: P25=%.3f P50=%.3f P75=%.3f max=%.3f" % (
    np.percentile(a, 25), p50, p75, a.max()))
for d in dists:
    if np.isnan(d):
        labels.append("unknown")
    elif d > p75:
        labels.append("outer")
    elif d > p50:
        labels.append("transition")
    else:
        labels.append("inner")
print("标签分布:", {x: labels.count(x) for x in set(labels)})

# ---------- 写列 ----------
with open(os.path.join(LIB, "compounds_tagged.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    head = list(rows[0].keys()) + ["knn_dist", "extrapolation"]
    w.writerow(head)
    for r, d, lb in zip(rows, dists, labels):
        w.writerow([r[k] for k in rows[0].keys()] + [round(d, 3) if d == d else "", lb])
print("compounds_tagged.csv 已更新")

# ---------- 交叉表: 外推 × 新颖性 ----------
nov = {r["smiles"]: r["category"] for r in csv.DictReader(open(os.path.join(LIB, "A2A" and "top50_novelty.csv"), encoding="utf-8-sig"))} if os.path.exists(os.path.join(LIB, "top50_novelty.csv")) else {}
cross = {}
for r, lb in zip(rows, labels):
    if r["smiles"] in nov:
        c = nov[r["smiles"]]
        cross[(lb, c)] = cross.get((lb, c), 0) + 1
print("Top50 交叉表 (外推标签 × 新颖性):", cross if cross else "novelty 文件在库目录缺失")
# Top50 外推分布
top50 = [lb for r, lb in zip(rows, labels) if r["smiles"] in nov]
print("Top50 外推分布:", {x: top50.count(x) for x in set(top50)})
