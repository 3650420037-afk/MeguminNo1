# -*- coding: utf-8 -*-
"""任务6 实验可测性约束引擎
模式:
  --regression  已知集回归测试(737 已知活性分子过规则, 校准规则严宽)
  --dry         全库 dry-run(统计+Top50命中)
  --label       正式打标(compounds.csv 副本 + assayability 列)
"""
import os, sys, json, csv
sys.path.insert(0, r"D:\MMModel\Pocket2Mol")
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

CFG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "constraints_config.json")
CFG = json.load(open(CFG_PATH, encoding="utf-8")) if os.path.exists(CFG_PATH) else {}

# ---------- 规则集 (初始版, 回归测试后校准) ----------
FLUOR_SMARTS = [
    ("coumarin/chromone", "c1ccc2c(c1)[O,S][CX3]=[OX1]"),
    ("xanthene", "c1ccc2c(c1)Oc1ccccc12"),
    ("polyaromatic>=4rings", None),  # 启发式: 芳香环数>=4
    ("boron_fluorophore", "[B]"),
    ("cyanine_like", "[N+]=C1C=CC(N)=CC1"),
    ("long_conjugated_chain", None),  # 启发式
]
REACTIVE_SMARTS = [
    ("aldehyde", "[CX3H2](=O)[#6]"),
    ("acyl_chloride", "[CX3](=O)(Cl)"),
    ("sulfonyl_chloride", "[SX4](=O)(=O)(Cl)"),
    ("anhydride", "[CX3](=O)[OX2][CX3](=O)"),
    ("isocyanate", "[NX2]=[CX2]=[OX1]"),
    ("epoxide", "[OX2]1~[CX4]~[CX4]~1"),
]
# 软警示: 毒理学提示而非绝对排除 (D6 决策: 已知集存在 enone 活性分子, PAINS 已单独过滤)
SOFT_SMARTS = [
    ("enone_michael", "[CX3]=[CX3][CX3](=[OX1])"),
]
FLUOR_P = [(n, Chem.MolFromSmarts(s)) for n, s in FLUOR_SMARTS if s]
REACT_P = [(n, Chem.MolFromSmarts(s)) for n, s in REACTIVE_SMARTS]
SOFT_P = [(n, Chem.MolFromSmarts(s)) for n, s in SOFT_SMARTS]

def esol_logS(mol):
    """Delaney (ESOL) 近似: LogS(mol/L)。失败返回 None。"""
    try:
        mw = Descriptors.MolWt(mol)
        clogp = Crippen.MolLogP(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        na = mol.GetNumAtoms()
        nap = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        ap = nap / max(1, na)
        return 0.261 - 0.742 * clogp - 0.00634 * mw + 0.0228 * rotb - 0.352 * ap
    except Exception:
        return None

def classify(mol):
    """返回 (color, reasons)。color: green/yellow/red"""
    reasons = []
    if mol is None:
        return "yellow", ["parse_fail"]
    # 荧光
    for name, p in FLUOR_P:
        if p is not None and mol.HasSubstructMatch(p):
            reasons.append("fluor:" + name)
    # 荧光多并苯: 全碳稠合芳香环系 (蒽/芘/花级)。杂芳环与联苯组合荧光弱, 不计。
    ri = mol.GetRingInfo()
    car_rings = []
    for ring in ri.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if len(ring) >= 5 and all(a.GetIsAromatic() and a.GetSymbol() == "C" for a in atoms):
            car_rings.append(set(ring))
    if car_rings:
        # 环邻接图 (共享>=1键 = 稠合), 找最大连通稠合系
        parent = list(range(len(car_rings)))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for i in range(len(car_rings)):
            for j in range(i + 1, len(car_rings)):
                shared = sum(1 for b in mol.GetBonds()
                             if b.GetBeginAtomIdx() in car_rings[i] and b.GetEndAtomIdx() in car_rings[i]
                             and b.GetBeginAtomIdx() in car_rings[j] and b.GetEndAtomIdx() in car_rings[j])
                if shared >= 1:
                    pi, pj = find(i), find(j)
                    if pi != pj: parent[pi] = pj
        from collections import Counter
        comp = Counter(find(i) for i in range(len(car_rings)))
        fused = max(comp.values()) if comp else 0
        if fused >= 4:
            reasons.append("fluor:polyaromatic4")
        elif fused == 3:
            reasons.append("fluor:polyaromatic3")
    if any(a.GetSymbol() == "B" for a in mol.GetAtoms()) and Descriptors.MolWt(mol) > 200:
        reasons.append("fluor:boron")
    # 反应性
    for name, p in REACT_P:
        if mol.HasSubstructMatch(p):
            reasons.append("reactive:" + name)
    if reasons:
        return "red", reasons
    for name, p in SOFT_P:
        if mol.HasSubstructMatch(p):
            return "yellow", ["soft:" + name]
    # 溶解度
    ls = esol_logS(mol)
    if ls is None:
        return "yellow", ["solubility_calc_fail"]
    g = CFG.get("assayability", {})
    green = g.get("logS_green", -4.5)
    yellow = g.get("logS_yellow", -5.0)
    if ls >= green:
        return "green", ["logS=%.2f" % ls]
    elif ls >= yellow:
        return "yellow", ["logS=%.2f" % ls]
    else:
        return "red", ["logS=%.2f" % ls]

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "--dry"
    # ---------- 回归测试 ----------
    if mode == "--regression":
        known = [s for s in open(r"D:\MMModel\已知药物库\A2A_known_set.txt", encoding="utf-8").read().splitlines() if s.strip()]
        cnt = {"green": 0, "yellow": 0, "red": 0}
        red_hits = {}
        for s in known:
            m = Chem.MolFromSmiles(s)
            c, rs = classify(m)
            cnt[c] += 1
            if c == "red":
                for r in rs:
                    key = r.split("=")[0]
                    red_hits[key] = red_hits.get(key, 0) + 1
        n = len(known)
        print("已知集回归测试: %d 分子" % n)
        print("  green %.1f%% | yellow %.1f%% | red %.1f%%" % (100*cnt["green"]/n, 100*cnt["yellow"]/n, 100*cnt["red"]/n))
        print("  红标原因分布:", red_hits)
        rate = cnt["red"] / n
        if rate > 0.30:
            print("  判定: 规则过严 (红标>30%), 需校准")
        elif rate > 0.15:
            print("  判定: 偏严, 关注 top 命中规则")
        else:
            print("  判定: 规则校准通过")
        return

    # ---------- dry-run / label ----------
    lib = r"D:\MMModel\化合物库"
    rows = list(csv.DictReader(open(os.path.join(lib, "compounds.csv"), encoding="utf-8")))
    cnt = {"green": 0, "yellow": 0, "red": 0}
    red_hits = {}
    out = []
    for r in rows:
        m = Chem.MolFromSmiles(r["smiles"])
        c, rs = classify(m)
        cnt[c] += 1
        if c == "red":
            for x in rs:
                red_hits[x.split("=")[0]] = red_hits.get(x.split("=")[0], 0) + 1
        out.append((c, rs, r))
    print("全库 dry-run: %d 分子" % len(rows))
    print("  green %.1f%% | yellow %.1f%% | red %.1f%%" % (100*cnt["green"]/len(rows), 100*cnt["yellow"]/len(rows), 100*cnt["red"]/len(rows)))
    print("  红标原因分布:", red_hits)

    top = list(csv.DictReader(open(os.path.join(lib, "top_candidates.csv"), encoding="utf-8-sig")))
    tcnt = {"green": 0, "yellow": 0, "red": 0}
    tred = []
    smap = {r["smiles"]: (c, rs) for c, rs, r in out}
    for r in top:
        c, rs = smap.get(r["smiles"], ("yellow", ["not_scored"]))
        tcnt[c] += 1
        if c == "red":
            tred.append((r["rank"], rs))
    print("Top50: green %d | yellow %d | red %d" % (tcnt["green"], tcnt["yellow"], tcnt["red"]))
    if tred:
        print("  Top50 红标明细:")
        for rk, rs in tred[:10]:
            print("    rank %s: %s" % (rk, rs))

    if mode == "--label":
        dst = os.path.join(lib, "compounds_tagged.csv")
        with open(dst, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            head = list(rows[0].keys()) + ["assayability", "assay_reasons"]
            w.writerow(head)
            for (c, rs, r) in out:
                w.writerow([r[k] for k in rows[0].keys()] + [c, ";".join(rs)])
        print("已写", dst)

if __name__ == "__main__":
    main()
