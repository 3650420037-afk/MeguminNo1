# -*- coding: utf-8 -*-
"""C11 ADMET 规则化筛查 (透明可审计, 无外部模型依赖)

在送验前对候选做早期排雷: 心脏毒性(hERG)/CYP 抑制/致突变(AMES)/BBB 穿透/
溶解度(ESOL)/血浆蛋白结合倾向, 输出每分子警示标签与可开发性分级。

设计原则:
  - 全部为可解释的结构警示规则 + 公开理化阈值, 不使用黑盒 QSAR 模型
    (夜间无人值守场景下零下载依赖、零外部服务)
  - 明确标注为"规则化筛查", 不是定量预测; 阳性结果只表示"需关注/需实验确认"
  - 规则命中不计入分子筛除, 只打标签与分级 (与项目"标签优先, 不物理剔除"原则一致)

用法: python scripts/admet_screen.py [--library DIR] [--top50]
"""
import os, sys, csv, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import RESULTS, src_on_path
src_on_path()
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
RDLogger.DisableLog("rdApp.*")

# ---------------- 结构警示规则 (SMARTS) ----------------
HERG_SMARTS = [
    ("tertiary_amine", "[NX3;H0;!$(N=*);!$(N#*);!$([N+]);!$(Na)]"),
    ("basic_secondary_amine", "[NX3;H1;!$(N=*);!$(Na);!$(NC=O)]"),
]
CYP_SMARTS = [
    ("furan", "c1ccoc1"),
    ("nitro_aromatic", "[$([NX3](=O)=O),$([NX3+](=O)[O-])][c]"),
    ("pyridine_like", "c1ccncc1"),
    ("hydrazine", "[NX3][NX3]"),
]
AMES_SMARTS = [
    ("aromatic_nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])][c]"),
    ("aromatic_amine", "[NX3;H2,H1;!$(NC=O)][c]"),
    ("epoxide", "[OX2r3]1[CX4r3][CX4r3]1"),
    ("nitrogen_mustard", "[NX3]([CX4][Cl,Br,I])[CX4][Cl,Br,I]"),
    ("polyhalogen", "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]"),
    ("azo", "[NX2]=[NX2]"),
]
_P = lambda lst: [(n, Chem.MolFromSmarts(s)) for n, s in lst]

def screen(mol):
    if mol is None:
        return None
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    na = mol.GetNumAtoms()
    ap = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic()) / max(1, na)
    logS = 0.261 - 0.742*logp - 0.00634*mw + 0.0228*rotb - 0.352*ap

    herg_hits = [n for n, p in _P(HERG_SMARTS) if p and mol.HasSubstructMatch(p)]
    # hERG 风险还要看两亲性: 碱性氮 + 芳环>=2 + LogP>3
    herg_risk = "high" if (herg_hits and arom >= 2 and logp > 3) else \
                ("medium" if herg_hits and logp > 2 else "low")
    cyp_hits = [n for n, p in _P(CYP_SMARTS) if p and mol.HasSubstructMatch(p)]
    cyp_risk = "high" if (len(cyp_hits) >= 2 and mw > 400) else ("medium" if cyp_hits else "low")
    ames_hits = [n for n, p in _P(AMES_SMARTS) if p and mol.HasSubstructMatch(p)]
    ames_risk = "high" if ("nitrogen_mustard" in ames_hits or "epoxide" in ames_hits or "azo" in ames_hits) \
                else ("medium" if ames_hits else "low")
    # BBB (CNS MPO 简化): TPSA<90 且 MW<450 且 1<=LogP<=4 且 HBD<=3 且 无酸性基团
    acidic = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H1]")) or \
             mol.HasSubstructMatch(Chem.MolFromSmarts("[SX4](=O)(=O)[OX2H1]"))
    bbb_permeable = (tpsa < 90 and mw < 450 and 1.0 <= logp <= 4.0 and hbd <= 3 and not acidic)
    pbb_risk = "high" if (logp > 4.5 or (logp > 3.5 and acidic)) else ("medium" if logp > 3.5 else "low")
    sol = "poor" if logS < -6 else ("moderate" if logS < -4.5 else "good")

    alerts = 0
    alerts += {"high": 2, "medium": 1, "low": 0}[herg_risk]
    alerts += {"high": 2, "medium": 1, "low": 0}[cyp_risk]
    alerts += {"high": 2, "medium": 1, "low": 0}[ames_risk]
    alerts += 1 if pbb_risk == "high" else 0
    alerts += 1 if sol == "poor" else 0
    grade = "A" if alerts == 0 else ("B" if alerts <= 2 else "C")
    return dict(mw=round(mw,1), logp=round(logp,2), tpsa=round(tpsa,1), logS=round(logS,2),
                herg=herg_risk, herg_hits=";".join(herg_hits),
                cyp=cyp_risk, cyp_hits=";".join(cyp_hits),
                ames=ames_risk, ames_hits=";".join(ames_hits),
                bbb="permeable" if bbb_permeable else "not_permeable",
                pbb=pbb_risk, solubility=sol, alerts=alerts, grade=grade)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", default=RESULTS)
    ap.add_argument("--top50", action="store_true", help="只处理 top_candidates.csv")
    args = ap.parse_args()

    LIB = args.library
    src = os.path.join(LIB, "top_candidates.csv" if args.top50 else "compounds.csv")
    rows = list(csv.DictReader(open(src, encoding="utf-8-sig" if args.top50 else "utf-8")))
    out_rows = []
    from collections import Counter
    grades = Counter(); risks = Counter()
    for r in rows:
        mol = Chem.MolFromSmiles(r["smiles"])
        s = screen(mol)
        if s is None: continue
        rec = {"smiles": r["smiles"]}
        rec.update({k: v for k, v in s.items()})
        out_rows.append(rec)
        grades[s["grade"]] += 1
        risks["herg_" + s["herg"]] += 1
        risks["cyp_" + s["cyp"]] += 1
        risks["ames_" + s["ames"]] += 1

    dst = os.path.join(LIB, "admet_screen_top50.csv" if args.top50 else "admet_screen.csv")
    with open(dst, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    n = len(out_rows)
    print("ADMET 规则化筛查: %d 个分子 -> %s" % (n, dst))
    print("可开发性分级: " + " | ".join("%s=%d(%.0f%%)" % (g, grades[g], 100.0*grades[g]/n) for g in ("A","B","C") if g in grades))
    print("风险分布: " + " | ".join("%s=%d" % (k, v) for k, v in sorted(risks.items())))
    if args.top50:
        print()
        print("Top50 明细:")
        for r in out_rows[:10]:
            print("  %s %s | hERG=%s CYP=%s AMES=%s BBB=%s sol=%s alerts=%d" % (
                r["grade"], r["smiles"][:34], r["herg"], r["cyp"], r["ames"], r["bbb"], r["solubility"], r["alerts"]))


if __name__ == "__main__":
    main()
