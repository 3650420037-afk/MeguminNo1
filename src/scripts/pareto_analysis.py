# -*- coding: utf-8 -*-
"""A3 离线多目标帕累托分析 (零 GPU 成本)

对已完成的多维标注做帕累托支配排序, 找出"在类药性/结合强度/选择性/可测性/
ADMET 风险上没有短板"的非支配前沿集合。

维度 (5 目标, 全部转为越大越好):
  dock   : -vina_score        (结合强度)
  qed    : QED                (类药性)
  si     : 选择性百分位差      (A2A vs A1, 来自 selectivity_full.csv)
  assay  : 可测性系数          (green 1.0 / yellow 0.9 / red 0.7)
  safety : -admet_alerts      (ADMET 警示越少越好)

用法: python scripts/pareto_analysis.py
输出: results/pareto_front.csv + 控制台汇总
"""
import os, csv, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import OUTPUTS, RESULTS, src_on_path
src_on_path()

LIB = RESULTS


def load():
    tag = {r["smiles"]: r for r in csv.DictReader(open(os.path.join(LIB, "compounds_tagged.csv"), encoding="utf-8"))}
    dock = {r["smiles"]: float(r["vina_score"]) for r in csv.DictReader(open(
        os.path.join(OUTPUTS, "docking_a2a_library", "summary.csv"), encoding="utf-8"))
        if r["status"] == "ok" and r["vina_score"]}
    sel = {r["smiles"]: float(r["SI百分位点"]) for r in csv.DictReader(
        open(os.path.join(LIB, "selectivity_full.csv"), encoding="utf-8-sig"))}
    adm = {r["smiles"]: int(r["alerts"]) for r in csv.DictReader(
        open(os.path.join(LIB, "admet_screen.csv"), encoding="utf-8-sig"))}
    top = {r["smiles"]: int(r["rank"]) for r in csv.DictReader(
        open(os.path.join(LIB, "top_candidates.csv"), encoding="utf-8-sig"))}
    acoef = {"green": 1.0, "yellow": 0.9, "red": 0.7}
    rows = []
    for smi, t in tag.items():
        if smi not in dock or smi not in sel or smi not in adm:
            continue
        rows.append(dict(
            smiles=smi,
            dock=-dock[smi],
            qed=float(t["qed"]),
            si=sel[smi],
            assay=acoef.get(t["assayability"], 0.9),
            safety=-float(adm[smi]),
            top_rank=top.get(smi, 0),
            assay_label=t["assayability"],
            extrapolation=t.get("extrapolation", "?"),
        ))
    return rows


def pareto_front(rows, keys):
    """返回非支配解下标 (向量化两两比较)"""
    M = np.array([[r[k] for k in keys] for r in rows], dtype=float)
    n = len(M)
    dominated = np.zeros(n, dtype=bool)
    # 分块比较避免内存爆炸
    for i in range(n):
        if dominated[i]:
            continue
        ge = (M >= M[i]).all(axis=1)      # 各维不劣于 i
        gt = (M > M[i]).any(axis=1)       # 至少一维优于 i
        dominators = ge & gt
        if dominators.any():
            dominated[i] = True
    return np.where(~dominated)[0]


def main():
    rows = load()
    keys = ["dock", "qed", "si", "assay", "safety"]
    print("参与分析分子: %d (同时具备 docking/选择性/ADMET/可测性标签)" % len(rows))
    idx = pareto_front(rows, keys)
    front = [rows[i] for i in idx]
    front.sort(key=lambda r: (-r["dock"], -r["qed"]))
    print("帕累托前沿规模: %d (%.1f%%)" % (len(front), 100.0 * len(front) / len(rows)))
    # 对照: 仅按综合分(0.6 dock_norm + 0.4 qed)排序的 Top50 有多少落在前沿
    fr = set(r["smiles"] for r in front)
    top50 = sorted([r for r in rows if r["top_rank"]], key=lambda r: r["top_rank"])[:50]
    hit = sum(1 for r in top50 if r["smiles"] in fr)
    print("原 Top50 中位于帕累托前沿的: %d/50" % hit)
    print()
    print("前沿成员特征 (中位): dock %.2f | QED %.3f | SI %+.1f | 可测 %s | alerts %.0f" % (
        np.median([r["dock"] for r in front]), np.median([r["qed"] for r in front]),
        np.median([r["si"] for r in front]),
        ", ".join(sorted(set(r["assay_label"] for r in front))),
        np.median([-r["safety"] for r in front])))
    dst = os.path.join(LIB, "pareto_front.csv")
    with open(dst, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["rank", "vina_kcal", "qed", "SI百分位", "可测性", "admet_alerts",
                    "自洽性", "原Top50排名", "smiles"])
        for i, r in enumerate(front[:200], 1):
            w.writerow([i, round(-r["dock"], 2), round(r["qed"], 3), round(r["si"], 1),
                        r["assay_label"], int(-r["safety"]), r["extrapolation"],
                        r["top_rank"] or "-", r["smiles"]])
    print("已写:", dst, "(前 200 名)")
    print()
    print("前沿前 10:")
    for i, r in enumerate(front[:10], 1):
        print("  #%-2d vina %.2f | QED %.3f | SI %+.1f | %s | alerts %d | Top%s" % (
            i, -r["dock"], r["qed"], r["si"], r["assay_label"], int(-r["safety"]), r["top_rank"] or "-"))


if __name__ == "__main__":
    main()
