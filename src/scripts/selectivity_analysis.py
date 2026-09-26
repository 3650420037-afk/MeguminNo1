# -*- coding: utf-8 -*-
"""全库亚型选择性谱重算 (A2A vs A1, 百分位归一)

此前由内联脚本执行且受 summary.csv BOM 丢样影响(2930 而非 2931 配对)。本脚本固化逻辑:
  SI = pct(A2A 打分, A2A 池) - pct(A1 打分, A1 池)     # 百分位归一双受体系统偏差
输出: results/selectivity_full.csv
"""
import os, csv, bisect, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import OUTPUTS, RESULTS, src_on_path
src_on_path()

LIB = RESULTS
OUTROOT = OUTPUTS


def load(f):
    if not os.path.exists(f):
        return {}
    return {r["smiles"]: float(r["vina_score"]) for r in csv.DictReader(open(f, encoding="utf-8-sig"))
            if r["status"] == "ok" and r["vina_score"]}


def main():
    a2a = load(os.path.join(OUTROOT, "docking_a2a_library", "summary.csv"))
    a1 = {}
    # 仅取对 A1 受体(5UEN)打分的两份结果; docking_known_reverse 是"已知 A2A 分子对 A2A
    # 受体"的阳性对照, 不能混入 A1 池
    for d in ("docking_a1_rest", "docking_a1_top200"):
        a1.update(load(os.path.join(OUTROOT, d, "summary.csv")))
    common = sorted(set(a2a) & set(a1))
    print("A2A 打分 %d | A1 打分 %d | 配对 %d" % (len(a2a), len(a1), len(common)))
    pool_a2a = [a2a[s] for s in common]
    pool_a1 = [a1[s] for s in common]

    def pct(sc, pool):
        s = sorted(pool)
        i = bisect.bisect_left(s, sc)
        return 100.0 * (len(s) - i) / len(s)

    tag = {}
    tp = os.path.join(LIB, "compounds_tagged.csv")
    if os.path.exists(tp):
        tag = {r["smiles"]: r for r in csv.DictReader(open(tp, encoding="utf-8-sig"))}
    top = {}
    tcp = os.path.join(LIB, "top_candidates.csv")
    if os.path.exists(tcp):
        top = {r["smiles"]: int(r["rank"]) for r in csv.DictReader(open(tcp, encoding="utf-8-sig"))}

    out = []
    for smi in common:
        si = pct(a2a[smi], pool_a2a) - pct(a1[smi], pool_a1)
        out.append((si, smi, a2a[smi], a1[smi],
                    tag.get(smi, {}).get("assayability", "?"), top.get(smi, 0)))
    out.sort(reverse=True)
    a = np.array([x[0] for x in out])
    t = np.array([x[0] for x in out if x[5] > 0])
    print("全库 SI: 均值 %+.1f 中位 %+.1f | SI>0 %.0f%% | SI>25: %d" % (
        a.mean(), np.median(a), 100 * (a > 0).mean(), (a > 25).sum()))
    if len(t):
        print("Top50 SI: 均值 %+.1f 中位 %+.1f | SI>0 %d/%d" % (t.mean(), np.median(t), (t > 0).sum(), len(t)))
    with open(os.path.join(LIB, "selectivity_full.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["SI百分位点", "A2A分", "A1分", "可测性", "Top50排名", "smiles"])
        for si, smi, va, v1, asy, trk in out:
            w.writerow([round(si, 1), va, v1, asy, trk, smi])
    print("已写: selectivity_full.csv (%d 行)" % len(out))


if __name__ == "__main__":
    main()
