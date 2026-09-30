# -*- coding: utf-8 -*-
"""把 v10 库的 ADMET 规则化筛查结果回填到候选表（results*.csv）。

- 分级 -> `关键预测指标_ADMET分级`
- 警示数 -> `ADMET警报`
- logS -> `LogS`（若为空）
只按 SMILES 精确匹配；未匹配的保持空白并在末尾报告。
"""
import csv, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = r"D:\MMModel"

JOBS = [
    (os.path.join(ROOT, "results", "admet_screen.csv"), os.path.join(ROOT, "results", "results.csv"), "A2A"),
    (os.path.join(ROOT, "results", "D3", "admet_screen.csv"), os.path.join(ROOT, "results", "D3", "results.csv"), "D3"),
    (os.path.join(ROOT, "results", "5HT2B", "admet_screen.csv"), os.path.join(ROOT, "results", "5HT2B", "results.csv"), "5HT2B"),
]

for admet_p, cand_p, tag in JOBS:
    admet = {}
    with open(admet_p, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            admet[r["smiles"].strip()] = r
    with open(cand_p, encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f)
        fields = list(rd.fieldnames)
        rows = list(rd)
    gcol = "关键预测指标_ADMET分级"
    acol = "ADMET警报"
    scol = "LogS"
    n_fill = 0
    for r in rows:
        a = admet.get((r.get("候选SMILES") or "").strip())
        if not a:
            continue
        r[gcol] = a.get("grade", "")
        r[acol] = a.get("alerts", "")
        if not (r.get(scol) or "").strip():
            r[scol] = a.get("logS", "")
        n_fill += 1
    with open(cand_p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    grades = {}
    for r in rows:
        grades[r[gcol]] = grades.get(r[gcol], 0) + 1
    print("%-6s 候选 %d 个，回填 %d 个；分级分布 %s -> %s" % (
        tag, len(rows), n_fill, dict(sorted(grades.items())), os.path.relpath(cand_p, ROOT)))
