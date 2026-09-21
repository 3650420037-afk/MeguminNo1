# -*- coding: utf-8 -*-
"""最终送验候选整合推荐 (五维 + 姿态置信度)

维度来源 (全部为已产出的实验数据, 零新增计算):
  dock    : Vina 打分 (docking_a2a_library)
  qed     : 类药性 (compounds_tagged)
  si      : A2A/A1 选择性百分位差 (selectivity_full)
  assay   : 实验可测性三色 (compounds_tagged)
  safety  : ADMET 警示数 (admet_screen)
  pose    : 结合模式置信度 = 最优-次优能量差 (pose_consistency)

分档规则 (透明可审计):
  第一梯队 (Tier 1): vina <= -11 且 QED >= 0.6 且 SI > 0 且 可测非 red 且 alerts <= 2 且 姿态能量差 >= 0.5
  第二梯队 (Tier 2): vina <= -10.5 且 可测非 red 且 alerts <= 3 (其余维度不作硬门槛)
  输出按 (姿态置信度, -vina) 排序 —— 优先送验"打分强且结合模式可信"的分子

用法: python scripts/final_recommendation.py
"""
import os, csv
LIB = r"D:\MMModel\化合物库"
DOCK = r"D:\MMModel\Pocket2Mol\outputs\docking_a2a_library\summary.csv"
POSEDIR = r"D:\MMModel\Pocket2Mol\outputs\docking_a2a_library"


def main():
    tag = {r["smiles"]: r for r in csv.DictReader(open(os.path.join(LIB, "compounds_tagged.csv"), encoding="utf-8"))}
    dock = {r["smiles"]: float(r["vina_score"]) for r in csv.DictReader(open(DOCK, encoding="utf-8"))
            if r["status"] == "ok" and r["vina_score"]}
    name2smi = {r["name"]: r["smiles"] for r in csv.DictReader(open(DOCK, encoding="utf-8"))}
    sel = {r["smiles"]: float(r["SI百分位点"]) for r in csv.DictReader(open(os.path.join(LIB, "selectivity_full.csv"), encoding="utf-8-sig"))}
    adm = {r["smiles"]: int(r["alerts"]) for r in csv.DictReader(open(os.path.join(LIB, "admet_screen.csv"), encoding="utf-8-sig"))}
    top = {r["smiles"]: int(r["rank"]) for r in csv.DictReader(open(os.path.join(LIB, "top_candidates.csv"), encoding="utf-8-sig"))}
    pose = {}
    pp = os.path.join(LIB, "pose_consistency.csv")
    if os.path.exists(pp):
        for r in csv.DictReader(open(pp, encoding="utf-8-sig")):
            gap = r.get("energy_gap", "")
            pose[r["smiles"]] = (float(gap) if gap not in ("", None) else None)

    rows = []
    for smi, v in dock.items():
        t = tag.get(smi)
        if t is None:
            continue
        si = sel.get(smi)
        al = adm.get(smi)
        pg = pose.get(smi)
        rows.append(dict(smiles=smi, vina=v, qed=float(t["qed"]), sa=float(t["sa"]), mw=float(t["mw"]),
                         si=(si if si is not None else None), assay=t["assayability"],
                         alerts=(al if al is not None else None), pose_gap=pg,
                         extrapolation=t.get("extrapolation", "?"), top_rank=top.get(smi, 0)))

    def tier(r):
        if (r["vina"] <= -11 and r["qed"] >= 0.6
                and (r["si"] is not None and r["si"] > 0)      # 不能用 `or -999`: SI 恰为 0.0 会被当假值
                and r["assay"] != "red" and (r["alerts"] is None or r["alerts"] <= 2)
                and (r["pose_gap"] is None or r["pose_gap"] >= 0.5)):
            return 1
        if r["vina"] <= -10.5 and r["assay"] != "red" and (r["alerts"] is None or r["alerts"] <= 3):
            return 2
        return 0

    for r in rows:
        r["tier"] = tier(r)
    t1 = sorted([r for r in rows if r["tier"] == 1], key=lambda r: (-(r["pose_gap"] or 0), r["vina"]))
    t2 = sorted([r for r in rows if r["tier"] == 2], key=lambda r: (-(r["pose_gap"] or 0), r["vina"]))

    print("参与整合: %d 个分子" % len(rows))
    print("第一梯队 (全维度达标 + 姿态可信): %d" % len(t1))
    print("第二梯队 (打分强 + 可测 + 低 ADMET 风险): %d" % len(t2))
    dst = os.path.join(LIB, "final_recommendation.csv")
    with open(dst, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["tier", "rank", "vina_kcal", "qed", "sa", "mw", "SI百分位", "可测性",
                    "admet_alerts", "姿态能量差", "自洽性(域外标签)", "原Top50排名", "smiles"])
        for ti, group in ((1, t1), (2, t2)):
            for i, r in enumerate(group[:100], 1):
                w.writerow([ti, i, round(r["vina"], 2), round(r["qed"], 3), round(r["sa"], 2),
                            round(r["mw"], 1),
                            (round(r["si"], 1) if r["si"] is not None else "-"),
                            r["assay"], (r["alerts"] if r["alerts"] is not None else "-"),
                            (round(r["pose_gap"], 2) if r["pose_gap"] is not None else "-"),
                            r["extrapolation"], r["top_rank"] or "-", r["smiles"]])
    print("已写:", dst)
    print()
    print("第一梯队前 10 (优先送验):")
    print("  %-4s %-7s %-6s %-7s %-6s %-5s %-6s %s" % ("#", "vina", "QED", "SI", "可测", "ADMET", "姿态差", "smiles"))
    for i, r in enumerate(t1[:10], 1):
        print("  %-4d %-7.2f %-6.3f %-+7.1f %-6s %-5s %-6s %s" % (
            i, r["vina"], r["qed"], r["si"] or 0, r["assay"],
            r["alerts"] if r["alerts"] is not None else "-",
            ("%.2f" % r["pose_gap"]) if r["pose_gap"] is not None else "-", r["smiles"][:40]))


if __name__ == "__main__":
    main()
