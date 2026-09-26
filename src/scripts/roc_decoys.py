# -*- coding: utf-8 -*-
"""C9 反筛/诱饵验证: 打分区分度 ROC/AUC (DUD-E 风格)

正类 = 已知 A2A 活性分子 (docking_known_reverse, 737 个已知集的对接结果)
负类 = 属性匹配诱饵 (docking_decoys, MW/LogP 分箱匹配已知集, 1:2 比例)

判读:
  AUC ~0.5  打分无区分度 -> 全部"强结合"结论需大打折扣
  AUC 0.6-0.7 弱区分度
  AUC >=0.7  可用区分度 (DUD-E 上 Vina 的典型水平为 0.6-0.7)
  EF1%/EF5% 富集因子: 打分前 1%/5% 中真实活性分子的富集倍数

用法: python scripts/roc_decoys.py
输出: results/roc_decoys.csv (阈值扫描) + 控制台汇总
"""
import os, csv, sys, tempfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import OUTPUTS, RESULTS, src_on_path
src_on_path()

LIB = RESULTS
POS = os.path.join(OUTPUTS, "docking_known_reverse", "summary.csv")
NEG = os.path.join(OUTPUTS, "docking_decoys", "summary.csv")


def load_scores(path):
    out = []
    if not os.path.exists(path):
        return out
    for r in csv.DictReader(open(path, encoding="utf-8")):
        if r.get("status") == "ok" and r.get("vina_score"):
            try:
                out.append(float(r["vina_score"]))
            except ValueError:
                pass
    return out


def load_scores_from_log(logpath):
    """对接仍在进行时的 fallback: 从实时日志解析已完成的 ok 分值"""
    import re, os as _os
    out = []
    if not _os.path.exists(logpath):
        return out
    pat = re.compile(r"score=(-?\d+\.\d+)\s+\(ok\)")
    for line in open(logpath, encoding="utf-8", errors="ignore"):
        m = pat.search(line)
        if m:
            out.append(float(m.group(1)))
    return out


def roc(pos, neg, thresholds):
    P, N = len(pos), len(neg)
    rows = []
    for t in thresholds:
        tp = sum(1 for s in pos if s <= t)
        fp = sum(1 for s in neg if s <= t)
        rows.append((t, tp / P, fp / N, tp, fp))
    return rows


def auc_from_roc(rows):
    """梯形法 (TPR 对 FPR 积分), rows 按阈值从松到紧"""
    trap = getattr(np, "trapezoid", None) or np.trapz   # 兼容 numpy<2 (trapz 已弃用)
    pts = sorted(set((f, t) for _, t, f, _, _ in rows))
    xs = [0.0] + [p[0] for p in pts]
    ys = [0.0] + [p[1] for p in pts]
    return float(trap(ys, xs))


def main():
    pos = load_scores(POS)
    neg = load_scores(NEG)
    if not neg:
        neg = load_scores_from_log(os.path.join(tempfile.gettempdir(), "dock_decoys.log"))
        if neg:
            print("(诱饵对接仍在进行: 使用实时日志中已完成的 %d 个分数)" % len(neg))
    print("正类(已知活性): %d | 负类(诱饵): %d" % (len(pos), len(neg)))
    if len(pos) < 50 or len(neg) < 50:
        print("样本不足, 等待对接完成")
        return
    base = len(pos) / (len(pos) + len(neg))
    print("基准阳性率: %.3f" % base)
    thr = np.arange(-14.0, -4.0, 0.1)
    rows = roc(pos, neg, thr)
    auc = auc_from_roc(rows)
    # 最佳阈值: Youden index = TPR - FPR
    best = max(rows, key=lambda r: r[1] - r[2])
    print()
    print("=== 核心指标 ===")
    print("AUC = %.3f" % auc)
    print("最佳阈值 (Youden): vina <= %.1f -> 灵敏度(TPR) %.1f%% | 假阳性率(FPR) %.1f%% | 富集倍数 %.2fx" % (
        best[0], 100 * best[1], 100 * best[2], best[1] / max(1e-9, best[2])))
    # 富集因子: 分数前 1% / 5% / 10% 中真实活性占比
    allsc = [(s, 1) for s in pos] + [(s, 0) for s in neg]
    allsc.sort()
    for frac in (0.01, 0.05, 0.10):
        k = max(1, int(len(allsc) * frac))
        hit = sum(lbl for _, lbl in allsc[:k]) / k
        print("前 %2d%% 打分中真实活性占比 %.1f%% -> 富集因子 %.2fx" % (100 * frac, 100 * hit, hit / base))
    # 各阈值区间对照
    print()
    print("阈值扫描 (vina <= t):")
    print("  %-8s %-9s %-9s %-9s" % ("t", "TPR", "FPR", "富集"))
    for t in (-12, -11, -10, -9, -8):
        r = next((x for x in rows if abs(x[0] - t) < 0.05), None)
        if r:
            print("  %-8.1f %-9.1f %-9.1f %-9.2f" % (r[0], 100 * r[1], 100 * r[2], r[1] / max(1e-9, r[2])))
    dst = os.path.join(LIB, "roc_decoys.csv")
    with open(dst, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["threshold", "TPR", "FPR", "TP", "FP"])
        for r in rows:
            w.writerow([round(r[0], 2), round(r[1], 4), round(r[2], 4), r[3], r[4]])
    print()
    print("已写:", dst)
    print()
    if auc >= 0.7:
        verdict = "打分有可用区分度 (DUD-E 上 Vina 典型水平 0.6-0.7), 打分结论可信"
    elif auc >= 0.6:
        verdict = "弱区分度 (与文献 Vina 水平一致), 打分可用于粗筛, 精细排序需实验/更严方法"
    else:
        verdict = "区分度不足, 全部'强结合'结论须谨慎表述"
    print("判读: AUC=%.3f -> %s" % (auc, verdict))


if __name__ == "__main__":
    main()
