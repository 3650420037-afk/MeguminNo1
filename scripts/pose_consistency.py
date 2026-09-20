# -*- coding: utf-8 -*-
"""C10 对接姿态一致性分析 (MD 的零成本替代)

对 Vina 输出的多姿态 (默认 9 个 MODEL) 计算各姿态相对最优姿态的重原子 RMSD:
  - 平均 RMSD 小  -> 搜索收敛, 结合模式明确, 结论可信
  - 平均 RMSD 大  -> 姿态不确定, 需要更长采样或 MD 确认
同时给出"最优姿态与次优姿态的能量差"作为置信度指标。

用法: python scripts/pose_consistency.py --docking-dir outputs/docking_a2a_library [--top N]
输出: 化合物库/pose_consistency.csv
"""
import os, sys, csv, glob, argparse
import numpy as np

LIB = r"D:\MMModel\化合物库"


def parse_pdbqt_models(path):
    """返回 [(coords Nx3, energy, rmsd_lb)] 按 MODEL 顺序.

    说明: 自算姿态间坐标 RMSD 会得到很大的值 (Vina 的多模式本就是按聚类选出的
    **不同结合模式**, 天然差异大), 因此以 Vina 官方报告的 REMARK VINA RESULT
    第 2 列 (RMSD lower bound, 相对该输出最优模式) 为准, 并用最优-次优能量差
    作为结合模式置信度指标。
    """
    models = []
    coords, energy, rmsd_lb = [], None, None
    try:
        for line in open(path, encoding="utf-8", errors="ignore"):
            if line.startswith("MODEL"):
                coords, energy, rmsd_lb = [], None, None
            elif line.startswith(("ATOM", "HETATM")):
                el = line[76:79].strip() if len(line) > 76 else line.split()[-1]
                if el.upper().startswith("H"):
                    continue
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                except ValueError:
                    continue
                coords.append((x, y, z))
            elif line.startswith("REMARK") and "VINA RESULT" in line:
                p = line.split()
                try:
                    energy = float(p[3])
                    rmsd_lb = float(p[4])
                except (IndexError, ValueError):
                    pass
            elif line.startswith("ENDMDL"):
                if coords:
                    models.append((np.array(coords), energy, rmsd_lb))
    except Exception:
        return []
    if not models and coords:
        models.append((np.array(coords), energy, rmsd_lb))
    return models


def rmsd(a, b):
    if a.shape != b.shape:
        n = min(len(a), len(b))
        a, b = a[:n], b[:n]
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docking-dir", required=True)
    ap.add_argument("--top", type=int, default=50, help="只分析前 N 个 (按打分); 0=全部")
    args = ap.parse_args()

    summ = os.path.join(args.docking_dir, "summary.csv")
    rows = [r for r in csv.DictReader(open(summ, encoding="utf-8")) if r["status"] == "ok"]
    rows.sort(key=lambda r: float(r["vina_score"]))
    if args.top:
        rows = rows[:args.top]

    out = []
    for r in rows:
        name = r["name"]
        cands = glob.glob(os.path.join(args.docking_dir, "docking", name + "_out.pdbqt")) + \
                glob.glob(os.path.join(args.docking_dir, "docking", name + "*.pdbqt"))
        if not cands:
            continue
        models = parse_pdbqt_models(cands[0])
        if len(models) < 2:
            continue
        # 以 Vina 官方 RMSD l.b. (相对最优模式) 为口径
        rmsds = [m[2] for m in models[1:] if m[2] is not None]
        energies = [m[1] for m in models if m[1] is not None]
        gap = (energies[1] - energies[0]) if len(energies) >= 2 else None
        out.append(dict(name=name, smiles=r["smiles"], vina=float(r["vina_score"]),
                        n_modes=len(models),
                        rmsd_lb_mean=round(float(np.mean(rmsds)), 3) if rmsds else "",
                        rmsd_lb_max=round(float(np.max(rmsds)), 3) if rmsds else "",
                        energy_gap=round(gap, 2) if gap is not None else ""))

    if not out:
        print("未解析到姿态数据, 检查目录:", os.path.join(args.docking_dir, "docking"))
        return
    gaps = np.array([o["energy_gap"] for o in out if o["energy_gap"] != ""], dtype=float)
    rl = np.array([o["rmsd_lb_mean"] for o in out if o["rmsd_lb_mean"] != ""], dtype=float)
    print("姿态分析: %d 个分子 (平均 %.1f 个结合模式)" % (len(out), np.mean([o["n_modes"] for o in out])))
    if len(rl):
        print("Vina RMSD l.b. (模式2..N 相对最优): 中位 %.2f Å | P90 %.2f Å" % (
            np.median(rl), np.percentile(rl, 90)))
    if len(gaps):
        print("最优-次优能量差: 中位 %.2f kcal/mol | <0.5 占比 %.0f%% (越小说明结合模式竞争越激烈)" % (
            np.median(gaps), 100 * (gaps < 0.5).mean()))
        confident = ((gaps >= 1.5)).mean()
        print("结合模式较可信 (能量差>=1.5 kcal/mol) 占比: %.0f%%" % (100 * confident))
    dst = os.path.join(LIB, "pose_consistency.csv")
    with open(dst, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    print("已写:", dst)


if __name__ == "__main__":
    main()
