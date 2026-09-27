# -*- coding: utf-8 -*-
"""在已有对接数据上训练**轻量奖励代理**, 用交叉验证检验它是否真的可用。

背景(上一轮的负面结果)
----------------------
手工几何接触特征(`fast_reward_proxy.py`)速度达标(1.4 ms/分子)但
Spearman(代理分, -Vina) 仅 **+0.116** —— 因为 `buried` 在 Vina 自己优化过的位姿上
恒为 1.0、其余接触量也饱和, 特征没有区分度。
诊断给出的下一步是: **用真实对接数据做监督**, 让模型自己学权重与非线性组合,
而不是继续手工调特征。

本脚本做什么
------------
1. 复用 `fast_reward_proxy` 的口袋/位姿解析与几何特征(已验证正确);
2. 再加一组**与 Vina 打分函数同族的**特征: 每原子接触数分布(不只均值)、
   接触原子的元素构成(疏水/极性/芳香近似)、配体大小与柔性的代理(可旋转键数、
   重原子数)、以及"埋藏深度"分位数 —— 饱和的是均值, 分布仍有信息;
3. 用 **K 折交叉验证** 拟合 ridge 回归并报 **折外 Spearman** ——
   只看折外指标, 避免用训练集自证。

纪律
----
- 折外 |Spearman| 达不到可用水平(经验上 <0.3)就**如实判失败**, 不调参硬凑;
- 该代理即便可用, **也不能替代 A/B 主指标**(Vina), 只能当筛选器;
- ⚠ 已知局限: 特征是在 **Vina 优化过的位姿**上算的, 而 ReST 要筛的是
  **模型自己生成的位姿** —— 分布不同, 迁移性需另行验证。

用法
----
    python src/scripts/learn_reward_proxy.py \
        --summary outputs/mass_dock/A2A/summary_merged.csv \
        --poses outputs/mass_dock/A2A/docking --target A2A --limit 3000 --folds 5
"""
import argparse
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from fast_reward_proxy import (load_receptor, load_pose, score, _spearman,  # noqa: E402
                              CONTACT, HBOND, CLASH, POLAR)


def features(lig_xyz, lig_e, rec_xyz, rec_e):
    """几何特征(含分布刻画, 不只均值) + 配体自身描述量。"""
    d = np.linalg.norm(lig_xyz[:, None, :] - rec_xyz[None, :, :], axis=2)
    dmin = d.min(axis=1)
    m = d.shape[0]
    n_contact = (d < CONTACT).sum(axis=1)            # 每个配体原子的接触数
    lp = np.isin(lig_e, list(POLAR))
    rp = np.isin(rec_e, list(POLAR))
    polar_pairs = (d < HBOND)[lp][:, rp].sum() if lp.any() else 0
    hydro_pairs = (d < CONTACT)[~lp][:, ~rp].sum()
    s = score(lig_xyz, lig_e, rec_xyz, rec_e)
    f = [
        float(m),                                        # 重原子数(尺寸)
        float(dmin.mean()), float(dmin.min()),           # 埋藏深度: 均值/最深
        float(np.percentile(dmin, 25)), float(np.percentile(dmin, 75)),
        float(n_contact.mean()), float(n_contact.max()), # 接触数: 均值/最大
        float(np.percentile(n_contact, 25)), float(np.percentile(n_contact, 75)),
        float((n_contact == 0).mean()),                  # 完全无接触的原子比例
        float(hydro_pairs), float(polar_pairs),
        float(polar_pairs) / max(1.0, m),                # 单位原子极性接触
        float(s["buried"]), float(s["contact"]), float(s["hydro"]), float(s["clash"]),
        float((d < CONTACT).sum()) / max(1.0, m * m),     # 接触密度
    ]
    return np.asarray(f, dtype=np.float64)


def ridge_fit(X, y, lam):
    Xb = np.hstack([X, np.ones((len(X), 1))])
    A = Xb.T @ Xb + lam * np.eye(Xb.shape[1])
    return np.linalg.solve(A, Xb.T @ y)


def main():
    ap = argparse.ArgumentParser(description="用对接数据训练轻量奖励代理并交叉验证")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--poses", required=True)
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--limit", type=int, default=3000)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--lam", type=float, default=1.0)
    args = ap.parse_args()

    rec_xyz, rec_e = load_receptor(args.target)
    rows = [r for r in csv.DictReader(open(args.summary, encoding="utf-8-sig"))
            if r.get("status") == "ok"][:args.limit]
    print("口袋原子 %d; 读入 %d 条 ok 记录" % (len(rec_xyz), len(rows)))

    X, y = [], []
    for r in rows:
        p = os.path.join(args.poses, r["name"] + "_out.pdbqt")
        if not os.path.exists(p):
            continue
        lx, le = load_pose(p)
        if lx is None:
            continue
        try:
            v = float(r["vina_score"])
        except ValueError:
            continue
        X.append(features(lx, le, rec_xyz, rec_e))
        y.append(-v)                    # 目标: -Vina(越大越好), 便于解释
    X = np.asarray(X)
    y = np.asarray(y)
    print("有效样本 %d, 特征维度 %d" % (len(X), X.shape[1]))

    # 手工特征的基线(上一轮的代理分)
    base = [score(load_pose(os.path.join(args.poses, r["name"] + "_out.pdbqt"))[0],
                  load_pose(os.path.join(args.poses, r["name"] + "_out.pdbqt"))[1],
                  rec_xyz, rec_e)["score"] for r in rows[:len(X)]]
    rho_base = _spearman(base, y)
    print("\n基线(手工几何分): Spearman = %+.3f" % rho_base)

    # K 折交叉验证: 只用折外预测算指标
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(X))
    folds = np.array_split(idx, args.folds)
    pred = np.zeros(len(X))
    for i, te in enumerate(folds):
        tr = np.concatenate([f for j, f in enumerate(folds) if j != i])
        mu, sd = X[tr].mean(axis=0), X[tr].std(axis=0) + 1e-9
        w = ridge_fit((X[tr] - mu) / sd, y[tr], args.lam)
        pred[te] = np.hstack([(X[te] - mu) / sd, np.ones((len(te), 1))]) @ w
    rho = _spearman(pred, y)
    print("学习代理(ridge, %d 折折外): Spearman = %+.3f" % (args.folds, rho))
    print("提升: %+.3f" % (rho - rho_base))
    print()
    if abs(rho) >= 0.30:
        print("判读: |rho| >= 0.30 —— 可作为**筛选器**的候选, 但仍需在")
        print("      模型自己生成的位姿上复验分布迁移性, 且**不得替代 A/B 主指标**。")
    else:
        print("判读: |rho| < 0.30 —— **仍不可用**。按纪律如实判失败, 不调参硬凑。")
        print("      下一步方向: 特征需更贴近 Vina 打分函数(静电/去溶剂化/扭转熵),")
        print("      或直接在生成位姿上做监督, 而不是在 Vina 优化过的位姿上。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
