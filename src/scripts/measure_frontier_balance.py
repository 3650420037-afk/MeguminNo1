# -*- coding: utf-8 -*-
"""量化 mask-fill 训练里 frontier 正类的占比 —— 决定 pos_weight 该给多大。

背景（M8 处方②）：v10term 用 frontier_pos_weight=2.0 抑制"过早停止"。
BCE 的 pos_weight 语义是"正类损失权重"；要把正/负类**拉平**，理论值是 (1-p)/p。
本脚本用**训练完全相同的 transform 链**统计 p（正类 = context 原子中"还有生长空间"的原子）。

用法（CPU 即可，不碰 GPU）：
    python src/scripts/measure_frontier_balance.py --config configs/train_gpcr_mass_v10term.yml --n 400
"""
import argparse
import os
import statistics as st
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch  # noqa: E402
from easydict import EasyDict  # noqa: E402
import yaml  # noqa: E402

from paths import ROOT  # noqa: E402
from utils.datasets import get_dataset  # noqa: E402
from utils.transforms import (  # noqa: E402
    AtomComposer, Compose, ContrastiveSample, EdgeSample, FeaturizeLigandAtom,
    FeaturizeProteinAtom, FocalBuilder, LigandCountNeighbors, RefineData, get_mask)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/train_gpcr_mass_v10term.yml")
    ap.add_argument("--n", type=int, default=400, help="抽样样本数（均匀取自 train split）")
    ap.add_argument("--split", default="train")
    args = ap.parse_args()

    with open(os.path.join(ROOT, args.config), encoding="utf-8") as fh:
        config = EasyDict(yaml.safe_load(fh))

    protein_featurizer = FeaturizeProteinAtom()
    ligand_featurizer = FeaturizeLigandAtom()
    masking = get_mask(config.train.transform.mask)
    composer = AtomComposer(protein_featurizer.feature_dim, ligand_featurizer.feature_dim,
                            config.model.encoder.knn)
    edge_sampler = EdgeSample(config.train.transform.edgesampler)
    cfg_ctr = config.train.transform.contrastive
    contrastive_sampler = ContrastiveSample(cfg_ctr.num_real, cfg_ctr.num_fake,
                                            cfg_ctr.pos_real_std, cfg_ctr.pos_fake_std,
                                            config.model.field.knn)
    transform = Compose([RefineData(), LigandCountNeighbors(), protein_featurizer,
                         ligand_featurizer, masking, composer, FocalBuilder(),
                         edge_sampler, contrastive_sampler])

    ds_root = os.path.join(ROOT, config.dataset.path.lstrip("./"))
    dataset, subsets = get_dataset(config=config.dataset, transform=transform)
    ds = subsets[args.split]
    n = min(args.n, len(ds))
    step = max(1, len(ds) // n)

    rates, n_atoms_all, lig_sizes = [], 0, []
    pos_all = 0
    for i in range(0, len(ds), step):
        try:
            data = ds[i]
        except Exception as e:
            print("跳过样本 %d: %s" % (i, e))
            continue
        fr = data.ligand_frontier
        # frontier 定义在 context(未掩码且被保留的)原子上: 与 y_frontier 监督的目标集一致
        n_ctx = int(fr.numel())
        if n_ctx == 0:
            continue
        pos = int(fr.sum().item())
        rates.append(pos / float(n_ctx))
        pos_all += pos
        n_atoms_all += n_ctx
        lig_sizes.append(int(data.ligand_element.numel()))
        if len(rates) >= n:
            break

    if not rates:
        print("没有取到样本；请检查 dataset.path 与 split")
        return 1
    p = st.mean(rates)
    p_micro = pos_all / float(n_atoms_all)   # BCE 按节点取 mean → 这才是相关量
    q = sorted(rates)
    print("数据集: %s (split=%s)" % (config.dataset.path, args.split))
    print("抽样样本数: %d（配体原子中位 %.0f；context 原子合计 %d）"
          % (len(rates), st.median(lig_sizes), n_atoms_all))
    print("frontier 正类占比：**微平均 p = %.4f**（按节点，BCE 相关量）" % p_micro)
    print("                  逐样本平均 %.4f / 中位 %.4f / 范围 %.3f–%.3f"
          % (p, st.median(rates), min(rates), max(rates)))
    print("  逐样本 p 的分位：p10 %.3f · p25 %.3f · p50 %.3f · p75 %.3f · p90 %.3f"
          % (q[int(0.10 * len(q))], q[int(0.25 * len(q))], q[int(0.50 * len(q))],
             q[int(0.75 * len(q))], q[min(len(q) - 1, int(0.90 * len(q)))]))
    print("→ 拉平正负类的理论 pos_weight = (1-p_micro)/p_micro = %.2f" % ((1 - p_micro) / p_micro))
    print("→ 当前 v10term 用 pos_weight = 2.0")
    for w in (1.0, 2.0, 3.0, (1 - p_micro) / p_micro):
        share = w * p_micro / (w * p_micro + (1 - p_micro))
        print("   pos_weight=%.2f → 正类占总损失权重 %.3f" % (w, share))

    # 逐样本平衡档（frontier_balance=true）在真实数据上会给出的权重分布
    ws = sorted(min(10.0, max(1.0, (1 - r) / r)) for r in rates if r > 0)
    if ws:
        print("   逐样本平衡档权重 clip((1-p)/p, 1, 10)：均值 %.2f · p10 %.2f · p50 %.2f · p90 %.2f"
              % (st.mean(ws), ws[int(0.10 * len(ws))], ws[int(0.50 * len(ws))],
                 ws[min(len(ws) - 1, int(0.90 * len(ws)))]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
