# -*- coding: utf-8 -*-
"""**纯 CPU** 探针：某个检查点的 frontier 头在验证集上"预测还会生长"的比例。

为什么需要它（§5.33 处方②的触发量很贵）：
处方②要求"每 250 步记录 `frac(HA≥28)`，连续两点 =0 视为停早失效"——
但 `frac(HA≥28)` 要**采样**（GPU，每点几分钟）。本探针用**一次前向**（无采样）算
`frac_pred_frontier = mean(sigmoid(frontier_logit) > 0)` 与平均概率，
检验它能否当作"停早模式"的**便宜代理量**（先在已知尺寸的检查点上做效度检验）。

做法：复用**训练时的同一批 transform 链**与 `model.get_loss(...)` 调用路径，
在 `frontier_pred` 上挂 forward hook 抓取**配体 frontier** 那一次的输出
（该模块在一次前向里被调用两次：先蛋白质表面、后配体上下文 → 取最后一次）。
**不采样、不建库、不碰 GPU**（`CUDA_VISIBLE_DEVICES` 为空时也能跑）。

用法：
    python src/scripts/probe_frontier_rate.py --ckpts <a.pt> <b.pt> --n 8 \\
        --out outputs/probe_frontier_rate.csv
"""
import argparse
import csv
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np  # noqa: E402
import torch  # noqa: E402
from easydict import EasyDict  # noqa: E402
import yaml  # noqa: E402

from paths import ROOT  # noqa: E402
from models.maskfill import MaskFillModelVN  # noqa: E402
from utils.datasets import get_dataset  # noqa: E402
from utils.transforms import (  # noqa: E402
    AtomComposer, Compose, ContrastiveSample, EdgeSample, FeaturizeLigandAtom,
    FeaturizeProteinAtom, FocalBuilder, LigandCountNeighbors, RefineData, get_mask)


def build_transform(config):
    pf = FeaturizeProteinAtom()
    lf = FeaturizeLigandAtom()
    composer = AtomComposer(pf.feature_dim, lf.feature_dim, config.model.encoder.knn)
    c = config.train.transform.contrastive
    return Compose([RefineData(), LigandCountNeighbors(), pf, lf,
                    get_mask(config.train.transform.mask), composer, FocalBuilder(),
                    EdgeSample(config.train.transform.edgesampler),
                    ContrastiveSample(c.num_real, c.num_fake, c.pos_real_std, c.pos_fake_std,
                                      config.model.field.knn)])


def load_config(ckpt):
    cfg = torch.load(ckpt, map_location="cpu", weights_only=False)["config"]
    return EasyDict(cfg)


def probe(ckpt, dataset, n, dims, mask_seed=1234, target_field="ligand_frontier"):
    """返回 (frac_pred>0, mean_prob, label_frac, n_atoms, n_used) 在该检查点上的平均。

    dims = (num_classes, num_bond_types, protein_dim, ligand_dim) —— 与 train.py 的
    MaskFillModelVN 构造方式一致（从同一套 featurizer/edge sampler 取维度）。
    """
    config = load_config(ckpt)
    model = MaskFillModelVN(config.model, *dims)
    sd = torch.load(ckpt, map_location="cpu", weights_only=False)["model"]
    missing, unexpected = model.load_state_dict(sd, strict=False)
    if unexpected:
        raise SystemExit("权重与模型不匹配: unexpected=%s" % list(unexpected)[:5])
    model.eval()

    captured = []

    def hook(_mod, _inp, out):
        captured.append(out.detach())

    h = model.frontier_pred.register_forward_hook(hook)
    fracs, probs, maxes, labels, n_atoms_total = [], [], [], [], 0
    try:
        with torch.no_grad():
            for i in range(min(n, len(dataset))):
                captured.clear()
                # **固定掩码**：不同检查点必须看到完全相同的掩码，否则"标签正类率"本身就不同,
                # 预测率不可比（首次实现踩过：39 vs 65 原子）。以样本序号为种子 → 可复现。
                torch.manual_seed(mask_seed + i)
                np.random.seed((mask_seed + i) % (2 ** 32))
                random.seed(mask_seed + i)
                batch = dataset[i]
                try:
                    model.get_loss(
                        pos_real=batch.pos_real, y_real=batch.cls_real.long(),
                        pos_fake=batch.pos_fake,
                        edge_index_real=torch.stack([batch.real_compose_edge_index_0,
                                                     batch.real_compose_edge_index_1], dim=0),
                        edge_label=batch.real_compose_edge_type,
                        index_real_cps_edge_for_atten=batch.index_real_cps_edge_for_atten,
                        tri_edge_index=batch.tri_edge_index, tri_edge_feat=batch.tri_edge_feat,
                        compose_feature=batch.compose_feature.float(),
                        compose_pos=batch.compose_pos,
                        idx_ligand=batch.idx_ligand_ctx_in_compose,
                        idx_protein=batch.idx_protein_in_compose,
                        y_frontier=batch.ligand_frontier,
                        idx_focal=batch.idx_focal_in_compose,
                        pos_generate=batch.pos_generate,
                        idx_protein_all_mask=batch.idx_protein_all_mask,
                        y_protein_frontier=batch.y_protein_frontier,
                        compose_knn_edge_index=batch.compose_knn_edge_index,
                        compose_knn_edge_feature=batch.compose_knn_edge_feature,
                        real_compose_knn_edge_index=torch.stack(
                            [batch.real_compose_knn_edge_index_0,
                             batch.real_compose_knn_edge_index_1], dim=0),
                        fake_compose_knn_edge_index=torch.stack(
                            [batch.fake_compose_knn_edge_index_0,
                             batch.fake_compose_knn_edge_index_1], dim=0))
                except Exception as e:
                    print("  样本 %d 失败: %s" % (i, str(e)[:90]))
                    continue
                if not captured:
                    continue
                logits = captured[-1].reshape(-1)          # 最后一次 = 配体 frontier
                p = torch.sigmoid(logits)
                fracs.append(float((p > 0.5).float().mean()))
                probs.append(float(p.mean()))
                maxes.append(float(p.max()))
                lab = batch[target_field].reshape(-1).float()
                labels.append(float(lab.mean()))
                n_atoms_total += int(lab.numel())
    finally:
        h.remove()
    if not fracs:
        return None
    has_frontier = sum(1 for m in maxes if m > 0.5) / float(len(maxes))
    return (sum(fracs) / len(fracs), sum(probs) / len(probs),
            sum(labels) / len(labels), n_atoms_total, len(fracs),
            sum(maxes) / len(maxes), has_frontier)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpts", nargs="+", required=True)
    ap.add_argument("--config", default="configs/train_gpcr_mass_v10term.yml",
                    help="只用于取数据集与 transform 定义（模型结构取自各检查点自身 config）")
    ap.add_argument("--n", type=int, default=8, help="验证集样本数")
    ap.add_argument("--split", default="test", help="用验证划分（train.py 里是 'test'）")
    ap.add_argument("--mask-seed", type=int, default=1234,
                    help="固定掩码种子：所有检查点必须用同一个值才可比")
    ap.add_argument("--out", default=os.path.join("outputs", "probe_frontier_rate.csv"))
    args = ap.parse_args()

    with open(os.path.join(ROOT, args.config), encoding="utf-8") as fh:
        config = EasyDict(yaml.safe_load(fh))
    transform = build_transform(config)
    # 维度与 train.py 一致：从同一套 featurizer / edge sampler 取
    pf = [t for t in transform.transforms if isinstance(t, FeaturizeProteinAtom)][0]
    lf = [t for t in transform.transforms if isinstance(t, FeaturizeLigandAtom)][0]
    es = [t for t in transform.transforms if isinstance(t, EdgeSample)][0]
    cs = [t for t in transform.transforms if isinstance(t, ContrastiveSample)][0]
    dims = (cs.num_elements, es.num_bond_types, pf.feature_dim, lf.feature_dim)
    dataset, subsets = get_dataset(config=config.dataset, transform=transform)
    ds = subsets[args.split]
    print("验证划分 %s: %d 个样本；每点取 %d 个" % (args.split, len(ds), args.n))

    rows = []
    for ck in args.ckpts:
        r = probe(ck, ds, args.n, dims, args.mask_seed)
        if r is None:
            print("%-28s 失败（无有效样本）" % os.path.basename(ck))
            continue
        frac_pred, mean_prob, label_frac, n_atoms, n_used, mean_max, has_fr = r
        print("%-28s frac_pred %.3f | 平均概率 %.3f | **每分子最大概率均值 %.3f** | "
              "有 frontier 的分子占比 %.2f | 标签正类率 %.3f | %d 原子/%d 样本"
              % (os.path.basename(ck), frac_pred, mean_prob, mean_max, has_fr,
                 label_frac, n_atoms, n_used))
        rows.append(dict(ckpt=os.path.abspath(ck), name=os.path.basename(ck),
                         frac_pred_frontier=round(frac_pred, 4),
                         mean_prob=round(mean_prob, 4),
                         mean_per_mol_max=round(mean_max, 4),
                         frac_mol_has_frontier=round(has_fr, 4),
                         label_frac=round(label_frac, 4),
                         n_atoms=n_atoms, n_samples=n_used))
    if rows:
        os.makedirs(os.path.dirname(os.path.join(ROOT, args.out)), exist_ok=True)
        with open(os.path.join(ROOT, args.out), "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print("明细: %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
