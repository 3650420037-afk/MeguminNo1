# -*- coding: utf-8 -*-
"""CPU 集成预检：V2（逐样本平衡档）配置在**真实训练路径**上能否跑通、损失是否有限。

为什么要有它：单元测试只验证了 pos_weight 的计算函数；而 V2 是要花 ~1.5 h GPU 的训练。
本预检在**纯 CPU**上把"配置 → transform 链 → 数据集 → 模型 → get_loss 一次前向"整条路走一遍，
把"配置写错/维度不匹配/损失出 NaN"这类问题在开训前暴露（本项目已真实踩过
"transform 在 iter 43 才崩，白跑一次训练"的坑，故 train.py 里也有 --verify-dataset 门禁）。

跑法（CPU，不占 GPU）：
    python tests/test_frontier_balance_integration.py
    # 也可指定别的配置：
    python tests/test_frontier_balance_integration.py configs/train_gpcr_mass_v10term.yml
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import torch  # noqa: E402
import yaml  # noqa: E402
from easydict import EasyDict  # noqa: E402

from models.maskfill import MaskFillModelVN  # noqa: E402
from utils.datasets import get_dataset  # noqa: E402
from utils.transforms import (  # noqa: E402
    AtomComposer, Compose, ContrastiveSample, EdgeSample, FeaturizeLigandAtom,
    FeaturizeProteinAtom, FocalBuilder, LigandCountNeighbors, RefineData, get_mask)


def build(config):
    pf, lf = FeaturizeProteinAtom(), FeaturizeLigandAtom()
    composer = AtomComposer(pf.feature_dim, lf.feature_dim, config.model.encoder.knn)
    c = config.train.transform.contrastive
    es = EdgeSample(config.train.transform.edgesampler)
    cs = ContrastiveSample(c.num_real, c.num_fake, c.pos_real_std, c.pos_fake_std,
                           config.model.field.knn)
    tr = Compose([RefineData(), LigandCountNeighbors(), pf, lf,
                  get_mask(config.train.transform.mask), composer, FocalBuilder(), es, cs])
    dims = (cs.num_elements, es.num_bond_types, pf.feature_dim, lf.feature_dim)
    return tr, dims


def run_one(config_path, split="test", idx=0):
    cfg_path = os.path.join(ROOT, config_path)
    with open(cfg_path, encoding="utf-8") as fh:
        config = EasyDict(yaml.safe_load(fh))
    tr, dims = build(config)
    dataset, subsets = get_dataset(config=config.dataset, transform=tr)
    ds = subsets[split]
    batch = ds[idx]

    model = MaskFillModelVN(config.model, *dims)
    model.eval()
    captured = {}
    orig = model._frontier_pos_weight

    def spy(target, device, dtype):
        w = orig(target, device, dtype)
        captured["pos_weight"] = None if w is None else float(w.reshape(-1)[0])
        captured["label_rate"] = float(target.mean())
        return w

    model._frontier_pos_weight = spy
    with torch.no_grad():
        out = model.get_loss(
            pos_real=batch.pos_real, y_real=batch.cls_real.long(), pos_fake=batch.pos_fake,
            edge_index_real=torch.stack([batch.real_compose_edge_index_0,
                                         batch.real_compose_edge_index_1], dim=0),
            edge_label=batch.real_compose_edge_type,
            index_real_cps_edge_for_atten=batch.index_real_cps_edge_for_atten,
            tri_edge_index=batch.tri_edge_index, tri_edge_feat=batch.tri_edge_feat,
            compose_feature=batch.compose_feature.float(), compose_pos=batch.compose_pos,
            idx_ligand=batch.idx_ligand_ctx_in_compose,
            idx_protein=batch.idx_protein_in_compose,
            y_frontier=batch.ligand_frontier, idx_focal=batch.idx_focal_in_compose,
            pos_generate=batch.pos_generate, idx_protein_all_mask=batch.idx_protein_all_mask,
            y_protein_frontier=batch.y_protein_frontier,
            compose_knn_edge_index=batch.compose_knn_edge_index,
            compose_knn_edge_feature=batch.compose_knn_edge_feature,
            real_compose_knn_edge_index=torch.stack([batch.real_compose_knn_edge_index_0,
                                                     batch.real_compose_knn_edge_index_1], dim=0),
            fake_compose_knn_edge_index=torch.stack([batch.fake_compose_knn_edge_index_0,
                                                     batch.fake_compose_knn_edge_index_1], dim=0))
    model._frontier_pos_weight = orig
    loss, loss_frontier = out[0], out[1]
    assert torch.isfinite(loss).all(), "总损失非有限：%s" % loss
    assert torch.isfinite(loss_frontier).all(), "frontier 损失非有限：%s" % loss_frontier
    return config, captured, float(loss), float(loss_frontier)


def main():
    paths = sys.argv[1:] or ["configs/train_gpcr_mass_v10term_bal.yml",
                             "configs/train_gpcr_mass_v10term.yml"]
    fails = 0
    for p in paths:
        config, cap, loss, lf = run_one(p)
        fb = bool(config.model.get("frontier_balance", False))
        pw = config.model.get("frontier_pos_weight", 1.0)
        print("%-46s frontier_balance=%-5s pos_weight=%-4s → 实算 pos_weight=%s "
              "标签正类率=%.3f 总损失=%.4f frontier损失=%.4f"
              % (os.path.basename(p), fb, pw, cap.get("pos_weight"),
                 cap.get("label_rate", float("nan")), loss, lf))
        # 期望：bal 档必须真的用上 >1 的权重；常数档按配置；默认档为 None
        got = cap.get("pos_weight")
        if fb:
            ok = got is not None and got >= 1.0
        elif abs(float(pw) - 1.0) < 1e-9:
            ok = got is None
        else:
            ok = got is not None and abs(got - float(pw)) < 1e-9
        if not ok:
            print("   !! 期望与实际不符：frontier_balance=%s, 配置 pos_weight=%s, 实算 %s" % (fb, pw, got))
            fails += 1
        else:
            print("   OK 权重路径符合配置语义")
    print("\nRESULT:", "ALL OK" if fails == 0 else "%d 项不符" % fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
