# -*- coding: utf-8 -*-
"""CPU 集成测试：LoRA 在**真实 MaskFillModelVN + 真实 get_loss 路径**上可用且数值透明。

要证明三件事（都写成可复现断言）：
① 单层层面：零初始化适配器的增量 **恰好为 0**（数学恒等）；
② 端到端层面：注入 LoRA 后的 encoder 输出差 ≤1e-5，且**把冻结撤销后差值恰好为 0.000e+00**
   —— 即这点微小差异**来源是 `requires_grad=False` 引起的 PyTorch 内核路径差异，不是 LoRA 的数学误差**
   （2026-09-29 实测定位：纯 deepcopy 差 0；只加 `requires_grad_(False)` 差 7.2e-06；
   加 LoRA 后再解除冻结差 0）。这条结论决定了"LoRA 开关可以做受控对照"是否成立。
③ 训练一步：只有适配器（与显式解冻的头部）改变，其余权重零变化；损失确实变化（梯度真的流经适配器）。

跑法：python tests/test_lora_integration.py
"""
import copy
import os
import sys

try:  # 无人值守/重定向下 stdout 可能是 GBK 管道：非 GBK 字符会直接崩（本项目已真实踩过）
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import torch  # noqa: E402
import yaml  # noqa: E402
from easydict import EasyDict  # noqa: E402

from models.maskfill import MaskFillModelVN  # noqa: E402
from utils.datasets import get_dataset  # noqa: E402
from utils.lora import (LoRALinear, apply_lora, count_trainable,  # noqa: E402
                        lora_parameters, merge_lora, unfreeze_patterns)
from utils.transforms import (  # noqa: E402
    AtomComposer, Compose, ContrastiveSample, EdgeSample, FeaturizeLigandAtom,
    FeaturizeProteinAtom, FocalBuilder, LigandCountNeighbors, RefineData, get_mask)

CFG = os.path.join(ROOT, "configs", "train_gpcr_mass_v10term.yml")
KERNEL_TOL = 1e-5      # requires_grad=False 造成的内核路径差异上界（实测 ~5e-06）


def build_all():
    with open(CFG, encoding="utf-8") as fh:
        config = EasyDict(yaml.safe_load(fh))
    pf, lf = FeaturizeProteinAtom(), FeaturizeLigandAtom()
    composer = AtomComposer(pf.feature_dim, lf.feature_dim, config.model.encoder.knn)
    c = config.train.transform.contrastive
    es = EdgeSample(config.train.transform.edgesampler)
    cs = ContrastiveSample(c.num_real, c.num_fake, c.pos_real_std, c.pos_fake_std,
                           config.model.field.knn)
    tr = Compose([RefineData(), LigandCountNeighbors(), pf, lf,
                  get_mask(config.train.transform.mask), composer, FocalBuilder(), es, cs])
    dims = (cs.num_elements, es.num_bond_types, pf.feature_dim, lf.feature_dim)
    dataset, subsets = get_dataset(config=config.dataset, transform=tr)
    return config, MaskFillModelVN(config.model, *dims), subsets["test"]


def loss_kwargs(batch):
    return dict(
        pos_real=batch.pos_real, y_real=batch.cls_real.long(), pos_fake=batch.pos_fake,
        edge_index_real=torch.stack([batch.real_compose_edge_index_0,
                                     batch.real_compose_edge_index_1], dim=0),
        edge_label=batch.real_compose_edge_type,
        index_real_cps_edge_for_atten=batch.index_real_cps_edge_for_atten,
        tri_edge_index=batch.tri_edge_index, tri_edge_feat=batch.tri_edge_feat,
        compose_feature=batch.compose_feature.float(), compose_pos=batch.compose_pos,
        idx_ligand=batch.idx_ligand_ctx_in_compose, idx_protein=batch.idx_protein_in_compose,
        y_frontier=batch.ligand_frontier, idx_focal=batch.idx_focal_in_compose,
        pos_generate=batch.pos_generate, idx_protein_all_mask=batch.idx_protein_all_mask,
        y_protein_frontier=batch.y_protein_frontier,
        compose_knn_edge_index=batch.compose_knn_edge_index,
        compose_knn_edge_feature=batch.compose_knn_edge_feature,
        real_compose_knn_edge_index=torch.stack([batch.real_compose_knn_edge_index_0,
                                                 batch.real_compose_knn_edge_index_1], dim=0),
        fake_compose_knn_edge_index=torch.stack([batch.fake_compose_knn_edge_index_0,
                                                 batch.fake_compose_knn_edge_index_1], dim=0))


def enc_out(model, batch, seed=0):
    box = []
    h = model.encoder.register_forward_hook(
        lambda _m, _i, o: box.append((o[0] if isinstance(o, (tuple, list)) else o).detach().clone()))
    try:
        model.eval()
        torch.manual_seed(seed)
        with torch.no_grad():
            model.get_loss(**loss_kwargs(batch))
    finally:
        h.remove()
    return box[0]


def main():
    config, base, ds = build_all()
    batch = ds[0]
    a_base = enc_out(base, batch)

    m = copy.deepcopy(base)
    n = apply_lora(m, targets=("encoder",), rank=4, alpha=8.0, verbose=False)
    tr, tot = count_trainable(m)
    print("LoRA 注入 %d 层 | 可训练 %d / %d (%.2f%%)" % (n, tr, tot, 100.0 * tr / tot))

    # ① 单层：零初始化 -> 增量恰好 0
    lin = [mod for mod in m.modules() if isinstance(mod, LoRALinear)][0]
    x = torch.randn(7, lin.base.in_features)
    with torch.no_grad():
        delta = (lin(x) - lin.base(x)).abs().max().item()
    assert delta == 0.0, delta
    print("PASS ① 单层零初始化增量恰好为 0（数学恒等）")

    # ② 端到端：差异 ≤ 内核容差；**解除冻结后必须恰好为 0** ← 决定性判据
    a_lora = enc_out(m, batch)
    d_frozen = (a_base - a_lora).abs().max().item()
    assert d_frozen <= KERNEL_TOL, d_frozen
    for nm, p in m.named_parameters():
        if 'lora_' not in nm:
            p.requires_grad_(True)
    a_unfrozen = enc_out(m, batch)
    d_unfrozen = (a_base - a_unfrozen).abs().max().item()
    assert d_unfrozen == 0.0, d_unfrozen
    print("PASS ② 端到端：冻结态差 %.2e（≤%.0e），**解除冻结后差 %.2e**"
          " -> 差异源于 requires_grad 的内核路径，不是 LoRA 数学误差" % (d_frozen, KERNEL_TOL, d_unfrozen))

    # ③ 训练一步：只改适配器与显式解冻的头部
    m2 = copy.deepcopy(base)
    apply_lora(m2, targets=("encoder",), rank=4, alpha=8.0, verbose=False)
    unfreeze_patterns(m2, ["frontier_pred"])
    snap = {k: v.detach().clone() for k, v in m2.state_dict().items()
            if "lora_" not in k and not k.startswith("frontier_pred")}
    before = enc_out(m2, batch)
    m2.train()
    opt = torch.optim.SGD([p for p in m2.parameters() if p.requires_grad], lr=1e-3)
    torch.manual_seed(0)
    loss = m2.get_loss(**loss_kwargs(batch))[0]
    opt.zero_grad()
    loss.backward()
    opt.step()
    changed = [k for k, v in m2.state_dict().items()
               if k in snap and not torch.equal(v, snap[k])]
    assert not changed, ("冻结部分不应被更新", changed[:5])
    after = enc_out(m2, batch)
    d_step = (before - after).abs().max().item()
    assert d_step > 0, "训练一步后 encoder 输出应变化（否则梯度没流到适配器）"
    assert any(float(p.detach().abs().max()) > 0 for p in lora_parameters(m2)), \
        "适配器在一步 SGD 后应当非零"
    print("PASS ③ 一步 SGD：冻结部分零变化（%d 键），encoder 输出变化 %.3e，适配器已更新"
          % (len(snap), d_step))

    n2 = merge_lora(m2)
    assert not any(isinstance(mod, LoRALinear) for mod in m2.modules())
    print("PASS ④ merge_lora 折入 %d 层，导出权重不再依赖 LoRA 代码" % n2)
    print("\nALL PASS (4/4)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
