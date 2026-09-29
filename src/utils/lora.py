# -*- coding: utf-8 -*-
"""LoRA / Adapter：参数高效微调（用户路线③），**默认关闭**，关闭时与原模型逐位一致。

为什么需要（对齐目标与已有证据）
----------------------------
- v3/v10/v11 都是**全量微调**：编码器（2.984 M，占 80.4 %）与头部一起动，容易把官方权重
  已有的"尺寸/终止先验"冲掉（§5.33 的尺寸漂移即发生在这种训练里）。
- LoRA 只训练低秩增量 `ΔW = B·A · α/r`，**基座权重全程冻结** —— 理论上更能保住官方先验，
  且训练参数量小 1–2 个数量级，适合 8 GB 显存与短预算。
- 本文件**不改变任何既有行为**：不调用 `apply_lora` 时模型与原来逐位相同（有测试保证）。

设计要点
--------
1. `LoRALinear` 包住既有 `nn.Linear`；`B` **零初始化** → 初始 `ΔW = 0`，
   包装前后前向输出**逐位相同**（这是"可回退/可对照"的前提，见 tests/test_lora.py）。
2. `apply_lora(model, targets=('encoder',), rank, alpha)` 按**模块名前缀**注入；
   基座参数 `requires_grad_(False)`，只有 `lora_A/lora_B` 可训练。
3. `lora_parameters(model)` 返回适配器参数（用于单独设学习率/优化器）。
4. `merge_lora(model)` 把适配器折进基座权重并**还原成普通 nn.Linear** —— 导出的权重文件
   不需要任何 LoRA 代码即可加载（对交付/可比性很重要）。
5. `lora_state_dict(model)` / `load_lora_state_dict(model, sd)` 便于只存取适配器（几 MB）。
"""
import torch
import torch.nn as nn


class LoRALinear(nn.Module):
    """把 nn.Linear 包成 base + 低秩增量；B 零初始化 ⇒ 初始输出与 base 逐位相同。"""

    def __init__(self, base: nn.Linear, rank: int = 8, alpha: float = 16.0):
        super().__init__()
        if not isinstance(base, nn.Linear):
            raise TypeError('LoRALinear 只能包 nn.Linear, 收到 %s' % type(base))
        if rank <= 0:
            raise ValueError('rank 必须 > 0')
        self.base = base
        for p in self.base.parameters():
            p.requires_grad_(False)
        self.rank = int(rank)
        self.alpha = float(alpha)
        self.scaling = self.alpha / self.rank
        self.lora_A = nn.Linear(base.in_features, self.rank, bias=False)
        self.lora_B = nn.Linear(self.rank, base.out_features, bias=False)
        nn.init.kaiming_uniform_(self.lora_A.weight, a=5 ** 0.5)
        nn.init.zeros_(self.lora_B.weight)          # ← 零初始化：初始 ΔW = 0
        # ⚠ 必须与 base 同 device/dtype：`apply_lora` 通常在模型 `.to(device)` **之后**调用，
        # 而 nn.Linear 默认建在 CPU —— 不搬的话首次前向必然
        # `RuntimeError: Expected all tensors to be on the same device`（独立代码审查实测发现，
        # 用 meta 设备可复现：base=meta、适配器=cpu）。
        self.lora_A.to(device=base.weight.device, dtype=base.weight.dtype)
        self.lora_B.to(device=base.weight.device, dtype=base.weight.dtype)

    def forward(self, x):
        return self.base(x) + self.lora_B(self.lora_A(x)) * self.scaling

    @torch.no_grad()
    def merged_weight(self):
        """返回折进适配器后的权重（不修改自身）。"""
        delta = (self.lora_B.weight @ self.lora_A.weight) * self.scaling
        return self.base.weight + delta

    def extra_repr(self):
        return 'in=%d, out=%d, rank=%d, alpha=%.1f, scaling=%.3f' % (
            self.base.in_features, self.base.out_features, self.rank, self.alpha, self.scaling)


def _qualified_names(model, targets):
    """返回所有位于 targets 前缀下、类型为 nn.Linear 的 (qualname, module)。"""
    out = []
    prefixes = tuple(targets or ())
    for name, mod in model.named_modules():
        if not name:
            continue
        if not isinstance(mod, nn.Linear):
            continue
        if prefixes and not name.startswith(prefixes):
            continue
        out.append((name, mod))
    return out


def _set_submodule(model, qualname, new_mod):
    parent = model
    parts = qualname.split('.')
    for p in parts[:-1]:
        parent = getattr(parent, p) if not p.isdigit() else parent[int(p)]
    setattr(parent, parts[-1], new_mod) if not parts[-1].isdigit() \
        else parent.__setitem__(int(parts[-1]), new_mod)


def apply_lora(model, targets=('encoder',), rank=8, alpha=16.0, verbose=True):
    """按模块名前缀注入 LoRA。返回被替换的层数。**基座冻结、只有适配器可训练。**

    ⚠ 若同一个 `nn.Linear` **实例**被多处引用（权重共享），本函数会**显式报错**：
    `merge_lora` 用的是去掉重复引用的 `named_modules()`，共享实例只会被折入一次，
    另一处引用仍指向旧包装器 → 会得到"看起来合并了、其实不一致"的模型。
    本仓库实测 `MaskFillModelVN` 无共享 Linear（262 实例 / 262 引用），但宁可失败也不要静默出错。
    """
    hits = _qualified_names(model, targets)
    if not hits:
        raise ValueError('没有匹配到任何 nn.Linear（targets=%s）—— 请检查前缀' % (targets,))
    seen_shared = {}
    for name, mod in model.named_modules(remove_duplicate=False):
        if isinstance(mod, nn.Linear):
            seen_shared.setdefault(id(mod), []).append(name)
    shared = {k: v for k, v in seen_shared.items() if len(v) > 1}
    hit_names = {n for n, _ in hits}
    bad = {k: v for k, v in shared.items() if any(n in hit_names for n in v)}
    if bad:
        raise NotImplementedError(
            '检测到模块级权重共享的 nn.Linear，apply_lora/merge_lora 不支持该情形（会静默漏合并）: %s'
            % list(bad.values())[:2])
    # **参数级**共享（两个不同 Linear 的 weight 是同一张量，如 `b.weight = a.weight`）也要拦：
    # merge_lora 的 `base.weight.copy_()` 是原地写，会把增量叠加两次（独立审查实测差 0.46）。
    pseen = {}
    for name, mod in model.named_modules(remove_duplicate=False):
        if isinstance(mod, nn.Linear) and name in hit_names:
            for pname, p in (('weight', mod.weight), ('bias', mod.bias)):
                if p is not None:
                    pseen.setdefault(id(p), []).append('%s.%s' % (name, pname))
    pshared = {k: v for k, v in pseen.items() if len(v) > 1}
    if pshared:
        raise NotImplementedError(
            '检测到参数级权重共享（同一张量被多个 Linear 使用），merge 会重复叠加: %s'
            % list(pshared.values())[:2])
    for name, mod in hits:
        _set_submodule(model, name, LoRALinear(mod, rank=rank, alpha=alpha))
    # 冻结所有非 LoRA 参数（含头部）；若要"适配器 + 头部"一起训，请随后自行解冻头部
    for name, p in model.named_parameters():
        p.requires_grad_('lora_' in name)
    if verbose:
        print('[LoRA] 注入 %d 层（targets=%s, rank=%d, alpha=%.1f）' % (
            len(hits), ','.join(targets), rank, alpha))
    return len(hits)


def unfreeze_patterns(model, patterns):
    """把名字含任一 pattern 的参数设为可训练（例如 'frontier_pred'）——适配器之外再解冻头部。"""
    n = 0
    for name, p in model.named_parameters():
        if any(pat in name for pat in patterns):
            p.requires_grad_(True)
            n += p.numel()
    return n


def lora_parameters(model):
    """返回可训练的 LoRA 适配器参数（用于单独给学习率）。"""
    return [p for name, p in model.named_parameters()
            if ('lora_A' in name or 'lora_B' in name) and p.requires_grad]


def lora_state_dict(model):
    """只含适配器的 state_dict（几 MB 量级，便于单独保存/比较）。"""
    return {k: v.detach().clone() for k, v in model.state_dict().items()
            if 'lora_A' in k or 'lora_B' in k}


def load_lora_state_dict(model, sd):
    """把适配器权重灌回（键必须完全匹配当前 LoRA 结构）。

    ⚠ 若模型**根本没有** LoRA 层（未注入），旧实现会"什么都不加载却返回成功计数"——
    独立审查实测到这一静默失败，现改为显式报错。
    """
    have = [k for k in model.state_dict() if 'lora_A' in k or 'lora_B' in k]
    if not have:
        raise RuntimeError('该模型没有 LoRA 层（未调用 apply_lora），load_lora_state_dict 无意义；'
                           '请先注入适配器')
    missing = [k for k in have if k not in sd]
    if missing:
        raise KeyError('适配器键缺失 %d 个，例如 %s' % (len(missing), missing[:3]))
    model.load_state_dict(sd, strict=False)
    return len(sd)


@torch.no_grad()
def merge_lora(model, trainable=False):
    """把 LoRA 折进基座并**还原为普通 nn.Linear**（导出的权重不再依赖本模块）。返回折入层数。

    trainable=False（默认）：折入后基座参数保持 `requires_grad=False`（导出/采样用，无需梯度）。
    trainable=True：折入后把**整个模型**的参数设为可训练（"合并后继续全量微调"）。
    注意：旧实现只解冻被包装的 base，encoder 内非 Linear 参数、field/pos_predictor 等仍是冻结的
    （独立审查实测只有 81.9% 可训练），与"全量微调"的语义不符 —— 现已改为全模型解冻。
    """
    hits = [(n, m) for n, m in model.named_modules() if isinstance(m, LoRALinear)]
    for name, mod in hits:
        base = mod.base
        base.weight.copy_(mod.merged_weight())
        _set_submodule(model, name, base)
    if trainable:
        for p in model.parameters():
            p.requires_grad_(True)
    return len(hits)


def count_trainable(model):
    tr = sum(p.numel() for p in model.parameters() if p.requires_grad)
    tot = sum(p.numel() for p in model.parameters())
    return tr, tot
