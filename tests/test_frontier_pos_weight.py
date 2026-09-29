# -*- coding: utf-8 -*-
"""CPU 单测：frontier 损失三档 pos_weight（默认档必须与官方逐位一致）。

跑法（CPU，不占 GPU）：
    python -m pytest tests/test_frontier_pos_weight.py -q
或直接：
    python tests/test_frontier_pos_weight.py
"""
import os
import sys

try:  # 无人值守/重定向下 stdout 可能是 GBK 管道：非 GBK 字符会直接崩（本项目已真实踩过）
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import types

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import torch  # noqa: E402
import torch.nn.functional as F  # noqa: E402

from models.maskfill import MaskFillModelVN  # noqa: E402


def _stub(**kw):
    """构造只带这些属性的假模型（方法按未绑定方式调用）。"""
    return types.SimpleNamespace(**kw)


def compute(target, **attrs):
    return MaskFillModelVN._frontier_pos_weight(_stub(**attrs), target,
                                                 torch.device("cpu"), torch.float32)


def test_default_branch_is_none():
    """① 默认配置 → None，且与官方 BCE 逐位相同。"""
    torch.manual_seed(0)
    logits = torch.randn(64, 1)
    target = (torch.rand(64, 1) > 0.7).float()
    w = compute(target, frontier_pos_weight=1.0)
    assert w is None, w
    a = F.binary_cross_entropy_with_logits(logits, target, pos_weight=w)
    b = F.binary_cross_entropy_with_logits(logits, target)          # 官方调用形式
    assert torch.equal(a, b), (a, b)
    print("PASS ① 默认档 pos_weight=None，与官方 BCE 逐位相同")


def test_constant_branch():
    """② 常数档 → 与显式 pos_weight=2.0 逐位相同。"""
    torch.manual_seed(1)
    logits = torch.randn(64, 1)
    target = (torch.rand(64, 1) > 0.7).float()
    w = compute(target, frontier_pos_weight=2.0)
    assert w is not None and w.shape == (1,) and abs(float(w) - 2.0) < 1e-12, w
    a = F.binary_cross_entropy_with_logits(logits, target, pos_weight=w)
    b = F.binary_cross_entropy_with_logits(logits, target, pos_weight=torch.tensor([2.0]))
    assert torch.equal(a, b), (a, b)
    print("PASS ② 常数档 = 2.0，与显式写法逐位相同")


def test_balance_branch_math():
    """③ 逐样本平衡档 → n_neg/n_pos 且被 clip。"""
    # 20 正 / 80 负 → 权重 4.0
    target = torch.cat([torch.ones(20, 1), torch.zeros(80, 1)])
    w = compute(target, frontier_balance=True, frontier_balance_clip=(1.0, 10.0))
    assert abs(float(w) - 4.0) < 1e-6, w
    # 极端不平衡 1 正 / 99 负 → 99 → 被 clip 到 10
    target2 = torch.cat([torch.ones(1, 1), torch.zeros(99, 1)])
    w2 = compute(target2, frontier_balance=True, frontier_balance_clip=(1.0, 10.0))
    assert abs(float(w2) - 10.0) < 1e-6, w2
    # 全正 → 权重被 clip 到下限 1.0（且不为 0/NaN）
    w3 = compute(torch.ones(30, 1), frontier_balance=True, frontier_balance_clip=(1.0, 10.0))
    assert abs(float(w3) - 1.0) < 1e-6, w3
    # 全负 → n_pos 被 clamp 到 1 → 权重 = n_neg → 再被 clip
    w4 = compute(torch.zeros(50, 1), frontier_balance=True, frontier_balance_clip=(1.0, 10.0))
    assert abs(float(w4) - 10.0) < 1e-6, w4
    print("PASS ③ 平衡档 n_neg/n_pos 与 clip 行为正确（含全正/全负边界）")


def test_balance_overrides_constant():
    """③ 优先于 ②：同时给两个键时按逐样本平衡算。"""
    target = torch.cat([torch.ones(10, 1), torch.zeros(30, 1)])
    w = compute(target, frontier_balance=True, frontier_balance_clip=(1.0, 10.0),
                frontier_pos_weight=2.0)
    assert abs(float(w) - 3.0) < 1e-6, w
    print("PASS ③ 平衡档覆盖常数档")


def test_default_attr_absent():
    """属性完全缺失（老权重/老配置）→ 仍返回 None（向后兼容）。"""
    target = torch.cat([torch.ones(5, 1), torch.zeros(5, 1)])
    assert compute(target) is None
    print("PASS 兼容 缺属性时 default=None（官方行为）")


if __name__ == "__main__":
    test_default_branch_is_none()
    test_constant_branch()
    test_balance_branch_math()
    test_balance_overrides_constant()
    test_default_attr_absent()
    print("ALL PASS (5/5)")
