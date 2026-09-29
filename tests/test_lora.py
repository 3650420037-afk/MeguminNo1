# -*- coding: utf-8 -*-
"""CPU 单测：LoRA 注入的正确性与"关闭/合并后逐位一致"。

关键不变量（决定它能不能用于受控实验）：
① 注入后 **初始前向输出与基座模型逐位相同**（B 零初始化）；
② 基座权重被冻结、只有适配器可训练；
③ `merge_lora` 之后**输出与 LoRA 模型一致**（在浮点容差内），且模型结构还原为普通 nn.Linear；
④ 适配器 state_dict 可存取往返。

跑法：python tests/test_lora.py
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
import torch.nn as nn  # noqa: E402

from utils.lora import (LoRALinear, apply_lora, count_trainable, load_lora_state_dict,  # noqa: E402
                        lora_parameters, lora_state_dict, merge_lora, unfreeze_patterns)


class Toy(nn.Module):
    """一个最小的"编码器 + 头"结构，用来复现注入/冻结/合并语义。"""

    def __init__(self, d=16):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(d, d), nn.ReLU(), nn.Linear(d, d))
        self.frontier_pred = nn.Linear(d, 1)

    def forward(self, x):
        return self.frontier_pred(self.encoder(x))


def test_init_is_bit_identical():
    torch.manual_seed(0)
    base = Toy()
    x = torch.randn(5, 16)
    with torch.no_grad():
        y0 = base(x).clone()
    m = copy.deepcopy(base)
    n = apply_lora(m, targets=('encoder',), rank=4, alpha=8.0, verbose=False)
    assert n == 2, n
    with torch.no_grad():
        y1 = m(x)
    assert torch.equal(y0, y1), (y0 - y1).abs().max()
    print("PASS ① 注入后初始前向与基座逐位相同（ΔW 零初始化）")


def test_base_frozen_adapters_trainable():
    torch.manual_seed(0)
    m = Toy()
    apply_lora(m, targets=('encoder',), rank=4, alpha=8.0, verbose=False)
    frozen = [n for n, p in m.named_parameters() if not p.requires_grad]
    trainable = [n for n, p in m.named_parameters() if p.requires_grad]
    # 除 LoRA 适配器外，一切参数都冻结（含被包装的 base 与未包装的头部）
    assert all('lora_' in n for n in trainable), ("不应有非适配器参数可训练", trainable[:5])
    assert any(n.endswith('base.weight') for n in frozen), ("base 权重应冻结", frozen[:5])
    assert 'frontier_pred.weight' in frozen, "未包装的头部也应冻结（要用再显式解冻）"
    tr, tot = count_trainable(m)
    assert 0 < tr < tot
    print("PASS ② 基座与未包装头部全部冻结、仅适配器可训练（可训练 %d / 总 %d = %.1f%%）"
          % (tr, tot, 100.0 * tr / tot))

    # 解冻头部（"适配器 + 头"一起训的常见配方）
    added = unfreeze_patterns(m, ['frontier_pred'])
    assert added > 0
    tr2, _ = count_trainable(m)
    assert tr2 == tr + added
    print("PASS ②b 可再解冻指定头部（新增可训练参数 %d）" % added)


def test_merge_matches_lora_output():
    torch.manual_seed(1)
    m = Toy()
    apply_lora(m, targets=('encoder',), rank=4, alpha=8.0, verbose=False)
    # 随机扰动适配器，使其偏离"零增量"
    with torch.no_grad():
        for p in lora_parameters(m):
            p.add_(torch.randn_like(p) * 0.1)
    x = torch.randn(5, 16)
    with torch.no_grad():
        y_lora = m(x).clone()
    n = merge_lora(m)
    assert n == 2, n
    assert not any(isinstance(mod, LoRALinear) for mod in m.modules()), "合并后仍残留 LoRALinear"
    with torch.no_grad():
        y_merged = m(x)
    diff = (y_lora - y_merged).abs().max().item()
    assert diff < 1e-5, diff
    print("PASS ③ 合并后输出与 LoRA 模型一致（最大差 %.2e），结构还原为普通 nn.Linear" % diff)


def test_adapter_state_dict_roundtrip():
    torch.manual_seed(2)
    m = Toy()
    apply_lora(m, targets=('encoder',), rank=4, alpha=8.0, verbose=False)
    with torch.no_grad():
        for p in lora_parameters(m):
            p.add_(torch.randn_like(p) * 0.05)
    sd = lora_state_dict(m)
    assert all('lora_' in k for k in sd), list(sd)[:3]
    x = torch.randn(3, 16)
    with torch.no_grad():
        y1 = m(x).clone()
    # 清空适配器 → 重新灌回 → 输出应恢复
    with torch.no_grad():
        for p in lora_parameters(m):
            p.zero_()
    n = load_lora_state_dict(m, sd)
    assert n == len(sd)
    with torch.no_grad():
        y2 = m(x)
    assert torch.allclose(y1, y2, atol=0, rtol=0), (y1 - y2).abs().max()
    print("PASS ④ 适配器 state_dict 存取往返无误差（%d 个张量）" % n)


def test_targets_restrict_injection():
    torch.manual_seed(3)
    m = Toy()
    n = apply_lora(m, targets=('frontier_pred',), rank=2, alpha=4.0, verbose=False)
    assert n == 1, n
    assert isinstance(m.encoder[0], nn.Linear), "未匹配的模块不应被包装"
    assert isinstance(m.frontier_pred, LoRALinear)
    print("PASS ⑤ targets 前缀限定生效（只包装匹配层）")


def test_no_match_raises():
    m = Toy()
    try:
        apply_lora(m, targets=('does_not_exist',), rank=2, alpha=4.0, verbose=False)
    except ValueError as e:
        print("PASS ⑥ 无匹配前缀时显式报错：%s" % str(e)[:48])
        return
    raise AssertionError("应当报错")


if __name__ == "__main__":
    test_init_is_bit_identical()
    test_base_frozen_adapters_trainable()
    test_merge_matches_lora_output()
    test_adapter_state_dict_roundtrip()
    test_targets_restrict_injection()
    test_no_match_raises()
    print("ALL PASS (6/6)")
