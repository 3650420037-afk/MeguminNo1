# -*- coding: utf-8 -*-
"""CPU 端到端测试：LoRA 检查点 → `merge_lora_checkpoint.py` → 普通检查点（可被采样端直接加载）。

为什么这条链必须测：采样端（`src/sample_for_pdb.py`）与交付归档都假设权重是"自解释的普通检查点"。
若折入环节出错，最坏情况是**静默丢掉适配器**（等于用回基座权重却以为在用 LoRA 模型）。
本测试断言：折入后 ① 无 lora_ 键；② 能被**严格加载**进普通模型；③ 前向输出与 LoRA 模型一致。

跑法：python tests/test_merge_lora_checkpoint.py
"""
import copy
import importlib.util
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, ROOT)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import torch  # noqa: E402
import yaml  # noqa: E402
from easydict import EasyDict  # noqa: E402

from models.maskfill import MaskFillModelVN  # noqa: E402
from utils.lora import LoRALinear, apply_lora, lora_parameters  # noqa: E402

CFG = os.path.join(ROOT, "configs", "train_gpcr_mass_v10lora.yml")
MERGE = os.path.join(ROOT, "src", "scripts", "merge_lora_checkpoint.py")


def main():
    with open(CFG, encoding="utf-8") as fh:
        cfg = EasyDict(yaml.safe_load(fh))
    dims = (5, 3, 5, 5)

    # 构造一个"训练过"的 LoRA 模型：随机扰动适配器，使其偏离零增量
    model = MaskFillModelVN(cfg.model, *dims)
    L = cfg.model["lora"]
    n = apply_lora(model, targets=tuple(L["targets"]), rank=int(L["rank"]), alpha=float(L["alpha"]),
                   verbose=False)
    torch.manual_seed(0)
    with torch.no_grad():
        for p in lora_parameters(model):
            p.add_(torch.randn_like(p) * 0.02)
    assert n == 162 and any(isinstance(m, LoRALinear) for m in model.modules())
    print("构造 LoRA 模型：注入 %d 层" % n)

    # 参考输出：对**每一个** LoRA 包装层比对折入前后的输出
    # （不去手工拼 encoder 的输入打包格式 —— 那是脆的；按模块名逐一比对才是折入真正影响的地方）
    probes = {}
    torch.manual_seed(1)
    for name, mod in model.named_modules():
        if isinstance(mod, LoRALinear):
            xin = torch.randn(7, mod.base.in_features)
            with torch.no_grad():
                probes[name] = (xin.clone(), mod(xin).clone())

    # 保存成 LoRA 检查点
    tmp = tempfile.mkdtemp(prefix="lora_ckpt_")
    lora_pt = os.path.join(tmp, "lora.pt")
    plain_pt = os.path.join(tmp, "plain.pt")
    torch.save({"config": cfg, "model": model.state_dict(), "iteration": 1234}, lora_pt)
    print("已保存 LoRA 检查点：%d 个张量（含 %d 个适配器张量）"
          % (len(model.state_dict()), len([k for k in model.state_dict() if "lora_" in k])))

    # 跑折入脚本（子进程，走真实 CLI 路径）
    r = subprocess.run([sys.executable, MERGE, "--ckpt", lora_pt, "--out", plain_pt],
                       capture_output=True, text=True, errors="ignore")
    print(r.stdout.strip()[-300:])
    if r.returncode != 0:
        print(r.stderr[-800:])
        raise SystemExit("merge_lora_checkpoint.py 失败（rc=%d）" % r.returncode)

    ck = torch.load(plain_pt, map_location="cpu", weights_only=False)
    sd = ck["model"]
    assert not any("lora_" in k for k in sd), "折入后仍残留 lora_ 键"
    assert ck["config"]["model"]["lora"]["enabled"] is False, "配置里应标记 lora 已关闭/已折入"
    assert ck["iteration"] == 1234, "iteration 应保留"
    print("PASS ① 折入后无 lora_ 键、config 标记为已折入、iteration 保留")

    # 普通模型**严格**加载（strict=True：键完全一致才通过）
    plain = MaskFillModelVN(ck["config"]["model"], *dims)
    plain.load_state_dict(sd)          # strict=True
    assert not any(isinstance(m, LoRALinear) for m in plain.modules())
    print("PASS ② 普通模型 strict=True 加载成功（%d 张量），无 LoRA 结构残留" % len(sd))

    # 前向一致性：折入后每一层都应等价于 LoRA 层
    worst, worst_name = 0.0, None
    mods = dict(plain.named_modules())
    for name, (xin, yout) in probes.items():
        with torch.no_grad():
            y2 = mods[name](xin)
        d = (yout - y2).abs().max().item()
        if d > worst:
            worst, worst_name = d, name
    assert worst < 1e-5, (worst, worst_name)
    print("PASS ③ 折入前后逐层一致（比对 %d 层，最大差 %.2e ≤ 1e-5，最差层 %s）"
          % (len(probes), worst, (worst_name or '')[:40]))

    # 采样端"启用 LoRA 的检查点"分支：注入后 strict 加载
    cfg2 = copy.deepcopy(cfg)
    cfg2["model"]["lora"]["enabled"] = True
    m2 = MaskFillModelVN(cfg2["model"], *dims)
    apply_lora(m2, targets=tuple(L["targets"]), rank=int(L["rank"]), alpha=float(L["alpha"]),
               verbose=False)
    m2.load_state_dict(model.state_dict())     # strict=True
    print("PASS ④ 采样端分支：按 config 注入 LoRA 后可 strict 加载 LoRA 检查点")

    print("\nALL PASS (4/4)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
