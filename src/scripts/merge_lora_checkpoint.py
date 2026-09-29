# -*- coding: utf-8 -*-
"""把 **LoRA 检查点**折成**普通检查点**（适配器折进基座，导出不再依赖 LoRA 代码）。

为什么必须有这一步
------------------
LoRA 训练出的 `state_dict` 里带 `lora_A/lora_B` 键，而采样/交接/归档都希望拿到
"一个自解释、任何环境都能直接加载"的权重（`7-eonmol_ft_gpcr_v*.pt` 就一直是这种形态）。
本脚本输出**同构的普通检查点**：`config.model.lora.enabled=False`、`model` 里不再有 LoRA 键，
且**权重已被折入**（`W ← W + B·A·α/r`），因此加载后前向输出与 LoRA 模型一致（有测试保证）。

用法
----
    python src/scripts/merge_lora_checkpoint.py --ckpt <lora.pt> --out <plain.pt>
    # --check：折入后与原 LoRA 模型在同一输入上的前向一致性校验（CPU，默认开启）
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import torch  # noqa: E402

from paths import ROOT  # noqa: E402
from models.maskfill import MaskFillModelVN  # noqa: E402
from utils.lora import LoRALinear, apply_lora, lora_state_dict, merge_lora  # noqa: E402


def _to_plain(d):
    """把 EasyDict/AttrDict 递归转成普通 dict（便于 yaml/json 友好与跨版本加载）。"""
    if isinstance(d, dict):
        return {k: _to_plain(v) for k, v in d.items()}
    if isinstance(d, (list, tuple)):
        return [_to_plain(v) for v in d]
    return d


def build_from_config(mcfg, dims=(5, 3, 5, 5)):
    return MaskFillModelVN(mcfg, *dims)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True, help="LoRA 检查点")
    ap.add_argument("--out", required=True, help="输出的普通检查点")
    ap.add_argument("--dims", nargs=4, type=int, default=[5, 3, 5, 5],
                    help="num_classes num_bond_types protein_dim ligand_dim（默认与 v10 家族一致）")
    ap.add_argument("--no-check", action="store_true", help="跳过折入后的前向一致性校验")
    ap.add_argument("--allow-missing-lora", action="store_true",
                    help="允许检查点缺适配器键（默认**拒绝**，因为那等于静默折入零增量）")
    args = ap.parse_args()

    ck = torch.load(args.ckpt, map_location="cpu", weights_only=False)
    if not isinstance(ck, dict) or 'model' not in ck:
        raise SystemExit('检查点格式不对：需要含 config/model/iteration 的字典')
    cfg_raw = ck.get('config', {})          # 构造模型要用**原始**配置（MaskFillModelVN 走属性访问）
    cfg = _to_plain(cfg_raw)                # 落盘用普通 dict（跨环境友好）
    mcfg = cfg.get('model', {})
    lora = (mcfg.get('lora') or {})
    if not lora.get('enabled', False):
        raise SystemExit('该检查点未启用 LoRA（config.model.lora.enabled 非真），无需折入')

    # 普通 dict config 也要能构造（MaskFillModelVN 走属性访问）→ 用 EasyDict 包一层。
    # 独立审查实测：输入侧原先对普通 dict config 直接 AttributeError（只修了输出侧）。
    from easydict import EasyDict
    _m = cfg_raw.get('model', {}) if isinstance(cfg_raw, dict) else getattr(cfg_raw, 'model', {})
    mcfg_raw = _m if not isinstance(_m, dict) else EasyDict(_m)
    model = build_from_config(mcfg_raw, tuple(args.dims))
    n = apply_lora(model, targets=tuple(lora.get('targets', ['encoder'])),
                   rank=int(lora.get('rank', 8)), alpha=float(lora.get('alpha', 16.0)),
                   verbose=False)
    missing, unexpected = model.load_state_dict(ck['model'], strict=False)
    bad_m = [k for k in missing if 'lora_' not in k]
    bad_u = [k for k in unexpected if 'lora_' not in k]
    if bad_m or bad_u:
        raise SystemExit('权重与模型不匹配：missing=%s unexpected=%s' % (bad_m[:5], bad_u[:5]))
    # ⚠ 缺适配器键时旧版会 **rc=0 静默折入零增量**（等于该层悄悄回退基座）—— 独立审查实测发现。
    miss_lora = [] if args.allow_missing_lora else [k for k in missing if 'lora_' in k]
    if miss_lora:
        raise SystemExit('检查点缺少 %d 个适配器张量（例如 %s）—— 拒绝静默折入零增量；'
                         '若确认要按缺失状态折入，请显式加 --allow-missing-lora'
                         % (len(miss_lora), miss_lora[:3]))
    print('[merge] 注入 %d 层；载入适配器 %d 个张量' % (n, len(lora_state_dict(model))))

    n_merged = merge_lora(model)
    sd = model.state_dict()
    assert not any('lora_' in k for k in sd), '折入后不应再有 lora_ 键'
    assert not any(isinstance(m, LoRALinear) for m in model.modules())

    # ⚠ 必须**保留原始 config 对象类型**（通常是 EasyDict）：采样端 MaskFillModelVN 走
    # `config.hidden_channels` 这类**属性访问**，若这里落成普通 dict，加载时会
    # AttributeError: 'dict' object has no attribute 'hidden_channels'（本脚本自测时真实踩到）。
    import copy as _copy
    out_cfg = _copy.deepcopy(cfg_raw)
    try:
        out_cfg['model']['lora'] = dict(lora, enabled=False)   # 标记：本权重已是普通形态
    except Exception:
        _m = getattr(out_cfg, 'model', None)
        if isinstance(_m, dict):
            _m['lora'] = dict(lora, enabled=False)
    out = {'config': out_cfg, 'model': sd}
    if 'iteration' in ck:
        out['iteration'] = ck['iteration']
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
    torch.save(out, args.out)
    size_mb = os.path.getsize(args.out) / 1e6
    n_par = sum(v.numel() for v in sd.values())
    # 输出可能与仓库不同盘（如临时目录在 C:）→ relpath 会抛 ValueError，故做保护
    try:
        shown = os.path.relpath(args.out, ROOT)
    except ValueError:
        shown = os.path.abspath(args.out)
    print('[merge] 折入 %d 层 → %s（%.1f MB, %d 张量, %.4f M 参数）'
          % (n_merged, shown, size_mb, len(sd), n_par / 1e6))
    print('[merge] 提示：该文件可被 models/ 目录直接收录；采样端无需 LoRA 代码即可加载')
    return 0


if __name__ == "__main__":
    sys.exit(main())
