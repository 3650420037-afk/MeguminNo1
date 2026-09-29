# -*- coding: utf-8 -*-
"""权重加载兼容性检查（CPU）：每个检查点都能按其自带 config 严格加载。

为什么要它：采样端（`src/sample_for_pdb.py`）刚加了 LoRA 分支，加载路径被改过；
而 `models/` 下几个权重是交付物，一旦加载出问题（键缺失/形状不符/静默随机初始化）后果严重。
本脚本把每个权重走一遍**与采样端相同的构建+加载路径**，并报告：
参数量、张量数、strict 加载结果、以及 config 里的关键字段（frontier 容量 / lora 开关）。

用法：python src/scripts/check_weight_compat.py [--dirs models other/weight_backups/*]
"""
import argparse
import glob
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

def _dims():
    """与采样端完全一致的维度取自 featurizer/edge sampler，不要手写猜（自测时猜错过一次）。"""
    from utils.transforms import (ContrastiveSample, EdgeSample, FeaturizeLigandAtom,
                                  FeaturizeProteinAtom)
    pf, lf = FeaturizeProteinAtom(), FeaturizeLigandAtom()
    import yaml
    with open(os.path.join(ROOT, 'configs', 'train_gpcr_mass_v10term.yml'), encoding='utf-8') as fh:
        cfg = yaml.safe_load(fh)
    from easydict import EasyDict
    c = EasyDict(cfg['train']['transform']['contrastive'])
    es = EdgeSample(EasyDict(cfg['train']['transform']['edgesampler']))   # 需要属性访问
    cs = ContrastiveSample(c.num_real, c.num_fake, c.pos_real_std, c.pos_fake_std,
                           cfg['model']['field']['knn'])
    return (cs.num_elements, es.num_bond_types, pf.feature_dim, lf.feature_dim)


DIMS = None      # 首次使用时惰性计算


def check(path):
    ck = torch.load(path, map_location="cpu", weights_only=False)
    if not isinstance(ck, dict) or 'model' not in ck or 'config' not in ck:
        return dict(path=path, ok=False, note='格式不符（缺 config/model）')
    cfg = ck['config']
    mcfg = cfg.get('model', {}) if isinstance(cfg, dict) else getattr(cfg, 'model', {})
    lora = (mcfg.get('lora', None) if isinstance(mcfg, dict)
            else getattr(mcfg, 'lora', None)) or {}
    global DIMS
    if DIMS is None:
        DIMS = _dims()
    model = MaskFillModelVN(mcfg, *DIMS)
    n_lora = 0
    if lora.get('enabled', False):
        from utils.lora import apply_lora
        n_lora = apply_lora(model, targets=tuple(lora.get('targets', ['encoder'])),
                            rank=int(lora.get('rank', 8)), alpha=float(lora.get('alpha', 16.0)),
                            verbose=False)
    sd = ck['model']
    try:
        model.load_state_dict(sd)                      # strict=True：键与形状必须完全一致
        strict = 'OK(strict)'
    except Exception as e:
        try:
            missing, unexpected = model.load_state_dict(sd, strict=False)
            strict = 'WARN: missing=%d unexpected=%d (%s)' % (
                len(missing), len(unexpected), str(e)[:60])
        except Exception as e2:
            return dict(path=path, ok=False, note='加载失败: %s' % str(e2)[:80])
    fr = mcfg.get('frontier', {}) if isinstance(mcfg, dict) else getattr(mcfg, 'frontier', {})
    fr = fr if isinstance(fr, dict) else {}
    n_par = sum(v.numel() for v in sd.values())
    return dict(path=os.path.relpath(path, ROOT), ok=('OK' in strict), note=strict,
                tensors=len(sd), params_M=round(n_par / 1e6, 4),
                frontier='%s/%s/%s' % (fr.get('hidden_sca', 128), fr.get('hidden_vec', 32),
                                       fr.get('layers', 1)),
                lora=('on(rank=%s, %d 层)' % (lora.get('rank'), n_lora)) if n_lora else 'off',
                iteration=ck.get('iteration'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="+", default=["models"])
    ap.add_argument("--recursive", action="store_true", help="递归含子目录（备份目录用）")
    args = ap.parse_args()

    files = []
    for d in args.dirs:
        pat = os.path.join(ROOT, d, "**", "*.pt") if args.recursive else os.path.join(ROOT, d, "*.pt")
        files += sorted(glob.glob(pat, recursive=args.recursive))
    files = [f for f in files if os.path.getsize(f) < 200 * 1e6]     # 跳过超大文件
    if not files:
        print("没有找到 .pt 权重")
        return 1

    bad = 0
    print("%-44s %6s %9s %-9s %-8s %s" % ("权重", "张量", "参数(M)", "frontier", "LoRA", "加载"))
    for f in files:
        r = check(f)
        if not r.get('ok'):
            bad += 1
            print("%-44s %s" % (r['path'][:44], r.get('note')))
            continue
        print("%-44s %6d %9.4f %-9s %-8s %s | iter=%s"
              % (r['path'][:44], r['tensors'], r['params_M'], r['frontier'], r['lora'],
                 r['note'], r['iteration']))
    print("\n结果：%d 个权重，%d 个有问题" % (len(files), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
