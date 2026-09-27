# -*- coding: utf-8 -*-
"""离线校验架构实验配置: 不训练, 只静态检查"以为改了、其实没改"的经典陷阱。

检查项
------
对每个传入的 yaml:
  1. 能否按 config.model 正确构建 MaskFillModelVN(含 frontier 扩容键);
  2. frontier_pred 实际参数量是多少(验证扩容真的生效);
  3. 从 init_checkpoint 加载是否成功, 且 strict=false 时**只有 frontier_pred** 有差异;
  4. get_optimizer 生成的参数组是否符合预期(每组的 lr 与参数量),
     尤其确认 encoder 真的落到了小 LR 组、frontier_pred 落到了大 LR 组;
  5. requires_grad=False 的参数总量(确认 freeze_* 行为)。

用法
----
    python src/scripts/verify_arch_configs.py configs/train_gpcr_arch_v7.yml \
                                               configs/train_gpcr_arch_v8.yml
"""
import os
import sys
import logging

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, CONFIGS, src_on_path  # noqa: E402

src_on_path()
import torch  # noqa: E402
from models.maskfill import MaskFillModelVN  # noqa: E402
from utils.misc import load_config  # noqa: E402
from utils.train import get_optimizer  # noqa: E402
from utils.transforms import (FeaturizeProteinAtom, FeaturizeLigandAtom,  # noqa: E402
                              EdgeSample, ContrastiveSample)

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger('verify_arch')

FAILS = []


def chk(name, ok, detail=""):
    print("  %s %s%s" % ("[OK]  " if ok else "[FAIL]", name, ("  -> " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)
    return ok


def build(cfg):
    """完全按 train.py 的方式取维度 —— 绝不硬编码。

    历史教训: 本脚本初版把 protein/ligand 特征维度写死成 32/43, 而实际是 27/13,
    于是每个配置都报"加载失败", 差点被误判成配置本身有问题。
    """
    pf = FeaturizeProteinAtom()
    lf = FeaturizeLigandAtom()
    es = EdgeSample(cfg.train.transform.edgesampler)
    ct = cfg.train.transform.contrastive
    cs = ContrastiveSample(ct.num_real, ct.num_fake, ct.pos_real_std,
                           ct.pos_fake_std, cfg.model.field.knn)
    print("  真实维度: protein=%d ligand=%d num_elements=%d num_bond_types=%d"
          % (pf.feature_dim, lf.feature_dim, cs.num_elements, es.num_bond_types))
    return MaskFillModelVN(cfg.model, num_classes=cs.num_elements,
                           num_bond_types=es.num_bond_types,
                           protein_atom_feature_dim=pf.feature_dim,
                           ligand_atom_feature_dim=lf.feature_dim)


def main(paths):
    for p in paths:
        p = p if os.path.isabs(p) else os.path.join(ROOT, p)
        print("=" * 74)
        print(os.path.relpath(p, ROOT))
        print("=" * 74)
        cfg = load_config(p)

        model = build(cfg)

        # --- 1/2) frontier 容量 ---
        fr = cfg.model.get('frontier', None) or {}
        fr_params = sum(x.numel() for n, x in model.named_parameters() if n.startswith('frontier_pred'))
        want_sca = int(fr.get('hidden_sca', 128))
        want_layers = int(fr.get('layers', 1))
        n_gvp = sum(1 for n, m in model.frontier_pred.net.named_children())
        print("  frontier 配置: hidden_sca=%s layers=%s  -> 实际子模块 %d 个, %.4f M 参数"
              % (want_sca, want_layers, n_gvp, fr_params / 1e6))
        chk("frontier_pred 子模块数 = layers+1", n_gvp == want_layers + 1,
            "得到 %d, 期望 %d" % (n_gvp, want_layers + 1))
        total = sum(x.numel() for x in model.parameters())
        print("  总参数 %.3f M; 各模块占比:" % (total / 1e6))
        import collections
        grp = collections.OrderedDict()
        for n, x in model.named_parameters():
            k = n.split('.')[0]
            grp[k] = grp.get(k, 0) + x.numel()
        for k, v in sorted(grp.items(), key=lambda kv: -kv[1]):
            print("      %-20s %8.4f M  %5.1f%%" % (k, v / 1e6, 100 * v / total))

        # --- 3) 加载 init_checkpoint ---
        ck = cfg.train.get('init_checkpoint', None)
        strict = cfg.train.get('init_strict', True)
        if ck:
            ckp = ck if os.path.isabs(ck) else os.path.join(ROOT, ck)
            if not chk("init_checkpoint 存在", os.path.exists(ckp), ckp):
                continue
            sd = torch.load(ckp, map_location='cpu', weights_only=False)['model']
            if strict:
                try:
                    model.load_state_dict(sd)
                    chk("strict 加载成功", True)
                except Exception as e:
                    chk("strict 加载成功", False, str(e)[:200])
            else:
                # 与 train.py 的 init_strict=false 保持同一套逻辑:
                # PyTorch strict=False 不容忍形状不匹配, 必须先显式剔除。
                sd = dict(sd)
                cur = model.state_dict()
                dropped = [k for k, v in sd.items()
                           if k in cur and tuple(v.shape) != tuple(cur[k].shape)]
                bad_d = [k for k in dropped if 'frontier_pred' not in k]
                chk("strict=false 剔除的形状不符张量仅限 frontier_pred", not bad_d,
                    str(bad_d[:6]))
                for k in dropped:
                    sd.pop(k)
                missing, unexpected = model.load_state_dict(sd, strict=False)
                bad = [k for k in list(missing) + list(unexpected) if 'frontier_pred' not in k]
                chk("strict=false 时差异仅限 frontier_pred", not bad,
                    "missing=%s unexpected=%s" % (list(missing)[:6], list(unexpected)[:6]))
                print("      剔除形状不符 %d 个, 未加载 %d 个, 忽略旧张量 %d 个"
                      % (len(dropped), len(missing), len(unexpected)))
        else:
            print("  (无 init_checkpoint, 跳过加载检查)")

        # --- 4/5) 冻结与参数组 ---
        if cfg.train.get('freeze_encoder', False):
            for x in model.encoder.parameters():
                x.requires_grad = False
        for pat in (cfg.train.get('freeze_patterns', []) or []):
            for n, x in model.named_parameters():
                if pat in n:
                    x.requires_grad = False
        tr = sum(x.numel() for x in model.parameters() if x.requires_grad)
        tot = sum(x.numel() for x in model.parameters())
        print("  可训练 %.3f M / %.3f M (%.1f%%)" % (tr / 1e6, tot / 1e6, 100 * tr / tot))
        chk("存在可训练参数", tr > 0)

        opt = get_optimizer(cfg.train.optimizer, model)
        print("  优化器参数组 %d 个:" % len(opt.param_groups))
        for i, g in enumerate(opt.param_groups):
            n = sum(x.numel() for x in g['params'])
            names = set()
            for x in g['params']:
                for pn, pm in model.named_parameters():
                    if pm is x:
                        names.add(pn.split('.')[0])
                        break
            print("      [%d] lr=%.2e  %.4f M  模块=%s" % (i, g['lr'], n / 1e6, sorted(names)))

        # 分层 LR 必须真的把 encoder 与 frontier_pred 分开
        specs = list(getattr(cfg.train.optimizer, 'param_groups', None) or [])
        if specs:
            pats = {str(s['pattern']): float(s['lr']) for s in specs}
            for pat, lr in pats.items():
                hit = None
                for g in opt.param_groups:
                    if abs(g['lr'] - lr) < 1e-12 and any(
                            pat in pn for x in g['params'] for pn, pm in model.named_parameters()
                            if pm is x):
                        hit = g
                        break
                chk("pattern '%s' 有独立参数组(lr=%.2e)" % (pat, lr), hit is not None)
        enc_lr = None
        for g in opt.param_groups:
            if any(pn.startswith('encoder') for x in g['params']
                   for pn, pm in model.named_parameters() if pm is x):
                enc_lr = g['lr']
        if cfg.train.get('freeze_encoder', False):
            chk("freeze_encoder=true 时 encoder 无可训练参数", enc_lr is None,
                "却有 lr=%s 的组" % enc_lr)
        else:
            chk("解冻后 encoder 有独立/默认学习率", enc_lr is not None, "lr=%s" % enc_lr)
        print("")

    print("=" * 74)
    if FAILS:
        print("结果: FAIL (%d 项): %s" % (len(FAILS), FAILS))
        return 1
    print("结果: PASS —— 配置行为与预期一致")
    return 0


if __name__ == '__main__':
    args = sys.argv[1:] or [os.path.join(CONFIGS, 'train_gpcr_arch_v7.yml'),
                            os.path.join(CONFIGS, 'train_gpcr_arch_v8.yml')]
    raise SystemExit(main(args))
