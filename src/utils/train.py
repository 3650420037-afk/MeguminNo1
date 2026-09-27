import copy
import warnings
import numpy as np
import torch
import torch.nn as nn
from torch_geometric.data import Data, Batch

def repeat_data(data: Data, num_repeat) -> Batch:
    datas = [copy.deepcopy(data) for i in range(num_repeat)]
    return Batch.from_data_list(datas)


def repeat_batch(batch: Batch, num_repeat) -> Batch:
    datas = batch.to_data_list()
    new_data = []
    for i in range(num_repeat):
        new_data += copy.deepcopy(datas)
    return Batch.from_data_list(new_data)


def inf_iterator(iterable):
    iterator = iterable.__iter__()
    while True:
        try:
            yield iterator.__next__()
        except StopIteration:
            iterator = iterable.__iter__()


def get_optimizer(cfg, model):
    """构造优化器, 支持分层学习率。

    为什么需要分层 LR
    -----------------
    微调时若整体解冻, 编码器(2.984 M, 占 80.4 %)会用与头部相同的步长更新, 极易把
    预训练得到的表征冲掉; 而若整体冻结, 模型表征能力根本不变(此前 5 次微调只训了
    0.728 M, 即只做了"头适配")。分层 LR 是这两者之间的正确折中: 编码器用小步长
    缓慢适配, 头部(尤其需要从零学的 frontier_pred)用大步长快速跟上。

    配置写法(yaml)
    --------------
        train:
          optimizer:
            type: adam
            lr: 2e-4                 # 未匹配任何 pattern 的参数用这个
            param_groups:
              - {pattern: encoder,          lr: 1e-5}   # 先匹配先归属
              - {pattern: frontier_pred,    lr: 5e-4}   # 扩容后需从零学, 给更大步长

    说明
    ----
    - pattern 是对参数名的**子串**匹配, 按列表顺序先匹配先归属, 不会重复分组;
    - 只收集 requires_grad=True 的参数, 与 freeze_encoder/freeze_patterns 正交;
    - 未写 param_groups 时行为与旧版完全一致(单组, lr=cfg.lr), 不影响既有实验。
    """
    if cfg.type != 'adam':
        raise NotImplementedError('Optimizer not supported: %s' % cfg.type)
    kw = dict(weight_decay=cfg.weight_decay, betas=(cfg.beta1, cfg.beta2, ))
    specs = list(getattr(cfg, 'param_groups', None) or [])
    if not specs:
        return torch.optim.Adam(model.parameters(), lr=cfg.lr, **kw)

    groups, owner, sizes = [], {}, {}
    for si, sp in enumerate(specs):
        d = sp if isinstance(sp, dict) else dict(sp)
        pat, lr = str(d['pattern']), float(d['lr'])
        params = []
        for name, p in model.named_parameters():
            if p.requires_grad and name not in owner and pat in name:
                owner[name], _ = si, params.append(p)
        if params:
            groups.append({'params': params, 'lr': lr})
            sizes[pat] = (len(params), sum(p.numel() for p in params), lr)
    rest = [p for n, p in model.named_parameters() if p.requires_grad and n not in owner]
    if rest:
        groups.append({'params': rest, 'lr': cfg.lr})
    if not groups:
        raise RuntimeError('param_groups 未匹配到任何 requires_grad=True 的参数, '
                           '请检查 pattern 是否写错(例如模块名拼写)')
    return torch.optim.Adam(groups, lr=cfg.lr, **kw)


def get_scheduler(cfg, optimizer):
    if cfg.type == 'plateau':
        return torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            factor=cfg.factor,
            patience=cfg.patience,
            min_lr=cfg.min_lr
        )
    else:
        raise NotImplementedError('Scheduler not supported: %s' % cfg.type)

