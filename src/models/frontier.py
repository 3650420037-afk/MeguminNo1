import torch
from torch.nn import Module, Sequential
from torch.nn import functional as F

from .invariant import GVPerceptronVN, GVLinear


class FrontierLayerVN(Module):
    """预测每个配体原子"是否仍是边界原子(frontier)"。

    它同时决定**长在哪里**与**何时停止**: 没有任何 frontier 时生成即终止
    (见 maskfill.py: sample_focal -> has_frontier=False)。官方实现把它写死成
    hidden 128/32, 只占总参数 1.4 %(0.053 M / 3.71 M), 容量与责任严重不匹配 ——
    实测其漂移会让分子永不终止(采样日志 `Finished k | Failed 0`, 完成数 63 -> 23/7)。

    因此这里把容量开放为可配置:
        frontier:
          hidden_sca: 256     # 默认 128
          hidden_vec: 64      # 默认 32
          layers: 2           # 默认 1

    **向后兼容**: n_layers=1 且 hidden 取 128/32 时, self.net 的子模块序列与参数名
    (`net.0.*`, `net.1.*`) 与官方实现完全一致, 可直接加载官方权重。
    一旦放大容量, 官方权重的 frontier_pred 张量形状不再匹配, 必须配合
    `train.init_strict: false` 加载(仅 frontier_pred 允许随机初始化, 其余仍是硬约束)。
    """

    def __init__(self, in_sca, in_vec, hidden_dim_sca, hidden_dim_vec, n_layers=1):
        super().__init__()
        mods = []
        cs, cv = in_sca, in_vec
        for _ in range(max(1, int(n_layers))):
            mods.append(GVPerceptronVN(cs, cv, hidden_dim_sca, hidden_dim_vec))
            cs, cv = hidden_dim_sca, hidden_dim_vec
        mods.append(GVLinear(hidden_dim_sca, hidden_dim_vec, 1, 1))
        self.net = Sequential(*mods)

    def forward(self, h_att, idx_ligans):
        h_att_ligand = [h_att[0][idx_ligans], h_att[1][idx_ligans]]
        pred = self.net(h_att_ligand)
        pred = pred[0]
        return pred
