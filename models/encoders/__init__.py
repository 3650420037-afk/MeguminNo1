from .cftfm import CFTransformerEncoderVN


def get_encoder_vn(config):
    """构建等变编码器。

    NOTE (架构澄清): 本编码器是**等变图消息传递网络 (MPNN)**, 不是 Transformer ——
    CFTransformerEncoderVN 内部为 num_interactions 层 AttentionInteractionBlockVN
    (距离/向量展开 + MessageModule + LayerNorm + 残差), **没有多头注意力**。
    下方的 key_channels 与 num_heads 属原版历史配置项, 传入后**未被使用**
    (cftfm.py 内亦标注 "not use"); 保留仅为兼容既有 configs, 调整它们不影响模型。
    项目中真正的注意力机制位于 field 的键预测头 (models/fields/classifier.py 的
    AttentionEdges, 含 tri-edge 键型偏置)。
    """
    if config.name == 'cftfm':
        return CFTransformerEncoderVN(
            hidden_channels = [config.hidden_channels, config.hidden_channels_vec],
            edge_channels = config.edge_channels,
            key_channels = config.key_channels,  # not use
            num_heads = config.num_heads,  # not use
            num_interactions = config.num_interactions,
            k = config.knn,
            cutoff = config.cutoff,
        )
    else:
        raise NotImplementedError('Unknown encoder: %s' % config.name)
