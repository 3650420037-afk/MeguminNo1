import torch
from torch.nn import Module
from torch.nn import functional as F

from .encoders import get_encoder_vn
from .fields import get_field_vn
from .common import *
from .embedding import AtomEmbedding
from .frontier import FrontierLayerVN
from .position import PositionPredictor
# from .debug import check_true_bonds_len, check_pred_bonds_len
from utils.misc import unique

class MaskFillModelVN(Module):

    def __init__(self, config, num_classes, num_bond_types, protein_atom_feature_dim, ligand_atom_feature_dim):
        super().__init__()
        self.config = config
        self.num_bond_types = num_bond_types
        
        self.emb_dim = [config.hidden_channels, config.hidden_channels_vec]
        self.protein_atom_emb = AtomEmbedding(protein_atom_feature_dim, 1, *self.emb_dim)
        self.ligand_atom_emb = AtomEmbedding(ligand_atom_feature_dim, 1, *self.emb_dim)

        self.encoder = get_encoder_vn(config.encoder)
        in_sca, in_vec = self.encoder.out_sca, self.encoder.out_vec
        self.field = get_field_vn(config.field, num_classes=num_classes, num_bond_types=num_bond_types, 
                                             in_sca=in_sca, in_vec=in_vec)
        # frontier_pred 决定"哪里还能生长 / 何时停止", 是生成能否终止的总开关。
        # 官方把它写死成 hidden 128/32(仅 0.053 M, 占 1.4%), 容量与责任不匹配。
        # 现开放为 config.frontier.{hidden_sca,hidden_vec,layers}; 缺省仍是
        # 128/32/1, 与官方参数名和形状完全一致(向后兼容)。
        _fr = config.get('frontier', None) or {}
        self.frontier_pred = FrontierLayerVN(in_sca=in_sca, in_vec=in_vec,
                                            hidden_dim_sca=int(_fr.get('hidden_sca', 128)),
                                            hidden_dim_vec=int(_fr.get('hidden_vec', 32)),
                                            n_layers=int(_fr.get('layers', 1)))
        # self.protein_frontier_pred = FrontierLayerVN(in_sca=in_sca, in_vec=in_vec,
        #                                                                     hidden_dim_sca=128, hidden_dim_vec=32)
        self.pos_predictor = PositionPredictor(in_sca=in_sca, in_vec=in_vec,
                                                                            num_filters=[config.position.num_filters]*2, n_component=config.position.n_component)

        self.smooth_cross_entropy = SmoothCrossEntropyLoss(reduction='mean', smoothing=0.1)
        self.bceloss_with_logits = nn.BCEWithLogitsLoss()

    def sample_init(self,
        compose_feature,
        compose_pos,
        idx_protein,
        compose_knn_edge_index,
        compose_knn_edge_feature,
        n_samples_pos=-1,
        n_samples_atom=-1,
        frontier_threshold=0,
        ):
        # frontier_threshold 默认 0 (= logit > 0, 等价 sigmoid > 0.5) 与官方行为一致;
        # 它是**采样侧校准参数**: 微调后 frontier 输出分布会偏移, 沿用预训练期的阈值
        # 会让模型过早判定"没有边界原子"从而终止, 生成分子偏小。
        idx_ligand = torch.empty(0).to(idx_protein)  # fake index of ligand
        focal_resutls = self.sample_focal(compose_feature, compose_pos, idx_ligand, idx_protein, compose_knn_edge_index, compose_knn_edge_feature, frontier_threshold=frontier_threshold)
        if focal_resutls[0]:  # has frontiers
            has_frontier, idx_frontier, p_frontier, idx_focal_in_compose, p_focal, h_compose = focal_resutls
            pos_generated, pdf_pos, idx_parent, abs_pos_mu, pos_sigma, pos_pi = self.sample_position(
                h_compose, compose_pos, idx_focal_in_compose, n_samples=n_samples_pos,
            )
            idx_focal_in_compose, p_focal = idx_focal_in_compose[idx_parent], p_focal[idx_parent]
            element_pred, element_prob, has_atom_prob, idx_parent = self.sample_init_element(
                pos_generated, h_compose, compose_pos, n_samples=n_samples_atom
            )
            pos_generated, pdf_pos, idx_focal_in_compose, p_focal = pos_generated[idx_parent], pdf_pos[idx_parent], idx_focal_in_compose[idx_parent], p_focal[idx_parent]
            return (has_frontier, idx_frontier, p_frontier,  # frontier
                        idx_focal_in_compose, p_focal,  # focal
                        pos_generated, pdf_pos, abs_pos_mu, pos_sigma, pos_pi,  # positions
                        element_pred, element_prob, has_atom_prob,  # element
                        )
        else:
            return (False, )

    def sample(self,
        compose_feature,
        compose_pos,
        idx_ligand,
        idx_protein,
        compose_knn_edge_index,
        compose_knn_edge_feature,
        ligand_context_bond_index,
        ligand_context_bond_type,
        n_samples_pos=-1,
        n_samples_atom=-1,
        # n_samples=5,
        frontier_threshold=0,
        ):
        focal_resutls = self.sample_focal(compose_feature, compose_pos, idx_ligand, idx_protein, compose_knn_edge_index, compose_knn_edge_feature, frontier_threshold=frontier_threshold)
        if focal_resutls[0]:  # has frontiers
            has_frontier, idx_frontier, p_frontier, idx_focal_in_compose, p_focal, h_compose = focal_resutls
            pos_generated, pdf_pos, idx_parent, abs_pos_mu, pos_sigma, pos_pi = self.sample_position(
                h_compose, compose_pos, idx_focal_in_compose, n_samples=n_samples_pos
            )
            idx_focal_in_compose, p_focal = idx_focal_in_compose[idx_parent], p_focal[idx_parent]
            element_pred, element_prob, has_atom_prob, idx_parent, bond_index, bond_type, bond_prob = self.sample_element_and_bond(
                pos_generated, h_compose, compose_pos, idx_ligand, ligand_context_bond_index, ligand_context_bond_type, n_samples=n_samples_atom
            )
            pos_generated, pdf_pos, idx_focal_in_compose, p_focal = pos_generated[idx_parent], pdf_pos[idx_parent], idx_focal_in_compose[idx_parent], p_focal[idx_parent]
            return (has_frontier, idx_frontier, p_frontier,  # frontier
                        idx_focal_in_compose, p_focal,  # focal
                        pos_generated, pdf_pos, abs_pos_mu, pos_sigma, pos_pi,  # positions
                        element_pred, element_prob, has_atom_prob,  # element
                        bond_index, bond_type, bond_prob  # bond
                        )
        else:
            return (False, )

    def sample_focal(self,
            compose_feature,
            compose_pos,
            idx_ligand,
            idx_protein,
            compose_knn_edge_index,
            compose_knn_edge_feature,
            n_samples=-1,
            frontier_threshold=0,
        ):
        # # 0: encode 
        h_compose = embed_compose(compose_feature, compose_pos, idx_ligand, idx_protein,
                                      self.ligand_atom_emb, self.protein_atom_emb, self.emb_dim)
        h_compose = self.encoder(
            node_attr = h_compose,
            pos = compose_pos,
            edge_index = compose_knn_edge_index,
            edge_feature = compose_knn_edge_feature,
        )
        # # For the initial atom
        if len(idx_ligand) == 0:
            idx_ligand = idx_protein
        # # 1: predict frontier
        y_frontier_pred = self.frontier_pred(
            h_compose,
            idx_ligand,
        )[:, 0]
        ind_frontier = (y_frontier_pred > frontier_threshold)
        has_frontier = torch.sum(ind_frontier) > 0
        frontier_scale = 1
        if has_frontier:
            # # 2: sample focal from frontiers
            idx_frontier = idx_ligand[ind_frontier]
            p_frontier = torch.sigmoid(y_frontier_pred[ind_frontier])
            if n_samples > 0:  # sample from frontiers
                p_frontier_in_compose = torch.zeros(len(compose_pos), dtype=torch.float32, device=compose_pos.device)
                p_frontier_in_compose_sf = torch.zeros_like(p_frontier_in_compose)
                p_frontier_in_compose_sf[idx_frontier] = F.softmax(p_frontier / frontier_scale, dim=0)
                p_frontier_in_compose[idx_frontier] = p_frontier
                idx_focal_in_compose = p_frontier_in_compose_sf.multinomial(num_samples=n_samples, replacement=True)
                p_focal = p_frontier_in_compose[idx_focal_in_compose]
            else:  # get all frontiers as focal
                idx_focal_in_compose = torch.nonzero(ind_frontier)[:, 0]
                p_focal = p_frontier

            return (has_frontier, idx_frontier, p_frontier,  # frontier
                        idx_focal_in_compose, p_focal,  # focal
                        h_compose)
        else:
            return (has_frontier, h_compose)

    def sample_position(self,
        h_compose,
        compose_pos,
        idx_focal_in_compose,
        n_samples=-1,
        ):
        n_focals = len(idx_focal_in_compose)
        # # 3: get position distributions and sample positions
        relative_pos_mu, abs_pos_mu, pos_sigma, pos_pi  = self.pos_predictor(
            h_compose,
            idx_focal_in_compose,
            compose_pos,
        )
        if n_samples < 0:
            pos_generated = self.pos_predictor.get_maximum(abs_pos_mu, pos_sigma, pos_pi,)  # n_focals, n_per_pos, 3
            n_candidate_samples = pos_generated.size(1)
            pos_generated = torch.reshape(pos_generated, [-1, 3])
            pdf_pos = self.pos_predictor.get_mdn_probability(
                mu=torch.repeat_interleave(abs_pos_mu, repeats=n_candidate_samples, dim=0),
                sigma=torch.repeat_interleave(pos_sigma, repeats=n_candidate_samples, dim=0),
                pi=torch.repeat_interleave(pos_pi, repeats=n_candidate_samples, dim=0),
                pos_target=pos_generated
            )
            idx_parent = torch.repeat_interleave(torch.arange(n_focals), repeats=n_candidate_samples, dim=0).to(compose_pos.device)
        else:
            # 原版遗留: n_samples>=0 (MDN 随机采样) 分支未实现, 旧代码会带着未定义变量
            # 落到下面的 return 抛 UnboundLocalError。此处显式报错以免误导;
            # 如需随机采样, 可接入 PositionPredictor.sample_batch (position.py 已备)。
            raise NotImplementedError(
                'sample_position(n_samples>=0) 未实现; 当前调用方均使用 n_samples=-1 (确定性取 mu)')

        return (pos_generated, pdf_pos, idx_parent, abs_pos_mu, pos_sigma, pos_pi)  # position

    def sample_element_and_bond(self, 
        pos_generated,
        h_compose,
        compose_pos,
        idx_ligand,
        ligand_bond_index,
        ligand_bond_type,
        n_samples
        ):
        # # 4: query positions 
        #NOTE: Only one parent batch (one compose graph) at a time (i.e. batch size = 1)
        n_query = len(pos_generated)
        n_context = len(idx_ligand)
        y_query_pred, edge_pred = self.query_position(
            pos_query = pos_generated,
            h_compose = h_compose,
            compose_pos = compose_pos,
            idx_ligand = idx_ligand,
            ligand_bond_index = ligand_bond_index,
            ligand_bond_type = ligand_bond_type
        )
        if n_samples < 0:
            # raise NotImplementedError('The following is not fixed (and for/ bond)')
            has_atom_prob =  1 - 1 / (1 + torch.exp(y_query_pred).sum(-1))
            y_query_pred = F.softmax(y_query_pred, dim=-1)
            element_pred = y_query_pred.argmax(dim=-1)  # multinomial(1)[:, 0]
            element_prob = y_query_pred[torch.arange(len(y_query_pred)), element_pred]
            idx_parent = torch.arange(n_query)
        else:
            has_atom_prob =  (1 - 1 / (1 + torch.exp(y_query_pred).sum(-1)))
            has_atom_prob = torch.repeat_interleave(has_atom_prob, n_samples, dim=0)  # n_query * n_samples
            y_query_pred = F.softmax(y_query_pred, dim=-1)
            element_pred = y_query_pred.multinomial(n_samples, replacement=True).reshape(-1)  # n_query * n_samples
            idx_parent = torch.repeat_interleave(torch.arange(n_query), n_samples, dim=0).to(compose_pos.device)
            element_prob = y_query_pred[idx_parent, element_pred]
        # # 5: determine bonds
        if n_samples < 0:
            all_edge_type = torch.argmax(edge_pred, dim=-1) # (num_generated, num_ligand_context)
            bond_index = torch.stack(torch.where(
                all_edge_type > 0,
            ), dim=0)
            bond_type = all_edge_type[bond_index[0], bond_index[1]]
            bond_prob = F.softmax(edge_pred, dim=-1)[bond_index[0], bond_index[1], bond_type]
        else:
            edge_pred = F.softmax(edge_pred, dim=-1)  # (num_query, num_context, 4)
            edge_pred_flat = edge_pred.reshape([n_query * n_context, -1])  # (num_query * num_context, 4)
            all_edge_type = edge_pred_flat.multinomial(n_samples, replacement=True)  # (num_query * num_context, n_samples)
            all_edge_type = all_edge_type.reshape([n_query, n_context, n_samples])  # (num_query, num_context, n_samples)
            all_edge_type = all_edge_type.transpose(1, 2)  # (num_query, n_samples, num_context)
            all_edge_type = all_edge_type.reshape([n_query * n_samples, n_context]) # (num_generated * n_samples, num_ligand_context)
            # drop duplicates
            id_element_and_bond = torch.cat([idx_parent.unsqueeze(-1), element_pred.unsqueeze(-1), all_edge_type], dim=1)
            id_element_and_bond, index_unique = unique(id_element_and_bond, dim=0)
            # all_edge_type = all_edge_type[index_unique]
            element_pred, element_prob, has_atom_prob, idx_parent = element_pred[index_unique], element_prob[index_unique], has_atom_prob[index_unique], idx_parent[index_unique]

            # get bond index
            all_edge_type = all_edge_type[index_unique]
            bond_index = torch.stack(torch.where(
                all_edge_type > 0,
            ), dim=0)
            bond_type = all_edge_type[bond_index[0], bond_index[1]]
            bond_prob = edge_pred[idx_parent[bond_index[0]], bond_index[1], bond_type]
        
            
        return (element_pred, element_prob, has_atom_prob, idx_parent, # element
                    bond_index, bond_type, bond_prob  # bond
                    )

    def sample_init_element(self, 
        pos_generated,
        h_compose,
        compose_pos,
        n_samples,
        ):
        # # 4: query positions 
        #NOTE: Only one parent batch (one compose graph) at a time (i.e. batch size = 1)
        n_query = len(pos_generated)
        query_compose_knn_edge_index = knn(x=compose_pos, y=pos_generated, k=self.config.field.knn, num_workers=16)
        y_query_pred, _ = self.field(
            pos_query = pos_generated,
            edge_index_query = [],
            pos_compose = compose_pos,
            node_attr_compose = h_compose,
            edge_index_q_cps_knn = query_compose_knn_edge_index,
        )
        if n_samples < 0:
            # raise NotImplementedError('The following is not fixed')
            has_atom_prob =  1 - 1 / (1 + torch.exp(y_query_pred).sum(-1))
            y_query_pred = F.softmax(y_query_pred, dim=-1)
            element_pred = y_query_pred.argmax(dim=-1)
            element_prob = y_query_pred[torch.arange(len(y_query_pred)), element_pred]
            idx_parent = torch.arange(n_query).to(compose_pos.device)
        else:
            has_atom_prob =  (1 - 1 / (1 + torch.exp(y_query_pred).sum(-1)))
            has_atom_prob = torch.repeat_interleave(has_atom_prob, n_samples, dim=0)  # n_query * n_samples
            y_query_pred = F.softmax(y_query_pred, dim=-1)
            element_pred = y_query_pred.multinomial(n_samples, replacement=True).reshape(-1)  # n_query, n_samples
            idx_parent = torch.repeat_interleave(torch.arange(n_query), n_samples, dim=0).to(compose_pos.device)
            element_prob = y_query_pred[idx_parent, element_pred]
            # drop duplicates
            identifier = torch.stack([idx_parent, element_pred], dim=1)
            identifier, index_unique = unique(identifier, dim=0)

            element_pred, element_prob, has_atom_prob, idx_parent = element_pred[index_unique], element_prob[index_unique], has_atom_prob[index_unique], idx_parent[index_unique]

        return (element_pred, element_prob, has_atom_prob, idx_parent) # element

    def get_loss(self, pos_real, y_real, pos_fake,  # query real positions,
                          index_real_cps_edge_for_atten, tri_edge_index, tri_edge_feat,  # for edge attention
                          edge_index_real, edge_label,  # edges to predict
                          compose_feature, compose_pos, idx_ligand, idx_protein,  # compose (protein and context ligand) atoms
                          y_frontier,  # frontier labels
                          idx_focal,  pos_generate,  # focal and generated positions  #NOTE: idx are in comopse
                          idx_protein_all_mask, y_protein_frontier,  # surface of protein
                          compose_knn_edge_index, compose_knn_edge_feature, real_compose_knn_edge_index,  fake_compose_knn_edge_index  # edges in compose, query-compose
        ):

        # # 边界守卫: 上下文为空时(配体被全部掩码), frontier 与位置预测没有任何可监督的目标,
        # # BCE / 负对数似然作用在空张量上会返回 NaN。train.py 的 loss_components 有限性检查
        # # 会因此跳过整批迭代 —— 表现为"每约 11 个 batch 就有 1 个被静默跳过"(实测 8.6%~9.4%,
        # # 且与数据集内容无关)。
        # # 触发原因在配置: train.transform.mask 的 max_ratio=1.1 且 min_num_unmasked=0 时,
        # # 掩码比例 = clip(U(0,1.1),0,1) 有 1/11 = 9.09% 的概率取到 1.0。
        # # 配置侧已统一改为 min_num_unmasked: 1; 这里再加一层兜底, 返回零损失而非 NaN,
        # # 语义为"本批无监督信号"。loss 保留计算图以便 backward 安全。
        if (idx_ligand is None or idx_ligand.numel() == 0
                or idx_focal is None or idx_focal.numel() == 0):
            zero_grad = compose_pos.sum() * 0.0
            zero = torch.zeros((), device=compose_pos.device, dtype=torch.float32)
            return (zero_grad, zero, zero, zero, zero, zero, zero, zero)

        # # emebedding
        h_compose = embed_compose(compose_feature, compose_pos, idx_ligand, idx_protein,
                                      self.ligand_atom_emb, self.protein_atom_emb, self.emb_dim)
        # # Encode compose
        h_compose = self.encoder(
            node_attr = h_compose,
            pos = compose_pos,
            edge_index = compose_knn_edge_index,
            edge_feature = compose_knn_edge_feature,
        )   # (N_p+N_l, H)    
        # # 0: frontier atoms of protein
        y_protein_frontier_pred = self.frontier_pred(
            h_compose,
            idx_protein_all_mask
        )
        # # 1: Fontier atoms
        y_frontier_pred = self.frontier_pred(
            h_compose,
            idx_ligand,
        )
        # # 2: Positions relative to focal atoms  `idx_focal`
        relative_pos_mu, abs_pos_mu, pos_sigma, pos_pi  = self.pos_predictor(
            h_compose,
            idx_focal,
            compose_pos,
        )

        # # 3: Element and bonds of the new position atoms
        y_real_pred, edge_pred = self.field(
            pos_query = pos_real,
            edge_index_query = edge_index_real,
            pos_compose = compose_pos,
            node_attr_compose = h_compose,
            edge_index_q_cps_knn = real_compose_knn_edge_index,

            index_real_cps_edge_for_atten = index_real_cps_edge_for_atten,
            tri_edge_index = tri_edge_index,
            tri_edge_feat = tri_edge_feat
        )   # (N_real, num_classes)

        # # fake positions
        y_fake_pred,  _ = self.field(
            pos_query = pos_fake,
            edge_index_query = [],
            pos_compose = compose_pos,
            node_attr_compose = h_compose,
            edge_index_q_cps_knn = fake_compose_knn_edge_index,
        )   # (N_fake, num_classes)

        # # loss
        loss_surf = F.binary_cross_entropy_with_logits(
            input=y_protein_frontier_pred,
            target=y_protein_frontier.view(-1, 1).float()
        ).clamp_max(10.)
        loss_frontier = F.binary_cross_entropy_with_logits(
            input = y_frontier_pred,
            target = y_frontier.view(-1, 1).float()
        ).clamp_max(10.)
        loss_pos = -torch.log(
            self.pos_predictor.get_mdn_probability(abs_pos_mu, pos_sigma, pos_pi, pos_generate) + 1e-16
        ).mean().clamp_max(10.)
        # loss_notpos = self.pos_predictor.get_mdn_probability(abs_pos_mu, pos_sigma, pos_pi, pos_notgenerate).mean()
        loss_cls = self.smooth_cross_entropy(y_real_pred, y_real.argmax(-1)).clamp_max(10.)    # Classes
        loss_edge = F.cross_entropy(edge_pred, edge_label).clamp_max(10.)
        # real and fake loss
        energy_real = -1 *  torch.logsumexp(y_real_pred, dim=-1)  # (N_real)
        energy_fake = -1 * torch.logsumexp(y_fake_pred, dim=-1)   # (N_fake)
        energy_real = torch.clamp_max(energy_real, 40)
        energy_fake = torch.clamp_min(energy_fake, -40)
        loss_real = self.bceloss_with_logits(-energy_real, torch.ones_like(energy_real)).clamp_max(10.)
        loss_fake = self.bceloss_with_logits(-energy_fake, torch.zeros_like(energy_fake)).clamp_max(10.)

        loss = (torch.nan_to_num(loss_frontier)
                    + torch.nan_to_num(loss_pos)
                    + torch.nan_to_num(loss_cls)
                    + torch.nan_to_num(loss_edge)
                    + torch.nan_to_num(loss_real)
                    + torch.nan_to_num(loss_fake)
                    + torch.nan_to_num(loss_surf)
        )
        return loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, torch.nan_to_num(loss_surf)  # loss_notpos
        
    def query_batch(self, pos_query_list, batch, limit=10000):
        pos_query, batch_query = concat_tensors_to_batch(pos_query_list)
        num_query = pos_query.size(0)
        if len(torch.unique(batch_query)) != 1:
            # 不用 assert: -O 运行时会被剥离, 保护会静默失效
            raise NotImplementedError('Modify get_batch_edge to support multiple batches')
        y_cls_all, y_ind_all = [], []
        for pos_query_partial, batch_query_partial in zip(split_tensor_to_segments(pos_query, limit), split_tensor_to_segments(batch_query, limit)):
            PM = batch_intersection_mask(batch.protein_element_batch, batch_query_partial)
            LM = batch_intersection_mask(batch.ligand_context_element_batch, batch_query_partial)
            ligand_context_bond_index, ligand_context_bond_type = get_batch_edge(
                batch.ligand_context_bond_index,
                batch.ligand_context_bond_type,
            )

            y_cls_partial, y_ind_partial, _ = self(
                # Query
                pos_query = pos_query_partial,
                batch_query = batch_query_partial,
                edge_index_query = [],
                # Protein
                protein_pos = batch.protein_pos[PM],
                protein_atom_feature = batch.protein_atom_feature.float()[PM],
                batch_protein = batch.protein_element_batch[PM],
                # Ligand
                ligand_pos = batch.ligand_context_pos[LM],
                ligand_atom_feature = batch.ligand_context_feature_full.float()[LM], 
                batch_ligand = batch.ligand_context_element_batch[LM],
                ligand_context_bond_index = ligand_context_bond_index,
                ligand_context_bond_type = ligand_context_bond_type,
            )
            y_cls_all.append(y_cls_partial)
            y_ind_all.append(y_ind_partial)
        
        y_cls_all = torch.cat(y_cls_all, dim=0)
        y_ind_all = torch.cat(y_ind_all, dim=0)

        lengths = [x.size(0) for x in pos_query_list]
        y_cls_list = split_tensor_by_lengths(y_cls_all, lengths)
        y_ind_list = split_tensor_by_lengths(y_ind_all, lengths)

        return y_cls_list, y_ind_list

    def query_position(self, pos_query, h_compose, compose_pos,
        idx_ligand, ligand_bond_index, ligand_bond_type):
        device = pos_query.device
        #NOTE: Only one parent batch at a time (i.e. batch size = 1)
        edge_index_query = torch.stack(torch.meshgrid(
                torch.arange(len(pos_query), dtype=torch.int64, device=device),
                torch.arange(len(idx_ligand), dtype=torch.int64, device=device),
                indexing=None
            ), dim=0).reshape(2, -1)
        query_compose_knn_edge_index = knn(x=compose_pos, y=pos_query, k=self.config.field.knn, num_workers=16)
        index_real_cps_edge_for_atten, tri_edge_index, tri_edge_feat = self.get_tri_edges(
            edge_index_query = edge_index_query,
            pos_query = pos_query,
            idx_ligand = idx_ligand,
            ligand_bond_index = ligand_bond_index,
            ligand_bond_type = ligand_bond_type
        )
        y_real_pred, edge_pred = self.field(
            pos_query = pos_query,
            edge_index_query = edge_index_query,
            pos_compose = compose_pos,
            node_attr_compose = h_compose,
            edge_index_q_cps_knn = query_compose_knn_edge_index,

            index_real_cps_edge_for_atten = index_real_cps_edge_for_atten,
            tri_edge_index = tri_edge_index,
            tri_edge_feat = tri_edge_feat
        )
        edge_pred = edge_pred.reshape(len(pos_query), len(idx_ligand), self.num_bond_types+1)
        return y_real_pred, edge_pred
        
    def get_tri_edges(self, edge_index_query, pos_query, idx_ligand, ligand_bond_index, ligand_bond_type):
        row, col = edge_index_query
        n_query = pos_query.size(0)
        n_context = len(idx_ligand)
        device = row.device
        # --- 向量化路径: 调用点为 query x context 全连接 (每个 query 的边数恒为 n_context) ---
        # 原实现逐 query 循环 + meshgrid, 复杂度 O(n_query * n_context^2); 此处用
        # repeat_interleave/repeat 一次性构造同一索引对集合, 与原输出逐位一致 (见
        # scripts/test_tri_edges_equiv.py 一致性测试)。不满足规则分组时回退原循环。
        use_fast = False
        if n_context > 0 and row.size(0) == n_query * n_context:
            expect = torch.arange(n_query, device=device).repeat_interleave(n_context)
            use_fast = bool(torch.equal(row, expect))
        if use_fast:
            offs = torch.arange(n_context, device=device)
            ii = offs.repeat_interleave(n_context)
            jj = offs.repeat(n_context)
            base = (torch.arange(n_query, device=device) * n_context).repeat_interleave(n_context * n_context)
            index_real_cps_edge_i = base + ii.repeat(n_query)
            index_real_cps_edge_j = base + jj.repeat(n_query)
        else:
            index_real_cps_edge_i, index_real_cps_edge_j = self._get_tri_edges_loop(row, n_query, device)

        node_a_cps_tri_edge = col[index_real_cps_edge_i]  # the node of tirangle edge for the edge attention (in the compose)
        node_b_cps_tri_edge = col[index_real_cps_edge_j]
        adj_mat = (torch.zeros([n_context, n_context], dtype=torch.long, device=device)
                   - torch.eye(n_context, dtype=torch.long, device=device))
        if ligand_bond_index.numel() > 0:
            adj_mat[ligand_bond_index[0], ligand_bond_index[1]] = ligand_bond_type
        tri_edge_type = adj_mat[node_a_cps_tri_edge, node_b_cps_tri_edge]
        tri_edge_feat = (tri_edge_type.view([-1, 1]) == torch.tensor([[-1, 0, 1, 2, 3]], device=device)).long()

        index_real_cps_edge_for_atten = torch.stack([
            index_real_cps_edge_i, index_real_cps_edge_j  # plus len(real_compose_edge_index_0) for dataloader batch
        ], dim=0)
        tri_edge_index = torch.stack([
            node_a_cps_tri_edge, node_b_cps_tri_edge  # plus len(compose_pos) for dataloader batch
        ], dim=0)
        return index_real_cps_edge_for_atten, tri_edge_index, tri_edge_feat

    @staticmethod
    def _get_tri_edges_loop(row, n_query, device):
        """回退实现: 与原版逐 query 循环语义一致 (非常规边分组时使用)."""
        acc_num_edges = 0
        index_real_cps_edge_i_list, index_real_cps_edge_j_list = [], []
        for node in range(n_query):
            num_edges = (row == node).sum()
            index_edge_i = torch.arange(num_edges, dtype=torch.long, device=device) + acc_num_edges
            index_edge_i, index_edge_j = torch.meshgrid(index_edge_i, index_edge_i, indexing='ij')
            index_real_cps_edge_i_list.append(index_edge_i.flatten())
            index_real_cps_edge_j_list.append(index_edge_j.flatten())
            acc_num_edges += num_edges
        return torch.cat(index_real_cps_edge_i_list, dim=0), torch.cat(index_real_cps_edge_j_list, dim=0)
