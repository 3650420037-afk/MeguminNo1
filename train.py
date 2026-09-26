# import sys
# sys.path.append('.')
import os
import sys

# 核心代码位于 src/ 下(models / utils / scripts), 需先让 src/ 可导入
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from paths import ROOT, SRC, CONFIGS, LOGS, ensure_dir  # noqa: E402

import shutil
import argparse
import json
from tqdm.auto import tqdm
import torch
from torch.nn.utils import clip_grad_norm_
import torch.utils.tensorboard
# import torch_geometric
# assert not torch_geometric.__version__.startswith('2'), 'Please use torch_geometric lower than version 2.0.0'
from torch_geometric.loader import DataLoader

from models.maskfill import MaskFillModelVN
from utils.datasets import *
from utils.transforms import *
from utils.misc import *
from utils.train import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='一键训练/微调 Pocket2Mol 掩码填充模型')
    parser.add_argument('--config', type=str, default=os.path.join(CONFIGS, 'train.yml'))
    parser.add_argument('--device', type=str, default='cuda')
    parser.add_argument('--logdir', type=str, default=LOGS)
    parser.add_argument('--verify-dataset', action='store_true',
                        help='只做数据集验收: 逐样本执行 transform, 全部通过则退出 0')
    args = parser.parse_args()

    # Load configs
    config = load_config(args.config)
    config_name = os.path.basename(args.config)[:os.path.basename(args.config).rfind('.')]
    seed_all(config.train.seed)
    if config.train.use_apex:
        from apex import amp

    # Logging
    ensure_dir(args.logdir)
    log_dir = get_new_log_dir(args.logdir, prefix=config_name)
    ckpt_dir = os.path.join(log_dir, 'checkpoints')
    os.makedirs(ckpt_dir, exist_ok=True)
    logger = get_logger('train', log_dir)
    writer = torch.utils.tensorboard.SummaryWriter(log_dir)
    logger.info(args)
    logger.info(config)
    shutil.copyfile(args.config, os.path.join(log_dir, os.path.basename(args.config)))
    # 归档的是模型架构代码(src/models), 不是 models/ 下的权重文件
    shutil.copytree(os.path.join(SRC, 'models'), os.path.join(log_dir, 'models'))

    # Transforms
    protein_featurizer = FeaturizeProteinAtom()
    ligand_featurizer = FeaturizeLigandAtom()
    masking = get_mask(config.train.transform.mask)
    composer = AtomComposer(protein_featurizer.feature_dim, ligand_featurizer.feature_dim, config.model.encoder.knn)
    
    edge_sampler = EdgeSample(config.train.transform.edgesampler)
    cfg_ctr = config.train.transform.contrastive
    contrastive_sampler = ContrastiveSample(cfg_ctr.num_real, cfg_ctr.num_fake, cfg_ctr.pos_real_std, cfg_ctr.pos_fake_std, config.model.field.knn)
    transform = Compose([
        RefineData(),
        LigandCountNeighbors(),
        protein_featurizer,
        ligand_featurizer,
        masking,
        composer,

        FocalBuilder(),
        edge_sampler,
        contrastive_sampler,
    ])

    # Datasets and loaders
    logger.info('Loading dataset...')
    dataset, subsets = get_dataset(
        config = config.dataset,
        transform = transform,
    )
    train_set, val_set = subsets['train'], subsets['test']

    # 训练前验收门: 把整个数据集过一遍特征化/掩码/构图转换, 任一失败即拒绝开训。
    # 动机: 实测过一次白跑 —— 数据集里有词表外元素的配体, transform 在 iter 43 才
    # 抛 AssertionError('Unexpected elements.') 崩掉训练。数据问题应在开训前暴露。
    if args.verify_dataset:
        logger.info('[Verify] 逐样本执行 transform ...')
        bad = []
        for name, ds in (('train', train_set), ('val', val_set)):
            for i in range(len(ds)):
                try:
                    ds[i]
                except Exception as e:
                    bad.append((name, i, '%s: %s' % (type(e).__name__, str(e)[:120])))
        logger.info('[Verify] 完成: train=%d, val=%d, 失败=%d'
                    % (len(train_set), len(val_set), len(bad)))
        for b in bad[:20]:
            logger.error('[Verify] 失败 %s[%d]: %s' % b)
        if bad:
            raise SystemExit('[Verify] 数据集未通过验收, 共 %d 个样本无法转换; 请先修复数据'
                             '(例如建库时剔除词表外元素/异常几何)' % len(bad))
        print('[Verify] 数据集全部 %d 个样本可正常转换' % (len(train_set) + len(val_set)),
              flush=True)
        raise SystemExit(0)

    follow_batch = []
    collate_exclude_keys = ['ligand_nbh_list']
    train_iterator = inf_iterator(DataLoader(
        train_set, 
        batch_size = config.train.batch_size, 
        shuffle = True,
        num_workers = config.train.num_workers,
        pin_memory = config.train.pin_memory,
        follow_batch = follow_batch,
        exclude_keys = collate_exclude_keys,
    ))
    val_loader = DataLoader(val_set, config.train.batch_size, shuffle=False, follow_batch=follow_batch, exclude_keys = collate_exclude_keys,)

    # Model
    logger.info('Building model...')
    if config.model.vn == 'vn':
        model = MaskFillModelVN(
            config.model, 
            num_classes = contrastive_sampler.num_elements,
            num_bond_types = edge_sampler.num_bond_types,
            protein_atom_feature_dim = protein_featurizer.feature_dim,
            ligand_atom_feature_dim = ligand_featurizer.feature_dim,
        ).to(args.device)
        init_checkpoint = config.train.get('init_checkpoint', None)
        if init_checkpoint:
            logger.info('Loading initialization checkpoint: %s' % init_checkpoint)
            checkpoint = torch.load(init_checkpoint, map_location=args.device)
            model.load_state_dict(checkpoint['model'])
        if config.train.get('freeze_encoder', False):
            for parameter in model.encoder.parameters():
                parameter.requires_grad = False
            logger.info('Frozen encoder parameters for local fine-tuning.')
        # 按参数名子串冻结任意模块(可叠加)。
        # 用途: frontier_pred 决定"哪里还可生长/何时停止", 它的漂移会让分子永不终止
        # (实测: 新数据微调后采样日志 Failed=0, 完成数 63 -> 23/7)。把编码器与
        # frontier_pred 一并冻结, 可让"何时停"保持预训练解, 只让化学头适配新靶点。
        freeze_pats = list(config.train.get('freeze_patterns', []) or [])
        if freeze_pats:
            n_frozen, names = 0, set()
            for pname, parameter in model.named_parameters():
                if any(pat in pname for pat in freeze_pats):
                    parameter.requires_grad = False
                    n_frozen += parameter.numel()
                    names.add(pname.split('.')[0] + '.' + (pname.split('.')[1] if '.' in pname else ''))
            logger.info('Frozen by patterns %s: %.2f M params (%s)'
                        % (freeze_pats, n_frozen / 1e6, sorted(names)))
    print('Num of parameters is', np.sum([p.numel() for p in model.parameters()]))

    # Optimizer and scheduler
    optimizer = get_optimizer(config.train.optimizer, model)
    scheduler = get_scheduler(config.train.scheduler, optimizer)
    if config.train.use_apex:
        model, optimizer = amp.initialize(model, optimizer, opt_level='O1')

    # 梯度累积: 每 grad_accum 个微批才更新一次参数, 等效放大 batch。
    # batch_size=1 时单样本梯度噪声极大(实测训练损失标准差 > 均值, 极差 35.8),
    # 累积是在 8 GB 显存下取得大 batch 梯度信号的最稳做法。
    grad_accum = max(1, int(config.train.get('grad_accum', 1)))

    def train(it):
        # model.train()  has been moved to the end of validation function
        if (it - 1) % grad_accum == 0:
            optimizer.zero_grad()
        batch = next(train_iterator).to(args.device)

        compose_noise = torch.randn_like(batch.compose_pos) * config.train.pos_noise_std
        loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, loss_surf = model.get_loss(

            pos_real = batch.pos_real,
            y_real = batch.cls_real.long(),
            # p_real = batch.ind_real.float(),    # Binary indicators: float
            pos_fake = batch.pos_fake,

            edge_index_real = torch.stack([batch.real_compose_edge_index_0, batch.real_compose_edge_index_1], dim=0),
            edge_label = batch.real_compose_edge_type,
            
            index_real_cps_edge_for_atten = batch.index_real_cps_edge_for_atten,
            tri_edge_index = batch.tri_edge_index,
            tri_edge_feat = batch.tri_edge_feat,

            compose_feature = batch.compose_feature.float(),
            compose_pos = batch.compose_pos + compose_noise,
            idx_ligand = batch.idx_ligand_ctx_in_compose,
            idx_protein = batch.idx_protein_in_compose,

            y_frontier = batch.ligand_frontier,
            idx_focal = batch.idx_focal_in_compose,
            pos_generate=batch.pos_generate,
            idx_protein_all_mask = batch.idx_protein_all_mask,
            y_protein_frontier = batch.y_protein_frontier,

            compose_knn_edge_index = batch.compose_knn_edge_index,
            compose_knn_edge_feature = batch.compose_knn_edge_feature,
            real_compose_knn_edge_index = torch.stack([batch.real_compose_knn_edge_index_0, batch.real_compose_knn_edge_index_1], dim=0),
            fake_compose_knn_edge_index = torch.stack([batch.fake_compose_knn_edge_index_0, batch.fake_compose_knn_edge_index_1], dim=0),
        )
        loss_components = (loss, loss_frontier, loss_pos, loss_cls, loss_edge, loss_real, loss_fake, loss_surf)
        if not all(torch.isfinite(value).all() for value in loss_components):
            logger.warning('[Train] Iter %d skipped non-finite loss batch.' % it)
            return
        # 按累积步数缩放后再反传; 日志里仍记录未缩放的损失, 便于与既有曲线对比
        scaled_loss_value = loss / grad_accum
        if config.train.use_apex:
            with amp.scale_loss(scaled_loss_value, optimizer) as scaled_loss:
                scaled_loss.backward()
        else:
            scaled_loss_value.backward()

        do_step = (it % grad_accum == 0) or (it == config.train.max_iters)
        if do_step:
            orig_grad_norm = clip_grad_norm_(model.parameters(), config.train.max_grad_norm, error_if_nonfinite=True)  # 5% running time
            optimizer.step()
        else:
            orig_grad_norm = float('nan')   # 累积中, 本步不更新

        logger.info('[Train] Iter %d | Loss %.6f | Loss(Fron) %.6f | Loss(Pos) %.6f | Loss(Cls) %.6f | Loss(Edge) %.6f | Loss(Real) %.6f | Loss(Fake) %.6f | Loss(Surf) %.6f  ' % (
            it, loss.item(), loss_frontier.item(), loss_pos.item(), loss_cls.item(), loss_edge.item(), loss_real.item(), loss_fake.item(), loss_surf.item()
        ))
        writer.add_scalar('train/loss', loss, it)
        writer.add_scalar('train/loss_fron', loss_frontier, it)
        writer.add_scalar('train/loss_pos', loss_pos, it)
        writer.add_scalar('train/loss_cls', loss_cls, it)
        writer.add_scalar('train/loss_edge', loss_edge, it)
        writer.add_scalar('train/loss_real', loss_real, it)
        writer.add_scalar('train/loss_fake', loss_fake, it)
        writer.add_scalar('train/loss_surf', loss_surf, it)
        writer.add_scalar('train/lr', optimizer.param_groups[0]['lr'], it)
        writer.add_scalar('train/grad', orig_grad_norm, it)
        writer.flush()

    def validate(it):
        sum_loss, sum_n = np.zeros(5 + 2 + 1), 0   # num of loss
        with torch.no_grad():
            model.eval()
            for batch in tqdm(val_loader, desc='Validate'):
                batch = batch.to(args.device)
                loss_list = model.get_loss(
                    pos_real = batch.pos_real,
                    y_real = batch.cls_real.long(),
                    pos_fake = batch.pos_fake,

                    edge_index_real = torch.stack([batch.real_compose_edge_index_0, batch.real_compose_edge_index_1], dim=0),
                    edge_label = batch.real_compose_edge_type,

                    index_real_cps_edge_for_atten = batch.index_real_cps_edge_for_atten,
                    tri_edge_index = batch.tri_edge_index,
                    tri_edge_feat = batch.tri_edge_feat,

                    compose_feature = batch.compose_feature.float(),
                    compose_pos = batch.compose_pos,
                    idx_ligand = batch.idx_ligand_ctx_in_compose,
                    idx_protein = batch.idx_protein_in_compose,
                    
                    y_frontier = batch.ligand_frontier,
                    idx_focal = batch.idx_focal_in_compose,
                    pos_generate = batch.pos_generate,
                    idx_protein_all_mask = batch.idx_protein_all_mask,
                    y_protein_frontier = batch.y_protein_frontier,

                    compose_knn_edge_index = batch.compose_knn_edge_index,
                    compose_knn_edge_feature = batch.compose_knn_edge_feature,
                    real_compose_knn_edge_index = torch.stack([batch.real_compose_knn_edge_index_0, batch.real_compose_knn_edge_index_1], dim=0),
                    fake_compose_knn_edge_index = torch.stack([batch.fake_compose_knn_edge_index_0, batch.fake_compose_knn_edge_index_1], dim=0),
                )
                sum_loss = sum_loss + np.array([torch.nan_to_num(l).item() for l in loss_list]) 
                sum_n += 1
        avg_loss = sum_loss / sum_n

        if config.train.scheduler.type == 'plateau':
            scheduler.step(avg_loss[0])
        elif config.train.scheduler.type == 'warmup_plateau':
            scheduler.step_ReduceLROnPlateau(avg_loss[0])
        else:
            scheduler.step()

        logger.info('[Validate]  Iter %d | Loss %.6f | Loss(Fron) %.6f | Loss(Pos) %.6f | Loss(Cls) %.6f | Loss(Edge) %.6f | Loss(Real) %.6f | Loss(Fake) %.6f  | Loss(Surf) %.6f' % (
            it, *avg_loss,
        ))
        writer.add_scalar('val/loss', avg_loss[0], it)
        writer.add_scalar('val/loss_fron', avg_loss[1], it)
        writer.add_scalar('val/loss_pos', avg_loss[2], it)
        writer.add_scalar('val/loss_cls', avg_loss[3], it)
        writer.add_scalar('val/loss_edge', avg_loss[4], it)
        writer.add_scalar('val/loss_real', avg_loss[5], it)
        writer.add_scalar('val/loss_fake', avg_loss[6], it)
        writer.add_scalar('val/loss_surf', avg_loss[7], it)
        writer.flush()
        return avg_loss

    # 按验证损失保留最优检查点 + 早停。
    # 背景: 原实现丢弃 validate() 的返回值并无条件保存, 而下游 phase4_retrain.ps1
    # 取"编号最大"的检查点, 于是用的正是最过拟合的那一个 —— 实测验证损失在
    # iter 1200 触底 0.774 后回升到 2000 的 1.898(2.45 倍), 被选中的却是 2000.pt。
    best_val = float('inf')
    best_it = -1
    bad_rounds = 0
    patience = int(config.train.get('early_stop_patience', 0) or 0)   # 0 = 关闭早停
    min_iters = int(config.train.get('early_stop_min_iters', 0) or 0)

    def save_ckpt(path, it):
        torch.save({
            'config': config,
            'model': model.state_dict(),
            'optimizer': optimizer.state_dict(),
            'scheduler': scheduler.state_dict(),
            'iteration': it,
        }, path)

    try:
        model.train()
        for it in range(1, config.train.max_iters+1):
            try:
                train(it)
            except RuntimeError as e:
                logger.error('Runtime Error ' + str(e))
            if it % config.train.val_freq == 0 or it == config.train.max_iters:
                avg_loss = validate(it)
                val_total = float(np.asarray(avg_loss).reshape(-1)[0])
                ckpt_path = os.path.join(ckpt_dir, '%d.pt' % it)
                save_ckpt(ckpt_path, it)
                if val_total < best_val - 1e-6:
                    best_val, best_it, bad_rounds = val_total, it, 0
                    save_ckpt(os.path.join(ckpt_dir, 'best.pt'), it)
                    logger.info('[Best] Iter %d | val_loss %.6f -> 更新 best.pt' % (it, val_total))
                else:
                    bad_rounds += 1
                    logger.info('[Best] Iter %d | val_loss %.6f (最优 %.6f @ iter %d, 连续未改善 %d/%s)'
                                % (it, val_total, best_val, best_it, bad_rounds,
                                   patience if patience else '-'))
                    if patience and it >= min_iters and bad_rounds >= patience:
                        logger.info('[EarlyStop] 连续 %d 次未改善, 于 iter %d 停止; '
                                    '最优 iter %d (val_loss %.6f), 权重见 %s'
                                    % (bad_rounds, it, best_it, best_val,
                                       os.path.join(ckpt_dir, 'best.pt')))
                        break
                # 供下游按验证损失挑选检查点(而不是按编号最大)
                with open(os.path.join(ckpt_dir, 'best.json'), 'w', encoding='utf-8') as f:
                    json.dump({'best_iter': best_it, 'best_val_loss': best_val,
                               'last_iter': it, 'last_val_loss': val_total,
                               'patience': patience, 'min_iters': min_iters,
                               'select_metric': 'val_loss(avg_loss[0])',
                               'note': '下游请使用 best.pt 或 iter=best_iter 的检查点'},
                              f, ensure_ascii=False, indent=1)
                # 逐次验证历史: 供 src/scripts/select_ckpt_by_generation.py 挑候选,
                # 用"生成质量"而非"掩码补全损失"做最终选择(二者可反向, 见技术全书 §4.4)
                with open(os.path.join(ckpt_dir, 'val_history.jsonl'), 'a', encoding='utf-8') as f:
                    f.write(json.dumps({
                        'iter': it,
                        'val_total': val_total,
                        'components': [float(x) for x in np.asarray(avg_loss).reshape(-1)],
                        'is_best': bool(abs(val_total - best_val) < 1e-12),
                    }, ensure_ascii=False) + '\n')
                model.train()
    except KeyboardInterrupt:
        logger.info('Terminating...')
        
