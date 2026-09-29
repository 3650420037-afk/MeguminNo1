import os
import sys

# 核心代码(models / utils)与本文件同在 src/ 下, 需先把 src/ 放入导入路径,
# 这样从仓库根执行 `python src/sample_for_pdb.py` 也能正确导入。
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from paths import (OUTPUTS, PRETRAINED, TARGETS, CONFIGS,  # noqa: E402
                   ensure_dir, missing_hint)

import argparse
import warnings
from easydict import EasyDict
from Bio import BiopythonWarning
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Selection import unfold_entities
from rdkit import Chem

from utils.protein_ligand import PDBProtein
from sample import *    # Import everything from `sample.py`
from utils.guidance import chem_score, logp_to_rank_prob_guided, diversity_penalty


def size_schedule_threshold(base_thr, n_ha, step_i, max_steps, sched):
    """轨迹内尺寸调度（S1）：返回该分子本步应使用的 frontier 阈值。

    ⚠ **机制已按实测修正（2026-09-29 20:40）**：原设计用"按步数线性外推期望值"
    （`target × step/max_steps`）判断落后 —— 但实测 **8 份日志 410 个分子，`重原子数 − 步数 ≡ 0`**：
    采样**每步恰好加 1 个重原子**，真实轨迹是 `1×step`，永远高于"线性外推期望"，
    故原规则在正常 40 步运行里**几乎永不触发**（只会在 max_steps 极小时误触发）。
    → 改为**绝对尺寸恒温器**：**分子只要还小于目标尺寸，就放松阈值鼓励它继续长**；
    超过目标（且开启 tighten）则收紧、让它更容易停。**分子大小 = 停止步数**，故这等价于
    直接把"停止步数"朝目标尺寸推。

    sched（dict；缺省/`target_ha`<=0 时**恒返回 base_thr**，即与历史行为逐位一致）：
      target_ha   目标重原子数（如数据先验 28–29）
      slack       容差（重原子数）：小于 target−slack 才算"偏小"
      relax       偏小时下调的幅度（阈值↓ = 更容易判为 frontier = 继续长）
      tighten     偏大时上调的幅度（默认 0 = 不干预偏大者）
      warmup_frac 前若干比例步数不干预（让模型自行起步）
      floor_steps 若 >0：**前 floor_steps 步一律用 floor_thr** —— 因为实测"每步恰好 1 个重原子"，
                  这等价于一条**硬尺寸下限**（分子至少长到 floor_steps 个重原子），
                  方向**由构造保证**（不像阈值微调那样可能反向，见 §5.39/§5.48）。
      floor_thr   下限期使用的阈值（如 −2.0 = 几乎总能找到 frontier，强制继续长）

    返回 (阈值, 动作)，动作用于统计/留痕：'off'/'warmup'/'floor'/'relax'/'tighten'/'keep'。
    """
    sched = sched or {}
    floor_steps = float(sched.get('floor_steps', 0) or 0)
    if floor_steps > 0 and step_i < floor_steps:
        return float(sched.get('floor_thr', -2.0)), 'floor'
    target = float(sched.get('target_ha', 0) or 0)
    if target <= 0:
        return base_thr, 'off'
    if step_i < float(sched.get('warmup_frac', 0.25)) * max(1, max_steps):
        return base_thr, 'warmup'
    slack = float(sched.get('slack', 4.0))
    if n_ha < target - slack:
        return base_thr - float(sched.get('relax', 0.35)), 'relax'
    tighten = float(sched.get('tighten', 0.0))
    if tighten > 0 and n_ha > target + slack:
        return base_thr + tighten, 'tighten'
    return base_thr, 'keep'


def pdb_to_pocket_data(pdb_path, center, bbox_size):
    center = torch.FloatTensor(center)
    warnings.simplefilter('ignore', BiopythonWarning)
    ptable = Chem.GetPeriodicTable()
    parser = PDBParser()
    model = parser.get_structure(None, pdb_path)[0]

    protein_dict = EasyDict({
        'element': [],
        'pos': [],
        'is_backbone': [],
        'atom_to_aa_type': [],
    })
    for atom in unfold_entities(model, 'A'):
        res = atom.get_parent()
        resname = res.get_resname()
        if resname == 'MSE': resname = 'MET'
        if resname not in PDBProtein.AA_NAME_NUMBER: continue   # Ignore water, heteros, and non-standard residues.

        element_symb = atom.element.capitalize()
        if element_symb == 'H': continue
        x, y, z = atom.get_coord()
        pos = torch.FloatTensor([x, y, z])
        if (pos - center).abs().max() > (bbox_size / 2): 
            continue

        protein_dict['element'].append( ptable.GetAtomicNumber(element_symb))
        protein_dict['pos'].append(pos)
        protein_dict['is_backbone'].append(atom.get_name() in ['N', 'CA', 'C', 'O'])
        protein_dict['atom_to_aa_type'].append(PDBProtein.AA_NAME_NUMBER[resname])
        
    if len(protein_dict['element']) == 0:
        raise ValueError('No atoms found in the bounding box (center=%r, size=%f).' % (center, bbox_size))

    protein_dict['element'] = torch.LongTensor(protein_dict['element'])
    protein_dict['pos'] = torch.stack(protein_dict['pos'], dim=0)
    protein_dict['is_backbone'] = torch.BoolTensor(protein_dict['is_backbone'])
    protein_dict['atom_to_aa_type'] = torch.LongTensor(protein_dict['atom_to_aa_type'])

    data = ProteinLigandData.from_protein_ligand_dicts(
        protein_dict = protein_dict,
        ligand_dict = {
            'element': torch.empty([0,], dtype=torch.long),
            'pos': torch.empty([0, 3], dtype=torch.float),
            'atom_feature': torch.empty([0, 8], dtype=torch.float),
            'bond_index': torch.empty([2, 0], dtype=torch.long),
            'bond_type': torch.empty([0,], dtype=torch.long),
        }
    )
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='在自定义受体口袋内生成分子(口袋引导束搜索); '
                    '一般无需直接调用, 建议用仓库根的 design.py / predict.py')
    parser.add_argument('--pdb_path', type=str,
                        default=os.path.join(TARGETS, '4EIY_A2A受体.pdb'),
                        help='受体 PDB 路径')
    parser.add_argument('--center', type=lambda s: list(map(float, s.split(','))),
                        default=[-0.42, 8.53, 17.13],
                        help='口袋包围盒中心, 格式 x,y,z')
    parser.add_argument('--bbox_size', type=float, default=23.0,
                        help='口袋包围盒边长(埃)')
    parser.add_argument('--config', type=str,
                        default=os.path.join(CONFIGS, 'sample_for_pdb.yml'))
    parser.add_argument('--device', type=str, default='cuda')
    parser.add_argument('--outdir', type=str, default=OUTPUTS)
    parser.add_argument('--save-snapshots', action='store_true',
                        help='额外写出逐步束状态快照 samples_init.pt / samples_<N>.pt。'
                             '单个 170-210 MB, 默认关闭(下游只读 SMILES.txt / SDF/ 与 '
                             'samples_all.pt); 仅在需要中断恢复或调试时打开。')
    args = parser.parse_args()

    miss = missing_hint([('受体结构', args.pdb_path), ('采样配置', args.config)])
    if miss:
        raise SystemExit('[sample_for_pdb] 缺少必需文件:\n' + miss)
    ensure_dir(args.outdir)

    # Load configs
    config = load_config(args.config)
    config_name = os.path.basename(args.config)[:os.path.basename(args.config).rfind('.')]
    seed_all(config.sample.seed)

    # Logging
    log_dir = get_new_log_dir(args.outdir, prefix='sample_for_pdb')
    logger = get_logger('sample', log_dir)
    logger.info(args)
    logger.info(config)
    shutil.copyfile(args.config, os.path.join(log_dir, os.path.basename(args.config)))    
    shutil.copyfile(args.pdb_path, os.path.join(log_dir, os.path.basename(args.pdb_path)))    

    # # Transform
    logger.info('Loading data...')
    protein_featurizer = FeaturizeProteinAtom()
    ligand_featurizer = FeaturizeLigandAtom()
    contrastive_sampler = ContrastiveSample(num_real=0, num_fake=0)
    masking = LigandMaskAll()
    transform = Compose([
        RefineData(),
        LigandCountNeighbors(),
        protein_featurizer,
        ligand_featurizer,
        masking,
    ])
    # # Data
    data = pdb_to_pocket_data(args.pdb_path, args.center, args.bbox_size)
    data = transform(data)

    # # Model (Main)
    logger.info('Loading main model...')
    ckpt = torch.load(config.model.checkpoint, map_location=args.device)
    model = MaskFillModelVN(
        ckpt['config'].model, 
        num_classes = contrastive_sampler.num_elements,
        protein_atom_feature_dim = protein_featurizer.feature_dim,
        ligand_atom_feature_dim = ligand_featurizer.feature_dim,
        num_bond_types = 3,
    ).to(args.device)
    # LoRA 检查点支持（用户路线③，2026-09-29）：若该检查点是用 model.lora 训练的，
    # 权重里会多出 lora_A/lora_B 键 —— 必须先按**同一配置**注入 LoRA，键名才对得上。
    # 注意 LoRA 训练里 `maskfill` 的前向是 base + ΔW，采样必须复现同样的结构，否则会静默丢适配器。
    # 推荐做法仍是先 `src/scripts/merge_lora_checkpoint.py` 把适配器折进基座再出货（导出不依赖 LoRA 代码）。
    _cfg = ckpt.get('config', {}) if isinstance(ckpt, dict) else {}
    _mcfg = _cfg.get('model', {}) if isinstance(_cfg, dict) else getattr(_cfg, 'model', {})
    _lora = ((_mcfg.get('lora', None) if isinstance(_mcfg, dict)
              else getattr(_mcfg, 'lora', None)) or {})
    if _lora.get('enabled', False):
        from utils.lora import apply_lora
        n_lora = apply_lora(model, targets=tuple(_lora.get('targets', ['encoder'])),
                            rank=int(_lora.get('rank', 8)), alpha=float(_lora.get('alpha', 16.0)),
                            verbose=False)
        logger.info('该检查点启用了 LoRA：已按配置注入 %d 层（rank=%d）后加载权重'
                    % (n_lora, int(_lora.get('rank', 8))))
    model.load_state_dict(ckpt['model'])

    # Sampling
    # The algorithm is the same as the one `sample.py`.

    pool = EasyDict({
        'queue': [],
        'failed': [],
        'finished': [],
        'duplicate': [],
        'smiles': set(),
    })
    # # Sample the first atoms
    logger.info('Initialization')
    pbar = tqdm(total=config.sample.beam_size, desc='InitSample')
    atom_composer = AtomComposer(protein_featurizer.feature_dim, ligand_featurizer.feature_dim, model.config.encoder.knn)
    data = transform_data(data, atom_composer)
    # frontier 判定阈值: 显式读取新键 sample.threshold.frontier_threshold, 缺省 0
    # (与 maskfill 的签名默认一致 -> 不写该键时行为与历史完全一致)。
    # 注意不要与 focal_threshold 混用: 后者在原版 sample.py 里是 **focal 概率**阈值,
    # 与 frontier 判定无关 —— 此前误以为改它有效, 实测产出逐位相同。
    _frontier_thr = float(config.sample.threshold.get('frontier_threshold', 0))
    if _frontier_thr:
        logger.info('frontier_threshold = %.4f (非默认, 会改变分子大小分布)' % _frontier_thr)

    # ==== 轨迹内尺寸调度（S1，2026-09-29 新增，默认关闭；关闭时行为与历史逐位一致）====
    # 动机（§5.33/§5.45）：v10 的残余差距 ~44% 由分子尺寸解释，"掉进小模式"是**轨迹早期**就发生的；
    # 固定的 frontier 阈值无法对"已经落后于目标尺寸"的个体做纠正。
    # 做法：每一步、对**每个在长分子**比较"当前重原子数"与"按步数线性外推的目标值"，
    # 落后超过 slack 就把该分子的阈值下调 relax（鼓励继续长）；超前可选上调 tighten。
    # 关键词：sample.threshold.frontier_threshold_schedule.{target_ha,slack,relax,tighten,warmup_frac}
    # 元素表只有 C/N/O/S/Se（重原子），故 ligand_context_element.numel() 即重原子数。
    _sched = config.sample.threshold.get('frontier_threshold_schedule', None) or {}
    _sch_target = float(_sched.get('target_ha', 0) or 0)
    _sch_floor = float(_sched.get('floor_steps', 0) or 0)
    _sch_slack = float(_sched.get('slack', 4.0))
    _sch_relax = float(_sched.get('relax', 0.35))
    _sch_tighten = float(_sched.get('tighten', 0.0))
    _sch_warmup = float(_sched.get('warmup_frac', 0.25))
    _sch_stats = {'relaxed': 0, 'tightened': 0, 'kept': 0, 'floored': 0}
    if _sch_target > 0 or _sch_floor > 0:
        logger.info('尺寸调度开启: target_ha=%.1f slack=%.1f relax=%.2f tighten=%.2f warmup=%.2f '
                    '| floor_steps=%.0f floor_thr=%.2f'
                    % (_sch_target, _sch_slack, _sch_relax, _sch_tighten, _sch_warmup,
                       _sch_floor, float(_sched.get('floor_thr', -2.0))))

    def _thr_for(data_i, step_i):
        """返回该分子本步应使用的 frontier 阈值（调度关闭时恒为 _frontier_thr）。"""
        if _sch_target <= 0 and _sch_floor <= 0:
            return _frontier_thr
        try:
            n_ha = int(data_i.ligand_context_element.numel())
        except Exception:
            return _frontier_thr
        thr_i, action = size_schedule_threshold(
            _frontier_thr, n_ha, step_i, config.sample.max_steps, _sched)
        _sch_stats[{'relax': 'relaxed', 'tighten': 'tightened',
                    'floor': 'floored'}.get(action, 'kept')] += 1
        return thr_i

    init_data_list = get_init(data.to(args.device),   # sample the initial atoms
            model = model,
            transform=atom_composer,
            threshold=config.sample.threshold,
            frontier_threshold=_frontier_thr
    )
    pool.queue = init_data_list
    if len(pool.queue) > config.sample.beam_size:
        pool.queue = init_data_list[:config.sample.beam_size]
        pbar.update(config.sample.beam_size)
    else:
        pbar.update(len(pool.queue))
    pbar.close()

    print_pool_status(pool, logger)
    logger.info('Saving samples...')
    # samples_init.pt 同样是整个束状态(百 MB 级), 且**没有任何下游读取它**
    # (下游只用 SMILES.txt / SDF/ 与 samples_all.pt)。因此默认不写。
    if args.save_snapshots:
        torch.save(pool, os.path.join(log_dir, 'samples_init.pt'))

    # # Sampling loop
    logger.info('Start sampling')
    global_step = 0

    try:
        while len(pool.finished) < config.sample.num_samples:
            global_step += 1
            if global_step > config.sample.max_steps:
                break
            queue_size = len(pool.queue)
            # # sample candidate new mols from each parent mol
            queue_tmp = []
            queue_weight = []
            for data in tqdm(pool.queue):
                nexts = []
                data_next_list = get_next(
                    data.to(args.device), 
                    model = model,
                    transform = atom_composer,
                    threshold = config.sample.threshold,
                    frontier_threshold = _thr_for(data, global_step)
                )

                for data_next in data_next_list:
                    if data_next.status == STATUS_FINISHED:
                        try:
                            rdmol = reconstruct_from_generated_with_edges(data_next)
                            # 可选: 保留模型姿态 + ETKDG 重建内部几何 + 力场优化 (config.sample.relax_output, 默认关)
                            if config.sample.get('relax_output', False):
                                rdmol, info = relax_mol_geometry(
                                    rdmol,
                                    keep_pose=config.sample.get('relax_keep_pose', True),
                                    max_pose_rmsd=config.sample.get('relax_max_pose_rmsd', 2.0),
                                )
                                if info is not None:
                                    data_next.pose_rmsd = info["pose_rmsd"]
                                    data_next.relax_energy = info["relax_energy"]
                                    data_next.strain_per_heavy = info["strain_per_heavy"]
                            data_next.rdmol = rdmol
                            mol = Chem.MolFromSmiles(Chem.MolToSmiles(rdmol))
                            smiles = Chem.MolToSmiles(mol)
                            data_next.smiles = smiles
                            if smiles in pool.smiles:
                                logger.warning('Duplicate molecule: %s' % smiles)
                                pool.duplicate.append(data_next)
                            elif '.' in smiles:
                                logger.warning('Failed molecule: %s' % smiles)
                                pool.failed.append(data_next)
                            else:   # Pass checks
                                logger.info('Success: %s' % smiles)
                                pool.finished.append(data_next)
                                pool.smiles.add(smiles)
                        except MolReconsError:
                            logger.warning('Ignoring, because reconstruction error encountered.')
                            pool.failed.append(data_next)
                    elif data_next.status == STATUS_RUNNING:
                        nexts.append(data_next)

                queue_tmp += nexts
                if len(nexts) > 0:
                    queue_weight += [1. / len(nexts)] * len(nexts)
            if len(queue_tmp) == 0:
                logger.warning('No candidates remain; stopping sampling.')
                break
            # # random choose mols from candidates
            guided_cfg = config.sample.get('guided', None)
            if guided_cfg is not None and guided_cfg.get('enabled', False):
                # 方向1: 引导束搜索 -- QED/SA 化学分数 + (可选)Tanimoto 多样性惩罚
                chem = [chem_score(p, guided_cfg.get('qed_w', 1.0), guided_cfg.get('sa_w', 1.0)) for p in queue_tmp]
                div_pen = diversity_penalty(queue_tmp, pool.smiles,
                                            guided_cfg.get('diversity_w', 0.0))
                chem = [c - d for c, d in zip(chem, div_pen)]
                prob = logp_to_rank_prob_guided(np.array([p.average_logp[2:] for p in queue_tmp]), chem,
                                                queue_weight, guided_cfg.get('lam', 1.0))
            else:
                prob = logp_to_rank_prob(np.array([p.average_logp[2:] for p in queue_tmp]), queue_weight)  # (logp_focal, logpdf_pos), logp_element, logp_hasatom, logp_bond
            n_tmp = len(queue_tmp)
            next_idx = np.random.choice(np.arange(n_tmp), p=prob, size=min(config.sample.beam_size, n_tmp), replace=False)
            pool.queue = [queue_tmp[idx] for idx in next_idx]

            print_pool_status(pool, logger)
            # 逐步快照**默认不写**: 每个快照会把整个束状态(含蛋白特征)序列化,
            # 实测单个 170-210 MB —— 一次 60 分子的采样就能产出数 GB, 累积曾吃掉
            # 仓库 43 GB(outputs 下 412 个 samples_<N>.pt), 有爆盘风险。
            # 官方 src/sample.py 的对应行本就是注释掉的, 这里对齐官方行为。
            # 确需中断恢复时用 --save-snapshots 显式打开。
            if args.save_snapshots:
                torch.save(pool, os.path.join(log_dir, 'samples_%d.pt' % global_step))
    except KeyboardInterrupt:
        logger.info('Terminated. Generated molecules will be saved.')

    if _sch_target > 0 or _sch_floor > 0:
        logger.info('尺寸调度统计: 下限期 %d 次 | 下调 %d 次 | 上调 %d 次 | 保持 %d 次 '
                    '(target_ha=%.1f relax=%.2f floor_steps=%.0f)'
                    % (_sch_stats['floored'], _sch_stats['relaxed'], _sch_stats['tightened'],
                       _sch_stats['kept'], _sch_target, _sch_relax, _sch_floor))

    # # Save sdf mols
    sdf_dir = os.path.join(log_dir, 'SDF')
    os.makedirs(sdf_dir)
    with open(os.path.join(log_dir, 'SMILES.txt'), 'a') as smiles_f:
        for i, data_finished in enumerate(pool['finished']):
            smiles_f.write(data_finished.smiles + '\n')
            rdmol = data_finished.rdmol
            Chem.MolToMolFile(rdmol, os.path.join(sdf_dir, '%d.sdf' % i))
            
    torch.save(pool, os.path.join(log_dir, 'samples_all.pt'))
