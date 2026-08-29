Pocket2Mol 竞赛项目 — 参赛代码 README

1. 项目概述
    本仓库包含"口袋引导的 GPCR 靶向分子生成"项目的全部代码：基线模型 Pocket2Mol（ICML 2022）、
    Windows/PyTorch 2.6 兼容补丁、GPCR 微调数据管线、方向1 引导束搜索（本项目核心创新）、
    构象审计/候选评估/化合物库构建等工程脚本。

2. 环境要求
    Windows 11 + NVIDIA GPU（8GB 显存即可，示例在 RTX 4070 Laptop 上验证）
    Miniconda，Python 3.10
    依赖：torch 2.6.0+cu124, torch_geometric 2.8.0, torch_cluster, torch_scatter,
          rdkit 2026.03, biopython, lmdb, easydict, pyyaml, tqdm, tensorboard
    预训练权重：ckpt/pretrained_Pocket2Mol.pt（42.84MB，见 Releases）

3. 兼容性补丁（相对原版 Pocket2Mol 的差异）
    utils/misc.py       PyTorch 2.6 torch.load weights_only 白名单 + PyG WITH_KNN 强制开启
    utils/transforms.py knn/radius 改用 torch_cluster（Windows pyg-lib 缺算子）
    models/common.py    同上
    utils/guidance.py   新增：QED/SA 引导分数 + Tanimoto 多样性惩罚 + rdkit.six 兼容 shim
    train.py            新增 init_checkpoint / freeze_encoder / 非有限 loss 跳过
    models/position.py  logsigma clamp 防数值爆炸
    models/sample.py    average_logp 五元组统一（修复排序阶段 ValueError）+ 空键保护

4. 快速开始（单口袋采样）
    python sample_for_pdb.py --pdb_path example/4yhj.pdb --center " 32.0,28.0,36.0" ^
        --config configs/sample_for_pdb.yml --outdir outputs/demo
    注意：--center 的值首字符必须是一个空格（与原版约定一致）。
    输出：<outdir>/sample_for_pdb_<时间戳>/{SMILES.txt, SDF/*.sdf, samples_all.pt, log.txt}

5. 方向1：引导束搜索（本项目核心算法创新）
    在 configs 中加入以下配置即启用（默认关闭时与原版行为完全一致）：
        sample:
          guided:
            enabled: true
            qed_w: 1.0        # QED 权重
            sa_w: 1.0         # 合成可及性权重（归一化后参与打分）
            lam: 3.0          # 引导强度（λ 扫描确定的最优值，见 docs/）
            diversity_w: 0.5  # 对成品库的 Tanimoto 最大相似度惩罚（0 为关闭）
    原理：束排序概率从 P ∝ exp(ΣlogP_model)+1 修改为
          P ∝ (exp(ΣlogP_model)+1) · exp(λ·(QED + (1-SA/10)) - λ·w_div·maxTanimoto)
    实验结果（A2A 受体 4EIY，50样本/100束宽/50步，两个种子）：
          基线      QED 0.667/0.688，闭合 62/52
          引导 λ=3  QED 0.803/0.821，闭合 53/53，SA 与 MW 不受损
    详见 scripts/guidance_ab_summary.md（λ 扫描与复验数据）。

6. GPCR 微调数据管线（方向2）
    scripts/build_gpcr_dataset.py  SMILES CSV + 受体 PDB → PDB/SDF/index/split 微调数据集
    scripts/merge_gpcr_datasets.py 多靶点数据集合并
    scripts/audit_conformers.py    构象质量审计（键长/原子重叠/共面塌缩）
    scripts/eval_candidates.py     未完成候选的质量评估
    configs/train_a2a_local.yml    冻结 encoder 的局部微调配置（init_checkpoint 指向预训练权重）
    python train.py --config configs/train_a2a_local.yml

7. 工程脚本
    scripts/build_library.py     产出过滤入库（PAINS/Brenk/理化过滤 + SQLite + SDF 库）
    scripts/docking_pipeline.py  AutoDock Vina 对接管线（方向3）
    scripts/overnight_a2a.ps1    批量过夜采样驱动

8. 目录结构
    sample.py / sample_for_pdb.py / train.py   官方入口（仅少量兼容性修改）
    models/          等变图神经网络（MaskFillModelVN）
    utils/           数据变换、重建、引导模块
    configs/         全部实验配置
    data/            GPCR 微调数据集（a2a/b2ar/d3/5ht2b/multitarget_v2）
    outputs/         采样产出（SMILES + SDF + 快照）
    evaluation/      官方评估工具（含 SA score）

9. 可复现性说明
    全部实验配置（seed/束宽/步数/λ）均在 configs/ 内固化；每个采样会话自动
    归档配置副本与完整日志（log.txt）；数据库构建脚本记录过滤漏斗与来源追踪。
