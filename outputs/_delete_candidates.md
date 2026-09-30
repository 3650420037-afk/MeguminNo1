# 待删候选清单（**未删除任何文件**）

## A. B2AR 相关（生成入口已禁用 → 库与中间产物可清）
- 条目数 **29**，文件 **11568**，逻辑大小 **10.82 GB**，实际占用 **0.07 GB**
- 性质：B2AR 药库/中间产物；**回归证据（CSV/日志）建议保留**，见 D
  - `data\gpcr_b2ar` — 0.00 GB
  - `data\gpcr_b2ar_name2id.pt` — 0.00 GB
  - `data\gpcr_b2ar_processed.lmdb` — 10.74 GB
  - `data\gpcr_b2ar_processed.lmdb-lock` — 0.00 GB
  - `other\library_backups\20260930_001605_v2_libraries\B2AR` — 0.00 GB
  - `other\library_backups\20260930_094311_after_swap\B2AR` — 0.00 GB
  - `outputs\ab_4t_B2AR_thrm05.csv` — 0.00 GB
  - `outputs\ab_4t_B2AR_thrm05.csv` — 0.00 GB
  - `outputs\ab_dock_v3_B2AR` — 0.01 GB
  - `outputs\ab_dock_v3_B2AR` — 0.01 GB
  - `outputs\ab_dock_v3_B2AR.csv` — 0.00 GB
  - `outputs\ab_dock_v3_B2AR.csv` — 0.00 GB
  - …另有 17 项

## B. 历史药库（官方权重库 5,864 分子 + v2 时代四靶点库备份）
- 条目数 **7**，文件 **56494**，逻辑大小 **0.29 GB**，实际占用 **0.40 GB**
- 性质：**交付历史与对照证据**；删除后 `official_weight_library` 的对账数字无法本地复算
  - `other\library_backups\20260930_001605_v2_libraries` — 0.01 GB
  - `other\library_backups\20260930_094311_after_swap` — 0.02 GB
  - `other\outputs_history\docking_5ht2b_library` — 0.04 GB
  - `other\outputs_history\docking_a2a_library` — 0.10 GB
  - `other\outputs_history\docking_b2ar_library` — 0.04 GB
  - `other\outputs_history\docking_d3_library` — 0.04 GB
  - `results\official_weight_library` — 0.04 GB

## C. 纯中间文件（采样快照 / 对接中间 / 训练检查点 / 缓存）
- 条目数 **267**，文件 **204054**，逻辑大小 **45.22 GB**，实际占用 **44.97 GB**
- 性质：**可重建或已无用**：`samples_*.pt` 是无下游读取的束状态快照；对接 pdbqt 可重跑；`__pycache__` 自动重建；训练检查点除已抽取的 v10@5500 外均非交付物（但为 §5.33 检查点分析的证据）
  - `__pycache__` — 0.00 GB
  - `logs\_smoke_lora_2026_09_29__19_48_15\checkpoints` — -
  - `logs\_smoke_lora_2026_09_29__19_48_15\models\__pycache__` — 0.00 GB
  - `logs\_smoke_lora_2026_09_29__19_48_15\models\encoders\__pycache__` — 0.00 GB
  - `logs\_smoke_lora_2026_09_29__19_48_15\models\fields\__pycache__` — 0.00 GB
  - `logs\train_gpcr_arch_v7_2026_09_27__10_44_03\checkpoints` — 0.99 GB
  - `logs\train_gpcr_arch_v7_2026_09_27__10_44_03\models\__pycache__` — 0.00 GB
  - `logs\train_gpcr_arch_v7_2026_09_27__10_44_03\models\encoders\__pycache__` — 0.00 GB
  - `logs\train_gpcr_arch_v7_2026_09_27__10_44_03\models\fields\__pycache__` — 0.00 GB
  - `logs\train_gpcr_arch_v8_2026_09_27__11_30_12\checkpoints` — 0.52 GB
  - `logs\train_gpcr_arch_v8_2026_09_27__11_30_12\models\__pycache__` — 0.00 GB
  - `logs\train_gpcr_arch_v8_2026_09_27__11_30_12\models\encoders\__pycache__` — 0.00 GB
  - …另有 255 项

## D. 巨型缓存/归档（收益最大，但需你确认）
- 条目数 **11**，文件 **92769**，逻辑大小 **296.05 GB**，实际占用 **189.52 GB**
- 性质：`*_processed.lmdb` = **训练时自动生成的数据集缓存**（lmdb map_size 固定 10 GB/个），删除后下次用该数据集会自动重建（需原始数据在库，均在）；`other/outputs_history` = 历史 outputs 归档
  - `data\gpcr_b2ar_processed.lmdb` — 10.74 GB
  - `data\gpcr_dock_v1_processed.lmdb` — 10.74 GB
  - `data\gpcr_ft_v3_processed.lmdb` — 10.74 GB
  - `data\gpcr_ft_v3c_processed.lmdb` — 10.74 GB
  - `data\gpcr_mass_v1_processed.lmdb` — 10.74 GB
  - `data\gpcr_multi_v1_processed.lmdb` — 10.74 GB
  - `data\gpcr_multitarget_processed.lmdb` — 10.74 GB
  - `data\gpcr_multitarget_v2_processed.lmdb` — 10.74 GB
  - `data\gpcr_v3_clean_processed.lmdb` — 10.74 GB
  - `data\gpcr_v3_merged_processed.lmdb` — 10.74 GB
  - `other\outputs_history` — 188.68 GB

本清单涉及逻辑总量约 **352.39 GB**（实际占用更少，lmdb 多为稀疏分配）。

**等你确认后我再删除**；建议顺序：C（纯中间）→ A（B2AR）→ D（缓存/归档）→ B（历史药库，最后）。
