# 随机种子记录

> 附件5 第二节要求提交「随机种子记录」。本文件汇总各类运行使用的随机种子，
> 便于逐项复现。种子值取自仓库内**实际使用的配置文件与脚本默认值**，非事后追认。

## 1. 数据集划分

| 用途 | 种子 | 出处 |
|---|---|---|
| 单靶点数据集划分（a2a / b2ar / d3 / 5ht2b） | **2021** | `src/scripts/build_gpcr_dataset.py --seed`（默认 2021），`--val-ratio 0.2` |
| 数据扩充建库（实验产物，已回退） | **2021** | `src/scripts/build_gpcr_v3.py --seed`（默认 2021） |
| 合并集划分（`gpcr_multitarget*`） | **无随机** | `src/scripts/merge_gpcr_datasets.py` 不做打乱，`train=前 80 %` / `val=后 20 %` |

> 合并集的验证集实为**按靶点留出**（5HT2B 全部在验证集），详见 `data/README.md` 第三节。
> 该项按位置切分，因此天然可复现，不需要种子。

## 2. 分子生成 / 采样

| 配置或入口 | 种子 | 说明 |
|---|---|---|
| `configs/sample_for_pdb_guided_l3.yml` 等 16 个引导配置 | **2024** | 本项目主用配置（λ=3） |
| `configs/sample.yml`、`configs/sample_for_pdb.yml` | **2020** | 测试集/PDB 采样默认 |
| A/B 双种子复现（`..._ab_*_s2025.yml`、`..._5010050_*_s2025.yml`） | **2025** | 与 2024 对照，验证 λ 增益不是单次随机 |
| 第三组 A/B（`..._ab_*_s2026.yml`） | **2026** | 同上 |
| `configs/sample_for_pdb_a2a_trial.yml` | **2021** | 早期试用配置 |
| `design.py` / `predict.py` 默认 `--seed` | **2024** | 与主用引导配置一致 |
| `configs/gen_batch.yml` | **8229** | 批量生成 |
| `configs/overnight_batch.yml` | **7113** | 夜间批跑 |
| `src/scripts/e2e_full_run.py` | **2024** | 端到端自检 |

## 3. 训练 / 微调

| 配置 | 种子 | 说明 |
|---|---|---|
| `configs/train_multitarget_v2.yml` | **2022** | *交付用*多靶点训练配置 |
| `configs/train.yml`、`configs/train_a2a_trial.yml` | **2021** | 早期训练 |
| `configs/train_multitarget_trial.yml` | **2023** | 多靶点试用 |
| `configs/abft_finetuned.yml`、`configs/re_abft_finetuned.yml` | **2024** | **已回退**的微调实验（未通过 A/B，见 ModelCard.md） |

## 4. 分析与诊断

| 脚本 | 种子 | 说明 |
|---|---|---|
| `src/scripts/diagnose_dataset.py`（掩码行为抽样） | **2021** | `simulate_mask_reality(..., seed=2021)` |
| `src/scripts/diagnose_dataset.py`（理论闭式概率） | **2021** | `p_empty_context_theory(..., rng)` |
| `src/scripts/build_target_registry.py` | 无随机 | 从结构文件确定性提取 |

## 5. 实现方式与确定性边界

- 统一种子设置函数：`src/utils/misc.py::seed_all(seed)`，训练与采样入口在
  加载配置后立即调用（`seed_all(config.train.seed)` / `seed_all(config.sample.seed)`）。
- 该函数设置 Python / NumPy / PyTorch（含 CUDA）的种子。

**确定性边界（如实说明）**：

1. **建库抽样**：`src/scripts/build_gpcr_v3.py` 的条目偏移原用内置 `hash(str)`，
   其进程随机盐会使同一份输入每次建出不同数据集；**已修复为 `zlib.crc32`**，
   修复后两次运行逐位一致（见 `docs/数据建库穿模修复_验证报告.md`）。
2. **生成过程**：束搜索的采样含随机性，固定种子可复现**同一次**运行，
   但不同 GPU / 驱动 / PyTorch 版本下浮点归约顺序不同，可能出现个体分子差异。
   因此对外报告的质量指标一律采用**中位数**等群体统计量，而非单个分子。
3. **对接（AutoDock Vina）**：Vina 自身默认随机，复现需显式传 `--seed`
   （见 `src/scripts/docking_pipeline.py`）；报告中对接数值为全库统计量。
4. **训练**：`seed_all` 固定种子后，同硬件同版本可复现；跨版本不保证逐位一致。

## 6. 与报告数字的对应

`README.md`「主要成果」中的关键数字对应的运行：

| 报告数字 | 配置 | 种子 |
|---|---|---|
| 引导 λ=3 使 QED 中位 0.667 → 0.803 | `sample_for_pdb_guided_l3.yml` vs 关闭引导 | 2024 |
| 同上，双种子复现 | 同上 | 2024 / 2025 |
| A2A 库 2,931 分子 | `overnight_a2a` 系列批次 | 见各批次配置 |
| 训练冒烟（`e2e_full_run.py` 第 7 阶段） | `train_multitarget_v2.yml`（`max_iters: 2`） | 2022 |
