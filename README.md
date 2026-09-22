# 7-eonmol · GPCR 靶向分子生成器

> **基于 Pocket2Mol 的引导束搜索改进** · 四靶点 6,864 个类药化合物 · 六层计算评估 · 全流程开源

[Pocket2Mol](https://arxiv.org/abs/2205.07249)（ICML 2022）用等变图神经网络做口袋条件的三维分子生成。
本仓库在此基础上做了**采样层面的改进**：用"引导束搜索"把生成分布推向类药区域，
并配套了从口袋结构到候选分子的**全自动、可追溯**流水线。

**核心数字**：引导束搜索使 QED 中位从 0.667 提升到 **0.803（+20%）**；
四靶点构建 **6,864** 个类药化合物；A2A 库 2,931 分子对接均值 **−9.68 kcal/mol**，
其中 1,137 个强于 −10；对已知活性分子对照集（737 个）富集 **5.80×**（ROC-AUC 0.685）。

> 本项目定位：**造引擎，不挖矿**——把"从口袋结构生成可用候选"做成可复用的引擎，
> 而不是重复筛一遍已有化合物库。

---

## 目录

- [一页速览](#一页速览)
- [1. 任务与定位](#1-任务与定位)
- [2. 主要成果](#2-主要成果)
- [3. 快速开始](#3-快速开始)
- [4. 方法：我们改了什么](#4-方法我们改了什么)
- [5. 端到端流程](#5-端到端流程)
- [6. 代码地图](#6-代码地图)
- [7. 复现步骤](#7-复现步骤)
- [8. 数据资产](#8-数据资产)
- [9. 常见问题](#9-常见问题)
- [10. 局限与未采纳路线](#10-局限与未采纳路线)
- [11. 许可与引用](#11-许可与引用)
- [12. 联系方式](#12-联系方式)

---

## 一页速览

| 维度 | 结论 |
|---|---|
| 基础模型 | Pocket2Mol 官方预训练权重（42.8 MB）；等变消息传递网络（MPNN，6 层 / 256 标量通道 / 64 方向通道）+ MDN 位置头（3 分量） |
| **核心改进** | **引导束搜索**：束排序概率乘化学分数增益 `exp(λ·chem)`；λ=3 时 QED 中位 0.667 → 0.803（**+20%**），双随机种子复现 |
| 伴随改进 | Tanimoto 多样性惩罚（MMR 式）、生成构象力场精修、入库漏斗、全自动对接与分析管线 |
| 主库（A2A） | **2,931** 分子 / **1,960** 唯一 Murcko 骨架 / QED 中位 **0.825** / SA 中位 **3.33** / MW 中位 **362 Da** |
| 四靶点合计 | 原始产出 10,004 → 过滤入库 **6,864** |
| 对接（A2A） | 全库对接成功；均值 **−9.68** kcal/mol；≤−10 有 **1,137** 个（38.8%）；最强 **−13.51** |
| 逆向阳性对照 | 737 个已知 A2A 活性分子走同一管线，均值 −9.32 → **生成库优于已知集** |
| 富集验证 | ROC-AUC **0.685**；对接分 ≤−12 富集 **5.80×**；top-1% 富集 2.24× |
| Top50 候选 | Vina −13.42 ~ −11.59；对已知集最大相似度仅 0.187–0.321；**全部 HIGH_NOVEL** |
| 成药性 | ADMET 规则分级 A **676** / B 1,940 / C 315 |
| 多目标 | 五目标 Pareto 前沿 **91** 个非支配解（旧 Top50 仅 27 个在沿上） |
| 选择性 | A2A vs A1 反向选择性 SI 中位 **+4.0**，59% 分子 >0，382 个 >25 |
| 姿态置信 | 64% 分子的最优/次优构象能量差 < 0.5 kcal/mol |
| 综合推荐 | 六维加权 Tier1 **166** 个 / Tier2 454 个 |
| 交付形态 | 桌面 GUI（免安装 exe，启动 1.3–2.6 秒）+ 5 GB 自包含部署包 + 全部脚本与验证报告 |

环境实测：Python 3.10.20 · PyTorch 2.6.0+cu124 · RDKit 2026.03.5 · PyG 2.8.0 · 单卡 RTX 4070 Laptop 8 GB。

---

## 1. 任务与定位

**赛道 3 · 离子通道/GPCR 小分子药物虚拟筛选 · 子任务 2**：
给定 GPCR 靶点口袋结构，生成**结构新颖、具备类药性**的小分子候选，并证明流程可执行、结果可复现。

拆成工程要求：① 口袋条件生成（分子必须"长在口袋里"）；② 类药性（QED/SA/理化达标，而非仅几何合理）；
③ 新颖性（与已知配体骨架不重合）；④ 可验证（每个数字都能用脚本重算）。

**技术路线取舍**

| 路线 | 代表 | 取舍 |
|---|---|---|
| 基于配体 | JT-VAE / REINVENT | 需要大量已知活性配体；GPCR 多靶点数据稀疏 |
| **基于结构** | **口袋条件 3D 生成（Pocket2Mol）** | 直接利用晶体结构，无需配体数据，原生输出三维构象 |
| 基于对接筛选 | 虚拟筛选已有库 | 依赖库本身质量，无法创造新化学空间 |

内置四个靶点（蛋白结构随仓库提供）：

| 靶点 | PDB | 口袋中心 (x,y,z) | 适应症 |
|---|---|---|---|
| A2A 腺苷受体 | 4EIY | −0.4, 8.5, 17.1 | 帕金森病 / 肿瘤免疫 / 炎症 |
| β2 肾上腺素受体 | 2RH1 | −29.5, 9.2, 6.9 | 哮喘 / 心血管疾病 |
| D3 多巴胺受体 | 3PBL | 0.085, −14.828, 10.432 | 精神分裂症 / 成瘾 |
| 5-HT2B 血清素受体 | 4IB4 | 22.448, 18.284, 11.726 | 偏头痛 / 肺动脉高压 |

---

## 2. 主要成果

### 2.1 生成质量：引导束搜索 +20% QED

严格 A/B（同种子、同 `50 样本 / 100 束宽 / 50 步`，仅切换引导开关）：

| 配置 | QED 中位 | 说明 |
|---|---|---|
| 基线（原版束搜索） | 0.667 | `guided.enabled: false`，行为与原版逐位一致 |
| **引导 λ=3** | **0.803** | +20%，两个独立随机种子（2024 / 2025）复现 |

### 2.2 化合物库

| 指标 | 数值 |
|---|---|
| A2A 库 | 2,931 分子 / 1,960 唯一 Murcko 骨架 |
| QED 中位（最高） | 0.825（0.946）；76% 分子 QED ≥ 0.8 |
| SA 中位 | 3.33（全部落在可合成区间 1.5–5.9） |
| MW 中位 | 362 Da（口服药物甜区） |
| 过滤漏斗 | 原始 10,004 → PAINS/BRENK → 理化 → QED → SA → 元素/环数 → **6,864** |

### 2.3 对接与新颖性

- **全库对接**：AutoDock Vina 1.2.5 + Open Babel 3.2.1 便携部署；受体自动去水/去共晶/转换；
  平均 **−9.68 kcal/mol**，≤−10 有 1,137 个（38.8%），≤−11 有 311 个，最强 −13.51。
- **逆向阳性对照**：737 个已知 A2A 高活性分子同管线均值 −9.32 → 管线可信且生成库更优。
- **新颖性**：Top50 对 362 个已知骨架零重合，最大 Tanimoto 相似度仅 0.187–0.321，全部判为 HIGH_NOVEL。

### 2.4 六层计算评估

| 评估层 | 方法 | 结论 |
|---|---|---|
| 富集能力 | decoy 对照 + ROC | AUC 0.685；≤−12 富集 5.80× |
| 姿态置信 | 多构象重打分，最优/次优能量差 | 64% 分子 < 0.5 kcal/mol（Top50 达 29% 高度可信 vs 全库 6%） |
| 成药风险 | 规则化 ADMET（hERG/肝毒/透膜等） | A 676 / B 1,940 / C 315，标注而非删除 |
| 多目标权衡 | 五目标 Pareto（打分/QED/SA/新颖性/ADMET） | 前沿 91 个非支配解 |
| 靶点选择性 | A2A vs A1 反向选择性 SI | 中位 +4.0，59% >0，382 个 >25（灵敏度检验 p=0.0102） |
| 综合推荐 | 六维加权（对接/QED/新颖/ADMET/选择性/可测性） | Tier1 166 / Tier2 454 |

### 2.5 交付形态

- **桌面应用**：Tkinter GUI（4 个标签页 + 浅/深色主题 + 内置环境自检），PyInstaller 打包为免安装 exe，
  启动 1.3–2.6 秒（已做启动优化：环境自检异步化、依赖裁剪、Tcl 数据裁剪）。
- **自包含部署包**：5 GB 环境 + 仓库 + exe，目标机无需安装 Python/CUDA 工具链。
- **全部脚本与报告**：从数据构建到 11+ 项分析的脚本化产物，可逐文件追溯。

---

## 3. 快速开始

### 方式一：桌面应用（推荐给非编程用户）

**便携版（启动最快）**

1. 解压 `7-eonmol_便携版.zip`，双击 `7-eonmol.exe`
2. 首次运行若提示"Windows 已保护你的电脑"，点「更多信息 → 仍要运行」

**自包含部署包（目标机无 Python 环境时）**

1. 整个 `7-eonmol` 文件夹拷到目标机（路径建议无中文）
2. 双击 `一键配置.bat`（把内置环境路径适配到本机，只需一次）
3. 双击 `7-eonmol.exe`

> GUI 自身零依赖；分子生成调用本机 Python 环境（默认 `D:\Miniconda3\envs\Pocket2Mol`），
> 可在「环境设置」页修改并保存。

### 方式二：源码运行

```bash
# 1) 环境（详见下方"复现步骤"）
conda create -n pocket2mol python=3.10
conda activate pocket2mol
pip install -r requirements.txt          # 含 torch / PyG / RDKit / PyYAML 等

# 2) 下载官方预训练权重到 ckpt/
#    https://drive.google.com/drive/folders/1KfdOczjUPITPhIvCuBmnj4xFTV-iI2xB
#    放到 ckpt/pretrained_Pocket2Mol.pt

# 3) 启动界面
python gui/app_pro.py                    # 可选 --tab 0-3 --theme dark
```

### 方式三：命令行采样（科研用法，全参数可控）

```bash
python sample_for_pdb.py \
  --pdb_path ./targets/4EIY_A2A受体.pdb \
  --center " -0.4,8.5,17.1" \
  --config ./configs/sample_for_pdb_guided_l3.yml \
  --outdir ./outputs
```

> 注意 `--center` 首个数值前的**前导空格必须保留**（负数参数解析约定），例如 `" -0.4,8.5,17.1"`。

---

## 4. 方法：我们改了什么

设计原则：**零侵入**——所有改进以"新增模块 + 配置开关"实现；
`config.sample.guided` 不存在时，程序行为与原版**逐位一致**。任何改进都可独立关闭、独立证伪。

### 4.1 引导束搜索（核心贡献）★

原版束搜索只用模型自身概率排序：

```
prob ∝ (exp(Σ log p_model) + 1) · weight
```

我们在每步束排序时乘上化学分数增益：

```
prob ∝ (exp(Σ log p_model) + 1) · weight · exp( λ · [ QED + (1 − SA/10) − w_div · maxTanimoto ] )
```

- `QED`：类药性定量估计；`SA`：合成可及性（越低越易合成）；
- `maxTanimoto`：与已生成分子的最大相似度 → 惩罚重复结构；
- `λ`：引导强度（默认 3.0，经 λ ∈ {0.3, 1, 2, 3, 5} 消融实验确定）；
- 实现：`utils/guidance.py`（`chem_score` / `logp_to_rank_prob_guided` / `diversity_penalty`）。

**与"事后过滤"的本质区别**：过滤是生成完再丢弃（浪费算力）；
引导是让**同等算力下的产出分布本身更优**。

### 4.2 多样性惩罚（MMR 式）

对每个候选计算其与当前池内分子的最大 Tanimoto 相似度，按 `exp(−w_div · maxTanimoto)` 降权，
防止搜索塌缩到单一骨架。成品库骨架数 1,960（/2,931 分子）即为该机制的间接证据。

### 4.3 生成构象精修

`utils/reconstruct.py::relax_mol_geometry`：在保持姿态的前提下，
用 ETKDG 内部几何 + MMFF94s 优化生成构象，并约束与原始姿态的 RMSD（`max_pose_rmsd`）、
输出应变能 `strain_per_heavy` 供后续审计。避免"几何看起来合理但力场下崩掉"的构象进入下游对接。

### 4.4 微调路线（方向 2）—— 诚实记录：未晋级

构建了四靶点 479 对微调数据集（`data/gpcr_multitarget_v2`），2000 步微调后 A/B 对比：

| 配置 | 产出分子数 | QED 中位 |
|---|---|---|
| 微调模型 | 50 | 0.408 |
| **官方权重 + 引导束搜索** | **63** | **0.597** |

结合官方权重训练数据（CrossDocked2020）的授权与收益比，**决定不把微调模型投入生产**，
采样仍使用官方预训练权重；该分支与数据全部保留在仓库中可复现。

### 4.5 工程可靠性（同样是改进）

- **失败要吵**：所有"静默丢弃"路径改为显式告警或计数（修复了三处对接管线静默丢样、
  配置文件末尾追加导致 YAML 非法等缺陷）；
- **可追溯**：每个对外数字都由脚本产出并落盘 CSV；
- **标注优先于删除**：对可疑分子打标签（`admet_alerts` / `可测性` / `tier`），而不是从库里删掉；
- **端到端自检**：`scripts/e2e_full_run.py` 覆盖 采样 → 入库 → 对接 → 训练 → 下游分析 八个阶段。

---

## 5. 端到端流程

```
靶点 PDB（4EIY / 2RH1 / 3PBL / 4IB4）
   │ ① 口袋提取：center + bbox_size(23.0 Å) 切出蛋白原子
   ▼
口袋张量（坐标 / 元素 / 残基类型 / kNN 图，cutoff 10 Å）
   │ ② 初始化采样：get_init() 采样首批原子（focal/pos/element 阈值过滤）
   ▼
束搜索池 pool{queue, finished, failed, duplicate}
   │ ③ 自回归生长（每步为每个候选预测下一个原子）
   │    sample_focal → sample_position(MDN) → sample_element_and_bond
   │ ④ 束排序 ★ 引导：prob ∝ (exp(Σlogp)+1)·weight·exp(λ·chem)
   ▼
完成分子 → 重建 RDKit mol
   │ ⑤ 可选构象精修（ETKDG + MMFF94s + 姿态保真）
   ▼
SMILES + 3D SDF
   │ ⑥ 入库漏斗：去重 → PAINS → BRENK → MW → LogP → QED → SA → 元素 → 环数
   ▼
化合物库（compounds.csv + library.db + 3D SDF，pdbqt_ready）
   │ ⑦ 下游分析（全脚本化）
   ▼
对接 → 新颖性 → 逆向对照 → ROC 富集 → 姿态置信 → ADMET → Pareto → 选择性 → 六维推荐
```

---

## 6. 代码地图

| 目录 | 规模 | 职责 |
|---|---|---|
| `models/` | 11 文件 / ~1,485 行 | 主干 `maskfill.py`、位置头 `position.py`、采样 `sample.py`、等变层 `invariant.py` |
| `utils/` | 9 文件 / ~1,686 行 | 特征化、训练工具、**引导算法 `guidance.py`**、**构象精修 `reconstruct.py`** |
| `evaluation/` | 5 文件 / ~556 行 | 打分函数、SA 评分、Vina 对接封装、批量评估 |
| `scripts/` | 38 文件 / ~3,846 行 | 数据构建、入库、对接、11+ 项分析、验证、批处理 |
| `gui/` | 8 文件 / ~3,086 行 | 桌面应用（`app_pro.py` + 组件库 `ui_kit.py` + 图标 + 打包脚本） |
| `configs/` | 37 文件 / ~875 行 | 训练与采样配置（含 A/B 矩阵、夜间批跑配置） |
| 顶层入口 | — | `sample.py`（测试集口袋）、`sample_for_pdb.py`（自定义 PDB 口袋）、`train.py` |

---

## 7. 复现步骤

### 7.1 环境

```bash
conda create -n pocket2mol python=3.10 -y
conda activate pocket2mol
# GPU 版 PyTorch（按 CUDA 版本调整）
pip install torch --index-url https://download.pytorch.org/whl/cu124
pip install torch-geometric torch-scatter torch-cluster torch-sparse
pip install rdkit pyyaml easydict lmdb biopython python-pptx
```

实测可用组合：Python 3.10.20 / torch 2.6.0+cu124 / RDKit 2026.03.5 / PyG 2.8.0 / Tk 8.6。

### 7.2 引导采样（λ=3 配置）

```bash
python sample_for_pdb.py --pdb_path ./targets/4EIY_A2A受体.pdb \
  --center " -0.4,8.5,17.1" \
  --config ./configs/sample_for_pdb_guided_l3.yml --outdir ./outputs/my_run
```

### 7.3 过滤入库

```bash
python scripts/build_library.py --runs ./outputs/my_run --library ./outputs/my_library
# 产出 compounds.csv / library.db / sdf/ / rejected.csv（带拒收原因） / summary.txt
```

### 7.4 对接与评估

```bash
python scripts/docking_pipeline.py --receptor ./targets/4EIY_A2A受体.pdb \
  --center=-0.4,8.5,17.1 --smiles-file ./outputs/my_library/smiles.txt --out ./outputs/docking
python scripts/roc_decoys.py            # 富集验证
python scripts/admet_screen.py          # ADMET 规则分级
python scripts/pareto_analysis.py       # 多目标前沿
python scripts/selectivity_analysis.py  # 选择性谱
python scripts/final_recommendation.py  # 六维综合推荐
```

### 7.5 端到端自检

```bash
python scripts/e2e_full_run.py          # 采样→入库→对接→训练→下游分析 八阶段
python scripts/verify_all_fixes.py      # 36 项工程校验
```

---

## 8. 数据资产

| 产物 | 说明 |
|---|---|
| `outputs/<批次>/` | 采样会话（`SMILES.txt`、`SDF/`、参数快照、日志） |
| `化合物库/compounds.csv` | 全库化合物与类药指标（QED/SA/MW/LogP/TPSA/HBD/HBA/骨架） |
| `化合物库/compounds_tagged.csv` | 带靶点标签与筛选标记 |
| `化合物库/top_candidates.csv` | 综合排序 Top 候选 |
| `化合物库/top50_novelty.csv` | Top50 与已知集对照（相似度/分类） |
| `化合物库/selectivity_full.csv` | 全库选择性谱 |
| `化合物库/admet_screen.csv` | ADMET 规则分级明细 |
| `化合物库/pareto_front.csv` | 五目标 Pareto 前沿 |
| `化合物库/pose_consistency.csv` | 对接姿态一致性 |
| `化合物库/final_recommendation.csv` | 六维 Tier 推荐 |
| `docs/7-eonmol_技术全书.md` | **项目技术全书**（12 章：模型解构 / 改进 / 评估体系 / 迭代史 / 缺陷审计） |

> 说明：大文件（预训练权重、数据集、输出目录）走 `.gitignore` 不进 git 历史，按需放置在对应目录。

---

## 9. 常见问题

| 现象 | 原因 / 处理 |
|---|---|
| `--center` 报参数错误 | 首个数值为负时**前导空格不能省**：`--center " -0.4,8.5,17.1"` |
| 采样很慢 | 建议 `taskset -c 0`（Linux）绑定单核，实测更快；或减小 `beam_size` |
| 结果列表为空 | 分子未生长完成（`max_steps` 太小）或未通过过滤漏斗；勾选 GUI「显示运行日志」查看后端输出 |
| 找不到靶点文件 | 「环境设置」页确认靶点 PDB 目录；部署模式下程序会自动探测 `structures/` 或 `targets/` |
| 引导开关无效 | 确认配置里 `sample.guided.enabled: true` 且 λ 生效（早期版本界面的 λ 未写入配置，已修复） |
| GUI 启动慢（>10 秒） | 旧版环境自检同步阻塞；现版本已异步化，冷启动约 2.6 秒、便携版 1.3 秒 |
| YAML 解析报 `ScannerError` | 不要手工在配置文件末尾追加缩进行；使用 `scripts/gen_sample_config.py` 生成配置 |
| 中文路径下 SDF 读取失败 | 早期 RDKit C++ 文件 API 的已知问题，已改为 `open()` + `MolFromMolBlock` 路径 |

---

## 10. 局限与未采纳路线

- **微调未晋级**：受限于训练数据规模与授权（CrossDocked2020 上游 PDBbind 对商用有限制），
  未能在本周期内产出优于官版权重的微调模型。
- **对接精度**：Vina 打分与真实亲和力相关性有限，故补充了 ROC 富集、姿态一致性、逆向对照三层交叉验证，
  但仍需湿实验确认。
- **选择性评估**：仅覆盖 A2A/A1 一对反向选择性，尚未做全家族选择性矩阵。
- **尚未实现**：动力学/MM-GBSA 精细重打分、主动学习闭环（需湿实验反馈）。
- **离线部分**：化合物送验、实验验证、答辩材料录制需线下完成。

---

## 11. 许可与引用

**许可**：本仓库基于 [Pocket2Mol](https://github.com/pengxingang/Pocket2Mol)（**MIT License**，
Copyright (c) 2022 Xingang Peng）构建，原许可与版权声明见 [`LICENSE`](./LICENSE)；
本项目新增的改进与脚本（引导束搜索、入库/对接/分析管线、桌面应用等）著作权归本项目作者。
第三方组件许可清单、数据来源与商用注意事项见 [`NOTICE.md`](./NOTICE.md)。

> 分发本软件时**必须保留** `LICENSE`（MIT 原文）——这是 MIT 的授权条件。

**引用**：使用本项目的改进方法（引导束搜索及相关流水线）时，请同时引用基础模型论文：

```bibtex
@inproceedings{peng2022pocket2mol,
  title     = {Pocket2Mol: Efficient Molecular Sampling Based on 3D Protein Pockets},
  author    = {Xingang Peng and Shitong Luo and Jiaqi Guan and Qi Xie and Jian Peng and Jianzhu Ma},
  booktitle = {International Conference on Machine Learning},
  year      = {2022}
}
```

---

## 12. 联系方式

- 邮箱：**2025158065@hrbmu.edu.cn**
- 问题反馈：欢迎提交 Issue；涉及数据或商用授权请邮件说明用途。

---

<sub>本 README 中的全部数字均由仓库内脚本对实际产物统计得出，可逐文件追溯；
更完整的技术细节见 `docs/7-eonmol_技术全书.md`。</sub>
