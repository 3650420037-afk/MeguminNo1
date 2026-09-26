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

## 提交信息速览

| 项 | 内容 |
|---|---|
| 赛道 | 赛道 3 · 离子通道/GPCR 小分子药物虚拟筛选 · 子任务 2 |
| 主运行入口 | `python predict.py --target A2A` → 生成 `results/results.csv` |
| 设计/生成入口 | `python design.py --target A2A` |
| 训练入口 | `python train.py --config configs/train_multitarget_v2.yml` |
| 可执行 Notebook | `notebooks/7-eonmol_流程演示.ipynb` |
| 模型 | Pocket2Mol 官方预训练权重（**未微调**）+ 本项目引导束搜索；见 [`ModelCard.md`](./ModelCard.md) |
| 数据 | RCSB PDB（CC0）+ ChEMBL（CC BY-SA 3.0）；来源/许可/划分见 [`data/README.md`](./data/README.md) |
| 许可 | Pocket2Mol 为 MIT；第三方归属见 [`NOTICE.md`](./NOTICE.md) |
| 联系方式 | 2025158065@hrbmu.edu.cn |

---

## 目录结构

按大赛《附件5·代码提交要求》的模板组织：

```
.
├── README.md              项目说明、环境配置、运行命令、输入输出与结果说明
├── requirements.txt       依赖库及锁定版本（环境/硬件说明见 README 第三节）
├── ModelCard.md           模型说明：适用范围、输入输出、已知局限、创新贡献
├── data/                  数据说明与获取方式（来源/许可/清洗/划分/泄漏防控）
├── src/                   核心源代码
│   ├── models/            等变 MPNN 主干、位置/元素/键预测头
│   ├── utils/             特征化、训练工具、引导算法、构象精修
│   ├── scripts/           数据构建、入库、对接、分析、验证脚本
│   ├── evaluation/        打分函数、SA 评分、对接封装
│   ├── gui/               桌面应用（可选的人机界面）
│   ├── paths.py           全项目唯一的路径来源（禁止硬编码绝对路径）
│   └── sample.py / sample_for_pdb.py   底层采样实现
├── models/                最终模型权重（pretrained_Pocket2Mol.pt）
├── notebooks/             可执行 Notebook（流程演示与结果核验）
├── train.py               一键训练脚本
├── design.py              一键设计/生成脚本
├── predict.py             一键生成最终结果文件（候选清单）
├── configs/               训练与采样配置（targets.json 为靶点登记表）
├── results/               运行结果与候选清单（results.csv / structures/ / run_manifest.json）
├── logs/                  训练与采样日志、随机种子记录
├── docs/                  技术全书、验证报告等文档
├── tools/                 第三方对接工具（AutoDock Vina / Open Babel，可不安装）
└── other/                 非提交内容（历史产物、备份、原始下载），不进 git
```

> 说明：`models/` 是**权重**，`src/models/` 是**模型架构代码**，两者按模板要求分开。
> `other/`、`tools/`、大数据集与运行产物均由 `.gitignore` 排除，说明见
> [`other/README.md`](./other/README.md)。

---

## 目录

- [一页速览](#一页速览)
- [1. 任务与定位](#1-任务与定位)
- [2. 主要成果](#2-主要成果)
- [3. 环境与依赖](#3-环境与依赖)
- [4. 快速开始](#4-快速开始)
- [5. 最终结果文件](#5-最终结果文件)
- [6. 端到端流程](#6-端到端流程)
- [7. 代码地图](#7-代码地图)
- [8. 复现步骤](#8-复现步骤)
- [9. 数据与模型资产](#9-数据与模型资产)
- [10. 常见问题](#10-常见问题)
- [11. 局限与未采纳路线](#11-局限与未采纳路线)
- [12. 许可与引用](#12-许可与引用)
- [13. 联系方式](#13-联系方式)

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

内置五个靶点（蛋白结构随仓库提供，登记于 `configs/targets.json`）：

| 靶点 | PDB | 口袋中心 (x,y,z) | 适应症 |
|---|---|---|---|
| A2A 腺苷受体 | 4EIY | −0.42, 8.53, 17.13 | 帕金森病 / 肿瘤免疫 / 炎症 |
| β2 肾上腺素受体 | 2RH1 | −29.52, 9.23, 6.94 | 哮喘 / 心血管疾病 |
| D3 多巴胺受体 | 3PBL | 0.09, −14.83, 10.43 | 精神分裂症 / 成瘾 |
| 5-HT2B 血清素受体 | 4IB4 | 22.45, 18.28, 11.73 | 偏头痛 / 肺动脉高压 |
| A1 腺苷受体 | 5UEN | 55.97, 58.90, 143.62 | 心血管 / 神经保护 |

> 口袋中心由 `src/scripts/build_target_registry.py` 从结构文件**自动提取**
> （= 共晶配体重原子的实验质心），不手工录入，避免写错；同时自动披露
> "沉积配体是否完整"。详见 `data/README.md`。

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

## 3. 环境与依赖

### 3.1 完整依赖清单

锁定版本见 [`requirements.txt`](./requirements.txt)（2026-09 实测环境）：

| 类别 | 包与版本 |
|---|---|
| 解释器 / 系统 | **Python 3.10.20** · Windows 11（Linux 同理） |
| 深度学习 | torch **2.6.0+cu124** · torch-geometric **2.8.0** · torch-scatter **2.1.2** · torch-cluster **1.6.3** · torch-sparse **0.6.18** |
| 化学信息学 | rdkit **2026.3.5** · biopython **1.88** |
| 数值/数据 | numpy **2.2.6** · scipy **1.15.3** · lmdb **2.3.0** · networkx **3.4.2** |
| 配置/工具 | easydict **1.13** · PyYAML **6.0.3** · tqdm **4.70.0** · tensorboard **2.21.0** |
| 报告生成 | python-pptx **1.0.2** |
| 第三方可执行（可选） | AutoDock Vina 1.2.5 · Open Babel 3.2.1（位于 `tools/`，分子对接用） |

### 3.2 GPU / 驱动 / 硬件

| 项 | 实测值 |
|---|---|
| GPU | NVIDIA GeForce RTX 4070 Laptop GPU，**8 GB 显存** |
| 驱动 | 610.62 |
| CUDA | 12.4（torch 自带运行时） |
| 纯 CPU 可用性 | 可运行（`--device cpu`），推理显著变慢 |

### 3.3 安装

```bash
# 1) 创建环境
conda create -n pocket2mol python=3.10 -y
conda activate pocket2mol

# 2) 安装依赖（torch 与 PyG 扩展需匹配本机 CUDA）
pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
pip install torch-scatter==2.1.2 torch-cluster==1.6.3 torch-sparse==0.6.18 \
    -f https://data.pyg.org/whl/torch-2.6.0+cu124.html
pip install -r requirements.txt

# 3) 环境自检
python -c "import torch,rdkit,torch_geometric;print(torch.__version__, torch.cuda.is_available())"
```

> ⚠️ 本工程含 **Windows / PyTorch 2.6 兼容补丁**（`src/utils/misc.py`、`src/utils/transforms.py`、
> `src/models/common.py`），请勿把 torch 降到 2.6 以下，也不要删除这些补丁。
> `env_cuda113.yml` 是**原版 Pocket2Mol 的 Linux 环境文件**（python 3.8 / torch 1.10.1+cu113），
> 与本项目实测环境不一致，仅供溯源，请勿用它复现环境。

### 3.4 预期运行时间与资源| 任务 | 配置 | 实测耗时（RTX 4070 Laptop 8 GB） |
|---|---|---|
| 快速自检生成 | 10 分子 / 束宽 50 / 20 步 | 2–4 分钟 |
| 标准生成 | 50 分子 / 束宽 100 / 50 步 | ≈ 10 分钟 |
| 一键出清单（`predict.py`，含入库与 3D 构象生成） | 同上 | ≈ 15–25 分钟 |
| 训练冒烟 | `max_iters: 2` | 2–5 分钟 |
| 完整训练/微调 | `configs/train_multitarget_v2.yml` | 数小时（视迭代数） |
| 分子对接（可选） | 每分子 1 次 Vina | 约 2–10 秒/分子，需已放置 `tools/` |

显存占用：标准生成约 2–4 GB；训练约 6–8 GB。
**提示**：生成耗时由**束宽**主导，其次是分子长度（单步耗时随队列原子数增长）；
`max_steps` 是**安全上限**而非目标步数，设得过小会提前结束导致完成分子偏少。

### 3.5 解释器与路径的环境变量

Python 入口（`design.py` / `predict.py` / `train.py` / `src/scripts/*.py`）会自动使用
**当前解释器**（`sys.executable`）并把仓库的 `src/` 加入导入路径，
因此只要用同一个环境执行即可，无需任何环境变量。

只有 **PowerShell / 批处理脚本**（`src/scripts/*.ps1`、`src/gui/build_exe.ps1`、
`src/scripts/run_sample.bat`）没有 `sys.executable` 可用，它们按以下顺序取解释器：

1. 环境变量 `EONMOL_PYTHON`（推荐显式设置，避免机器上有多个 Python 时取错）
2. 否则回退到 `PATH` 上的 `python`

```powershell
# 建议在会话开始时设置一次(WSL/CI 同样适用)
$env:EONMOL_PYTHON = 'D:\Miniconda3\envs\Pocket2Mol\python.exe'
```

同理，若把大数据集或第三方工具放在仓库之外，可用这些环境变量覆盖默认位置
（默认值全部由 `src/paths.py` 从仓库根推导）：

| 变量 | 覆盖的默认位置 |
|---|---|
| `EONMOL_DATA` | `data/` |
| `EONMOL_MODELS` | `models/` |
| `EONMOL_TOOLS` | `tools/` |
| `EONMOL_VINA` / `EONMOL_OBABEL` / `EONMOL_OBABEL_DATA` | `tools/` 下的对接工具 |
| `EONMOL_PYTHON` | 子进程使用的解释器 |

---

## 4. 快速开始

### 方式一：一键生成候选清单（推荐，评测最小闭环）

```bash
# 完整流程：生成 → 过滤入库 → 打分排序 → results/results.csv
python predict.py --target A2A

# 小规模自检（几分钟）
python predict.py --target A2A --num-samples 10 --beam 50 --max-steps 20

# 复用已生成的会话，只重跑筛选与排序
python predict.py --skip-design --runs outputs/design/A2A_20260101_120000 --target A2A
```

### 方式二：只做生成

```bash
python design.py --target A2A --num-samples 50 --beam 100 --seed 2024
# 关闭引导（复现原版基线用于 A/B 对照）
python design.py --target A2A --lam 0
```

### 方式三：Notebook（逐步演示与结果核验）

```bash
jupyter notebook notebooks/7-eonmol_流程演示.ipynb
```

### 方式四：桌面应用（给非编程用户）

```bash
python src/gui/app_pro.py            # 可选 --tab 0-3 --theme dark
```

便携版：解压 `7-eonmol_便携版.zip` → 双击 `7-eonmol.exe`（首次运行若提示
"Windows 已保护你的电脑"，点「更多信息 → 仍要运行」）。
自包含部署包：整个文件夹拷到目标机 → 双击 `一键配置.bat`（只需一次）→ 双击 `7-eonmol.exe`。

### 方式五：训练 / 微调

```bash
python train.py --config configs/train_multitarget_v2.yml
# 训练冒烟（2 次迭代，用于验证环境与数据链路）
python train.py --config configs/train_multitarget_v2.yml --logdir logs/smoke
```

---

## 5. 最终结果文件

`predict.py` 一键产出：

| 文件 | 说明 |
|---|---|
| `results/results.csv` | **最终候选清单**（UTF-8），字段见下表 |
| `results/structures/*.sdf` | 候选的三维结构，与清单「结构文件」列一一对应 |
| `results/run_manifest.json` | 本次运行的参数、随机种子、权重路径、环境、耗时（可溯源） |
| `<outdir>/library/compounds.csv` | 入库全量分子及描述符（未截断，供复核） |

### 清单字段

| 字段 | 含义 |
|---|---|
| `候选编号` | 形如 `C0001` |
| `所属赛道` | 赛道 3 · 离子通道/GPCR 小分子药物虚拟筛选 · 子任务 2 |
| `候选SMILES` | 候选结构（规范化 SMILES） |
| `结构文件` | 对应 SDF 的相对路径 |
| `关键预测指标_类药性综合分` | `QED + (1 − SA/10)`，即生成阶段的引导目标 |
| `关键预测指标_QED` / `_SA` | 类药性 / 合成可及性（1–10，越低越易合成） |
| `关键预测指标_新颖性最大Tanimoto` | 与靶点已知活性分子的最大相似度（越低越新） |
| `关键预测指标_对接分Vina(kcal/mol)` | AutoDock Vina 对接分（越负结合越强）；由 `src/scripts/docking_pipeline.py` 产出 |
| `关键预测指标_ADMET分级` | 规则化 ADMET 分级 A/B/C；`ADMET警报` 为命中的规则名 |
| `关键预测指标_选择性SI百分位` | A2A vs A1 反向选择性百分位；`可测性` 为可测性评估 |
| `关键预测指标_六维综合Tier` | 项目六维加权综合推荐的 Tier（1 最优） |
| `MW` `LogP` `TPSA` `LogS` `HBD` `HBA` `可旋转键` `环数` `Lipinski违反数` | 理化性质 |
| `姿态能量差(kcal/mol)` | 最优/次优构象能量差，越小姿态越可信 |
| `Murcko骨架` / `同骨架候选数` | 骨架与生成稳健性提示 |
| `排序` / `靶点` | 排序位次 / 靶点名 |
| `模型名称` `模型版本` `代码版本` | 例如 Pocket2Mol 官方权重（未微调）+ 引导束搜索 |
| `随机种子` `运行编号` | 复现所需 |
| `备注` | 自动标注：3D 构象未生成、与已知活性高度相似（Tanimoto≥0.8）、Lipinski 违反 >1、ADMET 分级 C、无对接分 |

### 排序逻辑与不确定性

`--rank-by` 控制排序依据，默认 `auto`：

| 取值 | 排序规则 |
|---|---|
| `auto`（默认） | 若库内有六维综合 Tier → 按 Tier（Tier1→Tier2→…，组内按类药性综合分）；否则若库内有对接分 → 按对接分升序（越负越强，缺失者置底）；再否则按类药性综合分降序 |
| `tier` / `vina` / `druglikeness` | 强制使用对应依据（数据缺失时自动退回类药性综合分并打印告警） |

- `类药性综合分 = QED + (1 − SA/10)`，与生成阶段的引导打分**同式**，
  保证"用来搜索的目标"与"用来排序的指标"一致。
- 清单会**自动并入库目录下已有的分析结果**（`admet_screen.csv`、`selectivity_full.csv`、
  `pose_consistency.csv`、`final_recommendation.csv`、`top_candidates.csv`），
  以规范化 SMILES 为键；文件缺失时对应列留空，不影响清单生成。
  用到了哪些文件会写进 `results/run_manifest.json` 的 `enrichment` 字段。
- **不确定性**：本模型为单点生成，**不输出置信区间**；因此清单给出
  「同骨架候选数」作为生成稳健性的弱提示、「姿态能量差」作为对接姿态可信度提示，
  并在 `备注` 中标注高相似、构象缺失与 ADMET 分级 C 的候选。
  清单中的全部指标均为**计算代理指标**，未经湿实验验证。

> 模型本身**不预测**结合亲和力与活性。清单中的对接分来自 **AutoDock Vina** 这一独立
> 第三方工具，类药性/新颖性/ADMET 指标来自 RDKit 与规则化评估，均与模型预测是两回事，
> README 与清单中分别标注。

---

## 6. 端到端流程

```
靶点 PDB（data/targets/，5 个 GPCR）
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
   │ ⑦ 打分排序 → results/results.csv（predict.py 自动完成）
   ▼
（可选）下游分析：对接 → 新颖性 → 逆向对照 → ROC 富集 → 姿态置信
                → ADMET → Pareto → 选择性 → 六维推荐
```

---

## 7. 代码地图

| 目录 | 规模 | 职责 |
|---|---|---|
| `src/models/` | 11 文件 / ~1,485 行 | 主干 `maskfill.py`、位置头 `position.py`、采样 `sample.py`、等变层 `invariant.py` |
| `src/utils/` | 9 文件 / ~1,686 行 | 特征化、训练工具、**引导算法 `guidance.py`**、**构象精修 `reconstruct.py`** |
| `src/evaluation/` | 5 文件 / ~556 行 | 打分函数、SA 评分、Vina 对接封装、批量评估 |
| `src/scripts/` | 40+ 文件 / ~4,000 行 | 数据构建、入库、对接、11+ 项分析、验证、批处理 |
| `src/gui/` | 8 文件 / ~3,086 行 | 桌面应用（`app_pro.py` + 组件库 `ui_kit.py` + 图标 + 打包脚本） |
| `src/paths.py` | — | **全项目唯一路径来源**，禁止在其它模块硬编码盘符路径 |
| `configs/` | 38 文件 | 训练与采样配置（含 A/B 矩阵、靶点登记表 `targets.json`） |
| 顶层入口 | — | `design.py`（一键生成）、`predict.py`（一键出清单）、`train.py`（一键训练） |
| 底层采样 | — | `src/sample.py`（测试集口袋）、`src/sample_for_pdb.py`（自定义 PDB 口袋） |

---

## 8. 复现步骤

```bash
# 0) 环境：见第 3 节
conda activate pocket2mol

# 1) 一键生成最终候选清单（评测主入口）
python predict.py --target A2A --num-samples 50 --beam 100 --seed 2024 --top 100
#    -> results/results.csv / results/structures/ / results/run_manifest.json

# 2) 生成（可选，独立于清单）
python design.py --target A2A --num-samples 50 --beam 100

# 3) 过滤入库
python src/scripts/build_library.py --runs outputs/design/A2A_<时间戳> \
    --library outputs/my_library
#    产出 compounds.csv / library.db / sdf/ / rejected.csv（带拒收原因） / summary.txt

# 4) 对接与下游分析（可选，需要 tools/ 下的 Vina 与 Open Babel）
python src/scripts/docking_pipeline.py \
    --receptor data/targets/4EIY_A2A受体.pdb \
    --center=-0.42,8.53,17.13 \
    --smiles-file outputs/my_library/smiles.txt --out outputs/docking
python src/scripts/roc_decoys.py            # 富集验证
python src/scripts/admet_screen.py          # ADMET 规则分级
python src/scripts/pareto_analysis.py       # 多目标前沿
python src/scripts/selectivity_analysis.py  # 选择性谱
python src/scripts/final_recommendation.py  # 六维综合推荐

# 5) 端到端自检
python src/scripts/e2e_full_run.py --fast   # 采样→入库→对接→训练→下游分析 八阶段
python src/scripts/verify_all_fixes.py      # 36 项工程校验

# 6) 靶点登记表重建（结构变动后）
python src/scripts/build_target_registry.py
```

**关于 `--center`**：`src/sample_for_pdb.py` 直接调用时，首个数值为负时**前导空格不能省**
（负数参数解析约定），例如 `--center " -0.4,8.5,17.1"`；用 `design.py` / `predict.py`
时由程序生成，无需手工处理。

---

## 9. 数据与模型资产

| 资产 | 位置 | 说明 |
|---|---|---|
| 数据集来源/许可/清洗/划分/泄漏防控 | `data/README.md` | **必读**；含各数据集样本数与划分记录 |
| 靶点结构 | `data/targets/` | 5 个 GPCR 受体坐标 + 共晶配体定义（RCSB PDB，CC0） |
| 最小示例 | `data/example/` | `4yhj.pdb` + 参考配体 SDF |
| 已知活性/诱饵参考集 | `data/known_drugs/` | 新颖性与富集评估用（ChEMBL，CC BY-SA） |
| 模型权重 | `models/pretrained_Pocket2Mol.pt` | 官方预训练权重（MIT）；下载说明见 `models/README.md` |
| 模型卡 | `ModelCard.md` | 适用范围、输入输出、已知局限、创新贡献 |
| 候选清单 | `results/results.csv` 等 | 见第 5 节 |
| 全库化合物 | `results/compounds.csv`、`results/*.csv` | QED/SA/MW/骨架/对接/ADMET/Pareto/选择性/推荐 |
| 三维结构 | `results/sdf/`、`results/structures/` | 每个分子一个 SDF |
| 采样会话 | `outputs/<批次>/` | `SMILES.txt`、`SDF/`、参数快照、日志 |
| 训练与采样日志 | `logs/` | 含随机种子与检查点 |
| 技术文档 | `docs/` | 技术全书、数据建库穿模修复验证报告、参赛代码说明 |

> **不进 git 的内容**：预训练数据的 GB 级 LMDB、历代采样产物、旧备份。
> 获取与重建方式见 `data/README.md` 第六节；非提交内容的清单见 `other/README.md`。

---

## 10. 常见问题

| 现象 | 原因 / 处理 |
|---|---|
| `ModuleNotFoundError: models / utils` | 必须从**仓库根**运行入口脚本（`design.py` / `predict.py` / `train.py`），它们会把 `src/` 加入导入路径 |
| 找不到权重 | 确认 `models/pretrained_Pocket2Mol.pt` 存在；下载说明见 `models/README.md` |
| `--center` 报参数错误 | 直接调用 `src/sample_for_pdb.py` 时，首个数值为负需保留前导空格：`--center " -0.4,8.5,17.1"`；改用 `design.py` 可避免 |
| 采样很慢 | 减小 `beam_size`（主要成本旋钮）；Linux 下可 `taskset -c 0` 绑核 |
| 结果列表为空 | 分子未生长完成（`max_steps` 太小）或未通过过滤漏斗；查看 `outputs/.../*.log` 与 `library/rejected.csv` |
| 引导开关无效 | 确认配置里 `sample.guided.enabled: true` 且 `λ` 生效；用 `design.py --lam 3` 会写入配置 |
| GUI 启动慢（>10 秒） | 旧版环境自检同步阻塞；现版本已异步化，冷启动约 2.6 秒、便携版 1.3 秒 |
| YAML 解析报 `ScannerError` | 不要手工在配置末尾追加缩进行；用 `src/scripts/gen_sample_config.py` 生成配置 |
| 中文路径下 SDF 读取失败 | 早期 RDKit C++ 文件 API 的已知问题，已改为 `open()` + `MolFromMolBlock` 路径 |
| 对接步骤跳过 | `tools/` 下缺少 Vina/Open Babel；不使用对接时属预期行为，清单中不含对接分数 |

---

## 11. 局限与未采纳路线

- **权重未针对本项目靶点微调**：两轮微调（冻结底层 + 微调头部，含验证集早停）**均未通过
  A/B 对照**——类药性中位 基线 0.780 vs 微调 v1 0.542 vs v2 0.578；QED≥0.8 占比
  基线 37.9% vs 1.1% vs 2.3%。按"不达标即如实记录并保留基线"的原则**整体回退**，
  交付仍使用官方预训练权重。失败原因分析见 `docs/7-eonmol_技术全书.md` 与 `ModelCard.md`。
- **元素词表受限**：含 Br/I/Na 等词表外元素的分子无法处理（实测占候选集 4%–19%）。
- **引导会牺牲多样性**：λ 越大 QED 越高但结构多样性下降，`w_div` 只能部分补偿。
- **对接精度**：Vina 打分与真实亲和力相关性有限，故补充了 ROC 富集、姿态一致性、逆向对照三层交叉验证，
  但仍需湿实验确认。
- **选择性评估**：仅覆盖 A2A/A1 一对反向选择性，尚未做全家族选择性矩阵。
- **无不确定性量化**：单点预测，无置信区间。
- **尚未实现**：动力学/MM-GBSA 精细重打分、主动学习闭环（需湿实验反馈）。
- **离线部分**：化合物送验、实验验证、答辩材料录制需线下完成。

---

## 12. 许可与引用

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

## 13. 联系方式

- 邮箱：**2025158065@hrbmu.edu.cn**
- 问题反馈：欢迎提交 Issue；涉及数据或商用授权请邮件说明用途。

---

<sub>本 README 中的全部数字均由仓库内脚本对实际产物统计得出，可逐文件追溯；
更完整的技术细节见 `docs/7-eonmol_技术全书.md`，模型说明见 `ModelCard.md`，
数据说明见 `data/README.md`。</sub>
