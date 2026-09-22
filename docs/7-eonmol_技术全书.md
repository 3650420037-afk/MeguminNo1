# Pocket2Mol-GPCR 项目技术全书

> 版本 1.0 · 2026-09-22 · 面向团队成员、答辩评审与后续接手者
> 仓库：`D:\MMModel\Pocket2Mol`（git 远端 `3650420037-afk/MeguminNo1`，73 次提交）
> 本文所有数字均由仓库内脚本对实际产物统计得出，可逐文件追溯。

---

## 目录

- [第 0 章 一页速览](#第-0-章-一页速览)
- [第 1 章 任务与问题定义](#第-1-章-任务与问题定义)
- [第 2 章 系统总览](#第-2-章-系统总览)
- [第 3 章 基础模型解构（Pocket2Mol）](#第-3-章-基础模型解构pocket2mol)
- [第 4 章 我们的改进（核心贡献）](#第-4-章-我们的改进核心贡献)
- [第 5 章 推理全流程](#第-5-章-推理全流程从口袋到候选库)
- [第 6 章 评估与筛选体系](#第-6-章-评估与筛选体系)
- [第 7 章 数据资产与结果汇总](#第-7-章-数据资产与结果汇总)
- [第 8 章 工程实现](#第-8-章-工程实现)
- [第 9 章 迭代进程史](#第-9-章-迭代进程史)
- [第 10 章 缺陷审计与修复方法论](#第-10-章-缺陷审计与修复方法论)
- [第 11 章 局限、未采纳路线与展望](#第-11-章-局限未采纳路线与展望)
- [第 12 章 附录](#第-12-章-附录)

---

## 第 0 章 一页速览

**一句话**：基于 Pocket2Mol（口袋条件自回归三维分子生成模型），用"引导束搜索"把生成分布推向类药区，
构建了四靶点共 6,864 个类药化合物库，并完成对接、新颖性、富集、ADMET、多目标、选择性六层计算评估。

| 维度 | 结论 |
|---|---|
| 生成模型 | Pocket2Mol 官方预训练权重（42.8 MB），等变消息传递网络（MPNN）+ 混合密度网络（MDN）位置头 |
| 我们的核心改进 | **引导束搜索**：束排序概率乘化学分数增益 `exp(λ·chem)`，λ=3 时 QED 0.667→0.803（+20%） |
| 伴随改进 | Tanimoto 多样性惩罚、生成构象力场精修、入库漏斗、全自动对接与分析管线 |
| 主库规模 | A2A 2,931 分子 / 1,960 唯一 Murcko 骨架 / QED 中位 0.825 / SA 中位 3.33 / MW 中位 362 |
| 四靶点合计 | 原始 10,004 → 过滤入库 6,864 |
| 对接（A2A） | 2,934 分子全部成功；均值 −9.68 kcal/mol；≤−10 有 1,137 个（38.8%）；最强 −13.51 |
| 逆向阳性对照 | 737 个已知 A2A 活性分子同管线均值 −9.32，生成库更优 |
| 富集验证 | ROC-AUC 0.685；≤−12 富集 5.80×；top-1% 富集 2.24× |
| Top50 候选 | Vina −13.42 ~ −11.59；对已知集最大相似度仅 0.187–0.321；全部 HIGH_NOVEL |
| 成药性 | ADMET 规则分级 A 676 / B 1,940 / C 315 |
| 多目标 | 五目标 Pareto 前沿 91 个非支配解（旧 Top50 仅 27 个在沿上） |
| 选择性 | A2A vs A1 反向选择性 SI 中位 +4.0，59% 分子 >0，382 个 >25 |
| 交付形态 | 桌面 GUI（免安装 exe，1.3–2.6 秒启动）+ 5 GB 自包含部署包 + 全部脚本与验证报告 |

---

## 第 1 章 任务与问题定义

### 1.1 赛事任务

**赛道 3 · 离子通道/GPCR 小分子药物虚拟筛选 · 子任务 2**：
给定 GPCR 靶点口袋结构，生成**结构新颖、具备类药性**的小分子候选，并证明流程可执行、结果可复现。

拆解成工程要求：
1. **口袋条件生成**：不是无条件生成，分子必须"长在口袋里"；
2. **类药性**：QED、合成可及性、理化性质须达标，而非仅满足几何；
3. **新颖性**：与已知配体骨架不能重合；
4. **可验证**：每一步产物、每个数字都要能用脚本重算。

### 1.2 技术路线的选择

| 路线 | 代表方法 | 取舍 |
|---|---|---|
| 基于配体 | 分子生成/优化（如 JT-VAE、REINVENT） | 需要大量已知活性配体，GPCR 多靶点数据稀疏 |
| 基于结构 | **口袋条件 3D 生成（Pocket2Mol）** | 直接利用晶体结构，无需配体数据；三维构象原生输出 |
| 基于对接筛选 | 虚拟筛选已有库 | 依赖库质量，无法创造新化学空间 |

**最终定位（贯穿全项目的一句话）**：**造引擎，不挖矿**。
我们不复现已有化合物库的筛选，而是把"从口袋结构生成可用候选"的能力做成一台可复用的引擎。

---

## 第 2 章 系统总览

### 2.1 端到端数据流

```
靶点 PDB（4EIY/2RH1/3PBL/4IB4）
   │  ① 口袋提取：以给定 center + bbox_size(23.0 Å) 切出蛋白原子
   ▼
口袋张量（原子坐标/元素/氨基酸类型/kNN 图）
   │  ② 初始化采样：get_init() 采样首批原子（focal/pos/element 阈值过滤）
   ▼
束搜索池 pool{queue, finished, failed, duplicate}
   │  ③ 自回归生长循环（每一步对每个候选预测下一个原子）
   │     ├─ sample_focal    ：在已有原子上选生长位点
   │     ├─ sample_position ：MDN 三分量预测新原子坐标
   │     └─ sample_element_and_bond ：预测元素与键级
   │  ④ 束排序（★ 我们的改进点）
   │     原版：prob ∝ (exp(Σlogp)+1)·weight
   │     引导：prob ∝ (exp(Σlogp)+1)·weight·exp(λ·[QED+(1−SA/10)−w·maxTanimoto])
   ▼
完成分子（STATUS_FINISHED）→ 原子/键重建 RDKit mol
   │  ⑤ 可选构象精修（ETKDG 内部几何 + MMFF94s + 姿态保真约束）
   ▼
SMILES + 3D SDF
   │  ⑥ 过滤入库漏斗：去重 → PAINS → BRENK → MW → LogP → QED → SA → 元素 → 环数
   ▼
化合物库（compounds.csv + library.db + 3D SDF，pdbqt_ready）
   │  ⑦ 下游分析（全部脚本化）
   ▼
对接 → 新颖性 → 逆向对照 → ROC 富集 → 姿态置信 → ADMET → Pareto → 选择性 → 可测性 → 自洽性 → 六维推荐
```

### 2.2 代码地图

| 目录 | 文件/行数 | 职责 |
|---|---|---|
| `models/` | 11 文件 / 1,485 行 | 模型定义：主干 `maskfill.py`、位置头 `position.py`、采样逻辑 `sample.py`、等变层 `invariant.py`、通用层 `common.py` |
| `utils/` | 9 文件 / 1,686 行 | 数据特征化、训练工具、**引导算法 `guidance.py`**、**构象精修 `reconstruct.py`**、蛋白配体解析 |
| `evaluation/` | 5 文件 / 556 行 | 打分函数、SA 评分 `sascorer.py`、对接封装 `docking.py`、批量评估 `evaluate.py` |
| `scripts/` | 38 文件 / 3,846 行 | 数据构建、入库、对接、11+ 项分析、验证、批处理 |
| `gui/` | 8 文件 / 3,086 行 | Tkinter 桌面应用（源码 + 组件库 + 图标 + 打包脚本） |
| `configs/` | 37 文件 / 875 行 | 训练与采样配置（含 A/B 实验矩阵、夜间批跑配置） |
| 顶层入口 | `sample.py` / `sample_for_pdb.py` / `train.py` | 官方采样、口袋文件采样（我们扩展）、训练 |

### 2.3 设计原则（贯穿全项目）

1. **零侵入**：所有改进以"新增模块 + 配置开关"实现，`config.sample.guided` 不存在时**行为与原版逐位一致**。
   这保证了任何改进都可以被独立关闭、独立证伪。
2. **可追溯**：每个对外数字都由脚本产出并落盘 CSV，报告中任意结论可回溯到文件行。
3. **标注优先于删除**：对可疑分子打标签（`admet_alerts`、`可测性`、`tier`），而不是直接从库里删掉。
4. **失败要吵**：凡是"静默丢弃"的代码路径都改造成显式告警或计数（后文第 10 章详述）。

---

## 第 3 章 基础模型解构（Pocket2Mol）

> 本章依据 `models/`、`utils/transforms.py`、`utils/protein_ligand.py` 与 `configs/train_multitarget_v2.yml` 撰写。

### 3.1 问题形式化

自回归三维分子生成：给定蛋白口袋上下文 $\mathcal{P}$，按"原子逐个生长"的方式生成分子
$M = \{(x_i, e_i)\}_{i=1}^{N}$ 与键集 $\mathcal{B}$：

$$
p(M \mid \mathcal{P}) = \prod_{t=1}^{N} \underbrace{p(f_t \mid S_{<t})}_{\text{生长位点}} \cdot
\underbrace{p(\Delta x_t \mid S_{<t}, f_t)}_{\text{新原子位置}} \cdot
\underbrace{p(e_t \mid S_{<t}, f_t, \Delta x_t)}_{\text{元素/键}}
$$

其中 $S_{<t}$ 是已生成的部分分子（部分图），$f_t$ 是被选中的"前沿原子"（focal）。
一步只加一个原子，因此**分子是一步步"长"出来的**，而不是一次性输出图。

### 3.2 特征化

| 对象 | 特征内容 |
|---|---|
| 蛋白原子 | 元素类型、氨基酸类型（残基种类）、是否主链、坐标（标量 + 方向向量） |
| 配体原子 | 元素类型、是否芳香、杂化、度数、形式电荷、手性、坐标 |
| 边 | 由 kNN 图构造：蛋白侧 `knn=48`、配体-蛋白界面 `knn=32`、截断 `cutoff=10.0 Å` |
| 几何 | 距离经 GaussianSmearing 展开；方向向量单独成通道（等变分支） |

关键实现：
- `utils/protein_ligand.py`：PDB/mmCIF 解析、配体提取、`parse_sdf_file`（含芳香键 warn 契约）；
- `utils/transforms.py`（760 行）：特征组合器 `AtomComposer`、掩码变换、对比采样、kNN 建图；
- `models/common.py`：`compose_context*`（蛋白+配体拼接）、`get_compose_knn_graph`、`embed_compose`。

### 3.3 主干网络（配置实测值）

```yaml
encoder:            # 名称 cftfm（CFTransformer 风格的消息传递网络）
  hidden_channels: 256          # 标量通道
  hidden_channels_vec: 64       # 方向（等变）通道
  edge_channels: 64
  key_channels: 128
  num_heads: 4                  # 注意力头
  num_interactions: 6           # 消息传递层数
  cutoff: 10.0                  # 相互作用截断（Å）
  knn: 48                       # 蛋白侧近邻数
```

**架构性质（需明确澄清）**：主干是 **等变消息传递网络（MPNN / 图神经网络）**，带 4 头注意力聚合；
**不是 Transformer 序列模型**。项目早期文档里"Transformer"的表述已在 hardening 阶段修正
（见提交 `47def1d`）。激活函数实际使用的是 **LeakyReLU**，而不是原论文仓库注释里写的 ShiftedSoftplus ——
这一条也经过代码级核对后写入了文档。

### 3.4 三个预测头

**(1) 生长位点头（focal）** —— `sample_focal()`
在已有原子的前沿集合中打分，输出"下一个原子挂在谁身上"，阈值 `focal_threshold=0.5`。

**(2) 位置头（位置预测）** —— `models/position.py::PositionPredictor`
用**混合密度网络（MDN）**建模连续三维坐标分布：

```yaml
position:
  num_filters: 128
  n_component: 3      # 3 个高斯分量
```

- 输出每个分量的 $\mu_k \in \mathbb{R}^3$、$\sigma_k$、混合权重 $\pi_k$；
- 训练时 `get_mdn_probability()` 计算目标位置在对数似然中的概率；
- 采样时 `sample_batch()` 按 $\pi$ 选分量、按 $\mathcal{N}(\mu_k,\sigma_k)$ 采坐标；
- 阈值 `pos_threshold=0.25` 过滤掉低概率位置。

**为什么用 MDN 而不是直接回归**：同一口袋位点存在多个合理结合姿态（多模态），
单峰回归会给出"平均姿态"这种物理上不存在的结果。

**(3) 元素与键头** —— `sample_element_and_bond()` / `sample_init_element()`
- 元素（atom type）分类，阈值 `element_threshold=0.3`；
- `hasatom` 二分类（是否还有原子可加），阈值 `hasatom_threshold=0.6` —— 这是**终止信号**；
- 键级分类（单/双/三键），阈值 `bond_threshold=0.4`。

### 3.5 训练目标

`MaskFillModelVN.get_loss()` 一次性返回 7 项损失，实测第 1 次迭代的值：

```
Loss 2.339936 | Fron 1.082749 | Pos 0.426216 | Cls 0.525598
              | Edge 0.001540 | Real 0.256227 | Fake 0.047607 | Surf 0.000000
```

| 项 | 含义 |
|---|---|
| `Fron` | 前沿/焦点预测（该在哪长） |
| `Pos` | 位置 MDN 负对数似然 |
| `Cls` | 元素分类 |
| `Edge` | 键级分类 |
| `Real` / `Fake` | **对比学习**：真实原子位置 vs 扰动位置的可分性 |
| `Surf` | 表面项（本配置下权重为 0） |

**对比学习细节**（`train.transform.contrastive`）：

```yaml
num_real: 20          # 每步采 20 个真实位置
num_fake: 20          # 采 20 个扰动位置
pos_real_std: 0.05    # 真实位置扰动很小
pos_fake_std: 2.0     # 假位置扰动很大（2 Å）
```

即：让模型学会"这个位置像真实原子位置"与"这个位置明显偏了"的区分，等价于在坐标空间做判别式正则。

### 3.6 训练数据构造：掩码策略

生成模型训练的本质是"挖空补全"：把部分原子遮住，让模型预测被遮住的部分。
`train.transform.mask` 使用 **mixed（混合）策略**：

```yaml
type: mixed
min_ratio: 0.0 / max_ratio: 1.1
min_num_masked: 1 / min_num_unmasked: 0
weights:  random 0.15 | BFS 0.60 | invBFS 0.25
```

| 掩码方式 | 权重 | 直觉 |
|---|---|---|
| random | 0.15 | 随机挖空，覆盖非连续缺失 |
| **BFS** | **0.60** | 从种子出发广度优先挖空 —— 模拟"生长中的前沿"，与采样时的自回归顺序一致 |
| invBFS | 0.25 | BFS 的逆序挖空，覆盖反向生长 |

生成式采样的顺序是"从中心向外长"，所以 BFS 掩码权重最高（0.6）—— **训练掩码分布与推理生长顺序对齐**，
这是该模型能逐步生成的关键设计。

边采样器 `edgesampler.k = 8`：每步在候选边中采样 8 条用于计算损失，控制显存开销。

### 3.7 优化与稳定性

```yaml
seed: 2022 ; init_checkpoint: ./ckpt/pretrained_Pocket2Mol.pt
batch_size: 1 ; num_workers: 0 ; pin_memory: false
max_iters: 2000 ; val_freq: 200 ; pos_noise_std: 0.1 ; max_grad_norm: 100.0
optimizer: adam  lr 2e-5  betas (0.99, 0.999)  weight_decay 0
scheduler: plateau  factor 0.6  patience 8  min_lr 1e-6
```

我们额外加入的稳定性措施（提交 `19f0b40`）：
`logsigma` 数值 clamp、`average_logp` 五元组格式对齐、空队列 break、ASCII 日志目录、
`init_checkpoint` 载入、`freeze_encoder` 选项、非有限 loss 跳过、`torch_cluster` 导入容错。

### 3.8 我们使用的权重

| 项 | 值 |
|---|---|
| 文件 | `ckpt/pretrained_Pocket2Mol.pt` |
| 大小 | 42.8 MB |
| 来源 | **Pocket2Mol 官方预训练权重**（我们没有替换它） |
| 使用处 | `configs/sample_for_pdb_guided_l3.yml → model.checkpoint`；`train_multitarget_v2.yml → train.init_checkpoint` |

**重要澄清（常被误解）**：`环境设置 → Pocket2Mol 仓库目录` 指向的是**我们的改进仓库**，
但仓库内部**加载的权重仍是官方基础权重**。我们的提升不来自权重，而来自**采样算法与后处理代码**
（第 4 章）。若把该路径指向官方原仓库，能加载同样的权重，但会失去 `guidance.py` 等全部改进。

我们自己微调出的权重位于 `logs/*/checkpoints/*.pt`（每份 19.9–42.9 MB），但**未被采样流程采用**，原因见 §4.4。

---

## 第 4 章 我们的改进（核心贡献）

### 4.1 引导束搜索（Guided Beam Search）★

#### 4.1.1 原版束搜索的缺陷

原版每一步从候选池中按模型概率随机保留 `beam_size` 个：

$$
p_i \propto \exp\!\Big(\sum \log p_{\text{model}}\Big) + 1
$$

问题：模型只知道自己"像不像真实分子"，并不知道"类不类药"。于是束搜索会稳定地
把化学空间探索导向模型训练集的平均偏好，**QED 等成药指标完全不受控**。

#### 4.1.2 我们的排序公式

在束排序概率上乘以一个**化学分数增益**：

$$
\boxed{\;p_i \;\propto\; \Big(\exp\big(\textstyle\sum \log p_{\text{model}}\big)+1\Big)\;\cdot\; w_i \;\cdot\; \exp\big(\lambda \cdot \text{chem}_i\big)\;}
$$

$$
\text{chem}_i = \underbrace{\text{QED}_i}_{w_{qed}=1}\cdot 1 \;+\; \Big(1-\frac{\min(SA_i,\,10)}{10}\Big)\cdot 1 \;-\; w_{div}\cdot \max_{j \in \text{ref}} \text{Tanimoto}(i, j)
$$

- $\lambda$：引导强度。**λ=0 时公式严格退化为原版**（$e^{0}=1$），这是可证伪性的基础；
- $\text{chem} \in [0, 2]$：QED（0–1）+ 合成可及性得分（0–1，SA 越小越高）；
- $w_i = 1/|\text{nexts}_i|$：父节点的候选数归一化权重（原版逻辑保留）；
- 最后一项是 §4.2 的多样性惩罚。

#### 4.1.3 实现要点（`utils/guidance.py`，130 行）

**难点**：束搜索中途的分子是"半成品"——价键不完整、可能超价。RDKit 默认 sanitize 会直接拒绝。

我们的做法 `build_partial_mol()`：

1. 用 `Chem.RWMol()` 手工建原子（`ligand_context_element`）、写构象（`ligand_context_pos`）；
2. 键类型映射 `{1: SINGLE, 2: DOUBLE, 3: TRIPLE}`，仅在 `i < j` 时加键，防止重复；
3. sanitize 时**显式关闭三项检查**：

```python
Chem.SanitizeMol(mol,
    Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE
                      ^ Chem.SANITIZE_SETAROMATICITY
                      ^ Chem.SANITIZE_PROPERTIES,
    catchErrors=True)
```

即：不要求芳香化、不要求价键守恒 —— 只要能算 QED/SA 即可。解析失败返回 `None`，
该候选得到一个中性的 fallback 分数（$0.5 \times (w_{qed}+w_{sa})$），**不惩罚也不奖励**，避免误伤。

**性能**：以 SMILES 为键的 `_SCORE_CACHE`（上限 200,000 条）缓存 QED/SA 结果；
SA 分数来自 `evaluation/sascorer.py`（RDKit Contrib 的合成可及性评分，1 = 易合成，10 = 极难）。

#### 4.1.4 关键工程决策：默认关闭

```python
guided_cfg = config.sample.get('guided', None)
if guided_cfg is not None and guided_cfg.get('enabled', False):
    ... 引导路径 ...
else:
    prob = logp_to_rank_prob(...)     # 原版路径，逐位一致
```

**为什么默认关闭**：① 保证原版行为可复现；② 让 A/B 实验中的"基线"是真的基线；
③ 一旦引导实现有 bug，用户去掉一行配置即可回到安全状态。

#### 4.1.5 参数定型：λ 消融与 A/B 证据

同一随机种子、同一口袋、同一束宽（50 样本 / 100 束宽 / 50 步）下的对照：

| 配置 | QED 中位 | 说明 |
|---|---|---|
| 基线（guided 关闭） | **0.667** | 原版行为 |
| 引导 λ=0.3 | 提升有限 | 增益过弱 |
| 引导 λ=1 | 中等提升 | — |
| **引导 λ=3** | **0.803（+20%）** | 定型参数 |
| 引导 λ=5 | 不再提升、多样性下降 | 过强导致塌缩 |

并做了**两个独立随机种子**（2024 / 2025）复现，确认提升不是种子偶然。
（配置矩阵实测记录见 `outputs/ab_5010050_*`、`guided_5010050_l*`。）

### 4.2 多样性惩罚（MMR 式）

**动机**：引导分数会让束搜索收敛到少数几个高分骨架，库的化学多样性崩塌。

`diversity_penalty(data_list, ref_smiles, div_w=0.5, max_ref=200)`：

1. 参照库 = **已完成分子集合**（`pool.smiles`）中的最近 200 个；
2. 每个参照分子算 Morgan 指纹（半径 2，`AllChem.GetMorganFingerprint`）；
3. 候选与参照库的**最大 Tanimoto 相似度**作为惩罚量，乘 `div_w` 后从 `chem` 中扣除；
4. 参照库为空或分子解析失败 → 惩罚 0（不误伤）；
5. 指纹缓存 `_REF_FP_CACHE`，仅当参照集合对象变化时重算。

这是最大边际相关（MMR）思想在**生成过程中**的在线版本——不是生成完再挑，而是**边生成边压制重复**。

**效果**：全库 2,931 个分子拥有 **1,960 个唯一 Murcko 骨架**（骨架/分子比 ≈ 0.67），
Top50 亦为 47 个唯一骨架（原版无惩罚时同类实验骨架集中度明显更高）。

### 4.3 生成构象精修 `relax_mol_geometry`

**动机**：模型逐步"贴"出来的原子，键长键角并不严格符合化学，直接交付会给出扭曲的三维结构。

`utils/reconstruct.py::relax_mol_geometry()` 的四步：

1. **ETKDG 重建内部几何**：以生成分子为拓扑，用 ETKDG 重新嵌入（`n_confs=4` 取最优）；
2. **姿态保真约束**：把重建构象与模型原姿态做对齐，`keep_pose=True` 时要求
   `pose_rmsd ≤ max_pose_rmsd`（默认 2.0 Å），否则回退到原姿态；
3. **MMFF94s 力场优化**（默认最多 500 次迭代）；
4. **保留显式氢交付**，并记录三个诊断量：`pose_rmsd`、`relax_energy`、`strain_per_heavy`（每重原子应变能）。

**设计取舍**：精修若完全放开，分子会"飞"出口袋——姿态保真阈值就是防止这一点。
该功能由 `config.sample.relax_output` 控制，**默认关闭**（与引导同样的零侵入原则）。

副产物：`strain_energy.csv`（全库应变能）成为后续筛选的一个独立维度。

### 4.4 微调路线（方向 2）—— 诚实记录：未晋级

**做了什么**：构建四靶点微调数据集 `gpcr_multitarget_v2`（**479 对**蛋白-配体复合物），
以官方权重为初始化，冻结编码器分阶段解冻，训练 2,000 步（`scripts/phase4_retrain.ps1`、`train.py`）。

**A/B 对照结果**：

| 指标 | 基线（官方权重） | 微调 2000 步 |
|---|---|---|
| 通过过滤的合格分子数 | 63 | 50 |
| QED 中位 | **0.597** | 0.408 |

**结论**：微调后**变差**，判定 `NOT_PROMOTED`，采样流程继续使用官方权重。

**为什么失败（技术归因）**：
1. 数据规模不足：479 对复合物对一个三维生成模型过少，模型在少量新分布上过拟合；
2. 官方预训练已在 CrossDocked2020（十万级）上完成，我们的数据量无法提供有效梯度信号；
3. 无法做官方等规模复训——仓库内 CrossDocked LMDB 是**空壳（0 条目）**，重训不可行。

**为什么值得写进全书**：这是项目可信度的证据。我们保留了完整判据
（`scripts/judge_promotion.py`）、原始日志与 `NOT_PROMOTED` 结论，并且**主动回滚**。
答辩材料中把它作为"科研诚信"的实证，而不是藏起来。

### 4.5 工程可靠性改进（同样属于"改进"，但不是算法）

| 改进 | 位置 | 收益 |
|---|---|---|
| `get_tri_edges` 向量化 | `models/maskfill.py` | 生成主循环 4× 提速，与原逐循环实现**逐位等价**（`scripts/test_tri_edges_equiv.py` 验证），保留回退分支 |
| `get_init` 非有限概率防护 | `sample_for_pdb.py` | 根除"全 NaN 概率 → 死循环"路径（重试上限 10，返回空列表） |
| `sample_position` 显式报错 | `models/maskfill.py` | 原先 `assert` 在 `-O` 优化下会被剥离而静默返回未定义变量，改为显式 `raise NotImplementedError` |
| 去裸 `except` | `utils/reconstruct.py`、`utils/transforms.py` | 把吞异常的 `except:` 改为 `except Exception`，异常不再被静默吞掉 |
| BOM 根除 | `utils/misc.py::load_config` | 配置读取统一 `utf-8-sig`，带 BOM 的 YAML 不再导致 `config['\ufeffmodel']` 崩溃 |
| 结构化配置生成 | `scripts/gen_sample_config.py` | 取代"在模板末尾追加缩进行"的危险做法（后者在模板结构变化后产出非法 YAML） |

---

## 第 5 章 推理全流程（从口袋到候选库）

> 依据 `sample_for_pdb.py`（225 行）逐段解读。

### 5.1 口袋提取

```python
data = pdb_to_pocket_data(args.pdb_path, args.center, args.bbox_size)
```

- `center`：口袋中心（四个靶点的中心已核实，见 §12.2）；
- `bbox_size`：口袋盒子边长，默认 **23.0 Å**（覆盖 GPCR 正构口袋的典型尺度）。

### 5.2 初始化采样

```python
pool = {'queue': [], 'failed': [], 'finished': [], 'duplicate': [], 'smiles': set()}
init_data_list = get_init(data, model=model, transform=atom_composer,
                          threshold=config.sample.threshold)
pool.queue = init_data_list[:config.sample.beam_size]
```

`get_init` 采样"第一批原子"（生长起点），并用三组阈值过滤：
`focal 0.5 / pos 0.25 / element 0.3`。此处即我们加入 NaN 防护的位置。

### 5.3 主循环（自回归生长 + 束裁剪）

```python
while len(pool.finished) < config.sample.num_samples:      # 目标：凑够成品数
    global_step += 1
    if global_step > config.sample.max_steps: break        # 安全上限，不是目标步数

    for data in pool.queue:                                # 对束中每个候选
        data_next_list = get_next(data, model, transform, threshold)
        for data_next in data_next_list:
            if data_next.status == STATUS_FINISHED:        # 该分子已终止
                rdmol = reconstruct_from_generated_with_edges(data_next)
                if config.sample.get('relax_output', False):
                    rdmol, info = relax_mol_geometry(rdmol, keep_pose=..., max_pose_rmsd=...)
                smiles = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rdmol)))
                if smiles in pool.smiles:      pool.duplicate.append(...)   # 去重
                elif '.' in smiles:            pool.failed.append(...)      # 断片，判失败
                else:                          pool.finished.append(...)    # 通过
            elif data_next.status == STATUS_RUNNING:
                nexts.append(data_next)
        queue_weight += [1./len(nexts)] * len(nexts)       # 父节点归一化

    if len(queue_tmp) == 0: break                          # 池枯竭，提前结束
    prob = <引导或原版束排序概率>
    next_idx = np.random.choice(len(queue_tmp), p=prob, size=min(beam_size, n), replace=False)
    pool.queue = [queue_tmp[i] for i in next_idx]
    torch.save(pool, 'samples_%d.pt' % global_step)         # 每步快照，可中断恢复
```

**四个要点**：

1. **`num_samples` 与 `max_steps` 的区别**（易被误读）：
   `num_samples` 是"希望凑够多少个成品"，`max_steps` 是"最多长多少步"的安全阀。
   把 `max_steps` 设小（如 15）会导致分子没长完就 break，成品数远低于目标 ——
   **这不是故障**，是参数语义。E2E 验证脚本即因早期误设此项而产生"每次只有 4 个分子"的困惑。
2. **四态状态机**：`queue`（生长中）/ `finished`（成品）/ `failed`（重建失败或含断片）/ `duplicate`（与已成品重复）。
   所有状态都被记录并计数，**不静默丢弃**。
3. **快照**：每一步存 `samples_N.pt`（含整个池），使长时运行可中断、可诊断、可续算。
4. **`np.random.choice(..., replace=False)`**：无放回抽样，避免同一父代被重复选中导致束内近亲繁殖。

### 5.4 输出

| 产物 | 说明 |
|---|---|
| `SMILES.txt` | 每行一个成品的规范 SMILES |
| `SDF/*.sdf` | 与 SMILES 行序对应的三维结构（0 基编号） |
| `samples_*.pt` | 每步束池快照（含 logp 明细，可复算排序概率） |
| `log.txt` | 后端日志（`[Pool] Queue/Finished/Failed` 行被 GUI 解析为进度） |
| `config.yml` / `<pdb>.pdb` | 本次运行参数与口袋文件的存档副本 |

### 5.5 过滤入库漏斗（`scripts/build_library.py`，566 行）

实测一次正式参数运行（A2A，50 样本 / 100 束宽 / 50 步）的漏斗：

```
raw SMILES lines read                       53
invalid/unparseable (dropped)                0
unique canonical SMILES (deduped)           53
after PAINS                                 51   (-2, 96.2%)
after BRENK                                 45   (-6, 84.9%)
after MW in [250, 500]                      45
after Crippen LogP in [1, 5]                45
after QED >= 0.4                            45
after SA score <= 6                         45
after elements CNOFPSCl                     45
after rings >= 1                            45
final library compounds                     45
with 3D SDF (pdbqt_ready)                   45
```

被拒分子连同原因写入 `rejected.csv`（**可审计**，而非直接消失）。
全库（2,931）层面的对应统计：PAINS/BRENK 合计剔除约 2,000 条原始记录，
最终 6,864（四靶点）进入化合物库。

---

## 第 6 章 评估与筛选体系

> 11 项脚本化评估，全部产物落盘于 `D:\MMModel\化合物库\`。核心思路：
> **计算层面的"好"必须能被多个相互独立的证据链交叉验证**。

### 6.1 分子对接（Vina）

`scripts/docking_pipeline.py`（361 行）：

| 项 | 设置 |
|---|---|
| 引擎 | AutoDock Vina 1.2.5 + Open Babel 3.2.1（便携部署于 `D:\MMModel\tools\`） |
| 口袋盒 | 22.5 Å 立方，中心由靶点给定 |
| 受体准备 | 自动去水、去共晶配体、加氢、转 PDBQT |
| 配体准备 | SMILES → 3D（ETKDG）→ 加氢 → PDBQT |
| 并行 | 多进程，单机批量 |

**A2A 结果（2,934 条记录，成功率 100%）**：

| 指标 | 值 |
|---|---|
| 均值 | **−9.68 kcal/mol** |
| 最强 | −13.51 |
| ≤ −9（强结合） | 1,142（38.8%） |
| ≤ −10 | **1,137** |
| ≤ −11（极强） | 311（10.6%） |
| Top1 | −13.42（QED 0.797） |

### 6.2 新颖性

- 对照集：**737 个已知 A2A 高活性分子**（ChEMBL，含 362 个 Murcko 骨架）+ 上市药物；
- 指标：Morgan(2) 最大 Tanimoto；
- Top50 结果：**最大相似度 0.187 – 0.321**，全部判为 `HIGH_NOVEL`，**零骨架重合**；
- 分类器：`KNOWN_SCAFFOLD`（骨架命中）/ `NEAR_KNOWN`（≥0.85）/ `mid_similarity`（≥0.4）/ `HIGH_NOVEL`（<0.4）。

### 6.3 逆向阳性对照（管线可信度）

把 **737 个已知活性分子**送入**同一条**对接管线：

| 集合 | 均值 |
|---|---|
| 已知 A2A 活性分子（阳性对照） | −9.32 |
| 我们生成的库 | **−9.68** |

**双重结论**：① 管线对已知活性分子给出合理强打分 → 对接设置正确（正向验证）；
② 生成库的均值**优于**已知活性分子集合 → 生成确实有效（不只是"能跑"）。

### 6.4 富集验证（ROC / 诱饵）

`scripts/roc_decoys.py`：以已知活性分子为正样本、**1,408 个诱饵**为负样本，
用对接分做排序，计算 ROC 与富集因子。

| 指标 | 值 |
|---|---|
| ROC-AUC | **0.685**（另一版统计 0.680，与文献 Vina 水平一致） |
| 富集倍数（≤ −12） | **5.80×** |
| 富集倍数（≤ −11） | 4.48× |
| top-1% 富集 | 2.24× |

意义：证明**对接打分本身具备区分力**——否则后文所有"强结合"结论都没有意义。

### 6.5 对接姿态置信度

`scripts/pose_consistency.py`：Vina 每个分子返回多个模式（mode），比较最优与次优模式的能量差
与 RMSD。

| 指标 | 值 |
|---|---|
| 全库中"最优/次优能量差 < 0.5 kcal/mol"的比例 | **64%** |
| Top50 中"可信姿态"（差值大、consensus 明确）占比 | 29%（全库仅 6%） |

**解读**：多数分子存在多个近简并姿态 → 刚性受体对接的固有局限；
但 Top50 的可信比例是全库的约 5 倍，说明**排序在挑选结合模式更明确的分子**，这是筛选有效性的又一证据。

### 6.6 ADMET 规则筛查

`scripts/admet_screen.py`：基于理化性质与结构规则的**风险分层**（不是预测模型）：

| 级别 | 含义 | 全库 | Top50 |
|---|---|---|---|
| A | 无警示 | **676** | 7 |
| B | 轻微警示 | 1,940 | 39 |
| C | 明显警示（如 hERG 风险片段） | 315 | 4 |

### 6.7 多目标 Pareto 分析

`scripts/pareto_analysis.py`：五个目标同时优化（对接分、QED、选择性、可测性、ADMET 警示数），
求非支配解集：

| 指标 | 值 |
|---|---|
| 前沿成员 | **91** |
| 旧 Top50 中落在前沿上的 | 仅 **27**（即旧排序漏掉了大量更优权衡点） |
| 前沿成员的选择性中位 | +21.3 百分位点 |

**意义**：证明"单一加权排序（0.6×打分+0.4×QED）"确实次优，多目标是必要的。

### 6.8 选择性谱（A2A vs A1）

`scripts/selectivity_analysis.py`：对全库分子同时对接 A2A（4EIY）与 A1（5UEN），
以百分位差定义选择性指数 SI：

| 指标 | 值 |
|---|---|
| 全库 SI 中位 | **+4.0** |
| SI > 0 的比例 | 59% |
| SI > 25 的分子 | **382** |
| Top50 中倾向 A2A 的比例 | 88% |

意义：GPCR 药物最大的失败原因是**选择性不足**（A2A 与 A1 同源度高）。
在生成库层面给出选择性分布，是把"能不能成药"往前推了一步。

### 6.9 实验可测性约束

`scripts/assay_rules.py`：荧光干扰、反应性片段、ESOL 溶解度等规则，
并用已知集回归校准（红标比例从 63.8% 降到 6.8%，即规则不是"随便打标签"）。
输出 `top50_assayranked.csv`。

### 6.10 模型自洽性 / 域外检测

`scripts/selfconsistency.py`：以参照库做 kNN 距离，把全库分为
inner 1,466 / transition 732 / outer 733 —— 落在训练分布外侧的分子标 `outer`，
提示其可靠性较低。这是**对自己模型的不确定性建模**。

### 6.11 六维综合推荐

`scripts/final_recommendation.py`：综合六个维度（对接分、QED、SA、MW、选择性、可测性/ADMET）
加权打分并分级：

| 级别 | 全库数量 |
|---|---|
| Tier 1 | 166 |
| Tier 2 | 454 |

导出清单 `final_recommendation.csv` 保存 **200 条**（Tier1 前 100 + Tier2 前 100），
供湿实验送验使用。

---

## 第 7 章 数据资产与结果汇总

### 7.1 微调数据集（`data/`）

| 数据集 | 规模 | 用途 | 备注 |
|---|---|---|---|
| `gpcr_a2a` | 136 | A2A 单靶点微调 | 拓扑孤立原子/构象重叠清洗后 |
| `gpcr_b2ar` | 179 | β2-AR 单靶点 | 来自 2RH1 |
| `gpcr_d3` | 91 | D3 单靶点 | 来自 3PBL |
| `gpcr_5ht2b` | 73 | 5-HT2B 单靶点 | 来自 4IB4 |
| `gpcr_multitarget_v2` | **479** | 四靶点联合微调 | 主力微调集（清洗后） |
| `gpcr_multitarget` | — | 第一版（被 v2 取代） | 保留供对照 |

构建链路：`scripts/整理配体数据.py` → `build_gpcr_dataset.py` → `merge_gpcr_datasets.py`
→ `purge_orphans.py`（清洗）→ `verify_gpcr_datasets.py`（校验 index↔split 一致性）
→ `kekulize_dataset.py` / `audit_bondtypes.py`（键类型审计）。

每个数据集目录含一个 LMDB 与 `split_by_name.pt`（**按名称划分**，避免同一复合物泄漏到验证集）。

### 7.2 化合物库总表（四靶点）

| 靶点 | PDB | 原始产出 | **入库** | 对接均值 | 最强 | 强结合占比 | Top50 新颖 |
|---|---|---|---|---|---|---|---|
| A2A 腺苷受体 | 4EIY | 4,231\* | **2,931** | −9.68 | −13.51 | 38.8% | 50/50 |
| β2 肾上腺素受体 | 2RH1 | 1,577 | **917** | −9.78 | −13.97 | 38.9% | 49/50 |
| D3 多巴胺受体 | 3PBL | 1,571 | **951** | −8.40 | −11.39 | 5.5% | 49/50 |
| 5-HT2B 血清素受体 | 4IB4 | 1,525 | **1,065** | −9.34 | −12.34 | 24.9% | 50/50 |
| **合计** | — | **10,004** | **6,864** | — | −13.97 | — | 198/200 |

\* A2A 的原始产出按首轮正式批次统计为 4,231；若把后续补充批次全部计入，全量重算为 4,809。
入库数 2,931 为权威值（与 `compounds_tagged.csv` 行数一致）。

**注意 D3 的 5.5%**：这反映受体口袋的可成药性差异（D3 正构口袋更浅、更暴露），
是一个**真实信号**而非流程故障——它的对接均值（−8.40）同样明显弱于其它三个靶点。

### 7.3 主库质量画像（2,931）

| 指标 | 值 |
|---|---|
| QED 中位 / 最大 | **0.825** / 0.947 |
| SA 中位（范围） | **3.33**（1.52 – 5.95，全部可合成区间） |
| MW 中位 | 362 Da（口服药物甜区） |
| 唯一 Murcko 骨架 | **1,960** |
| QED ≥ 0.8 占比 | 76% |
| QED 分布 | 0.4–0.6: 56 / 0.6–0.8: 1,098 / 0.8–1.0: 1,770 |

### 7.4 关键数字速查

```
生成：QED 中位 0.825 · SA 3.33 · MW 362 · 骨架 1,960 · 库 2,931（四靶点 6,864）
引导：λ=3 → QED 0.667→0.803（+20%），双种子复现
对接：均值 −9.68 · ≤−10 共 1,137（38.8%）· ≤−11 共 311 · 最强 −13.51 · Top1 −13.42
对照：已知集 737 个 → −9.32（生成库更优）
富集：AUC 0.685 · ≤−12 富集 5.80× · top-1% 2.24×
姿态：全库 64% 存在近简并姿态 · Top50 可信率 29% vs 全库 6%
成药：ADMET A676 / B1940 / C315 · Top50 A7 / B39 / C4
多目标：Pareto 前沿 91 · 旧 Top50 仅 27 在沿上
选择性：SI 中位 +4.0 · 59% >0 · 382 个 >25 · Top50 88% 倾向 A2A
推荐：Tier1 166 / Tier2 454 · 导出清单 200 条
工程：E2E 自检 8/8 阶段通过 · exe 启动 22s → 1.3s
```

---

## 第 8 章 工程实现

### 8.1 桌面应用（`gui/`）

**技术选型**：Python 标准库 Tkinter（零第三方依赖）。原因：需要打包成免安装 exe，
引入 Qt/PySide 会让体积与打包复杂度上升一个量级。

| 文件 | 行数 | 职责 |
|---|---|---|
| `app_pro.py` | ~1,350 | v2.0 Pro 主程序（4 标签页） |
| `ui_kit.py` | ~350 | 自绘组件库：圆角按钮、靶点卡片、分区标题、双主题、字体度量 |
| `app.py` / `app_easy.py` | 950 / 390 | 早期版本（保留） |
| `make_icon.py` | 120 | 生成多尺寸 ico/png/base64 图标 |
| `build_exe.ps1` | 95 | 双模式打包脚本 |

**四个标签页**：

| 标签页 | 功能 |
|---|---|
| 生成分子 | 靶点卡片选择 → 开始生成 → 实时进度 → 候选表（排序/阈值筛选/搜索）→ 2D 结构预览 → 导出 SDF/CSV |
| 候选库浏览 | 打开化合物库任意 CSV（26 个数据文件），排序、看结构、导出视图 |
| 分析结果 | 汇总化合物库 CSV、对接 summary、桌面报告，双击用系统程序打开 + 前 60 行预览 |
| 环境设置 | 路径可改并保存到 `gui_config.json`；一键环境自检（11 项，含 GPU） |

顶栏：`?` 帮助弹窗（随主题着色、跟随主窗口居中）、`☾/☀` 浅色/深色切换、环境状态灯。

**交互实现细节（踩过的坑）**：

1. **所有自绘控件的宽度按字体实测值计算**（`ui_kit.measure` 用 `tkfont.Font.measure`）。
   固定像素宽在 131%/175% 缩放下会把按钮文字裁掉——这是 v1 界面的主要缺陷。
2. **Canvas 属性不能用 `self._w`/`self._h`**：这两个名字是 tkinter 内部 widget 路径，
   覆盖后报 `invalid command name "96"`。
3. **`TreeviewSelect` 回调必须接收 event 参数**：`def _show_mol(self, event=None)`。
   否则点击行查看 2D 结构会静默失败（异常被事件循环吞掉）——v1 与 v2 早期都存在此缺陷。
4. **Label 显示图片前要清除字符单位 width/height**：`configure(image=img, width=img.width(),
   height=img.height())`，否则 440×360 的结构图只露出左上角一角。
5. **DPI 策略**：源码运行声明 per-monitor 感知并按 DPI 放大窗口；打包版交由 Windows 统一缩放
   （避免字体二次放大）。
6. **环境自检后台执行**：自检里含 `import torch` 的子进程（5–15 秒），
   原先同步执行会把窗口显示一起卡住——改为后台线程 + `root.after` 回填，是启动优化的主因。
7. **关窗协议**：采样中关闭窗口二次确认并终止子进程，避免孤儿进程占 GPU 继续写盘。

### 8.2 打包与分发

**双模式产物**：

| 模式 | 路径 | 启动 | 体积 | 用途 |
|---|---|---|---|---|
| 单文件 | `gui\dist\7-eonmol.exe` | 2.6 s | 7.7 MB | 只发一个文件的场景 |
| 快速启动 | `gui\dist\fast\7-eonmol\` | **1.3 s** | 171 文件 / 15.6 MB | 日常使用 / 便携 zip |
| 自包含部署包 | `桌面\7-eonmol\` | 2.6 s | 5.09 GB | 目标机无 Python 环境 |

**启动优化过程（22 秒 → 1.3 秒）**：

| 措施 | 效果 |
|---|---|
| 环境自检异步化（根因） | −15 s |
| 排除无关模块（numpy/torch/rdkit/PIL/ssl/hashlib/sqlite3/setuptools…） | 单文件 9.9 → 7.7 MB |
| 裁剪 Tcl/Tk 冗余数据（tzdata 时区 500+ 文件、msgs 多语言、多余 encoding） | 快速版 983 → 171 文件 |
| onedir 免解压 | 再次启动 1.31 s |

> 实测（模拟分发到新目录，各 3 次）：快速启动版 3.51 / 1.33 / 1.32 s；
> 单文件版 3.43 / 2.42 / 2.46 s。单文件每次要把 7.7 MB 压缩包解压成 **50 MB** 到 `%TEMP%`，
> 这就是那 1.1 秒差异的来源。**结论：分发推荐用便携目录版 + zip。**

**自包含部署包结构**：

```
7-eonmol\
├─ 7-eonmol.exe        ← 优化后的单文件版
├─ 一键配置.bat              ← 首次运行：conda-unpack 适配路径
├─ README.md                 ← 部署与使用说明
├─ repo\                     ← 改进版仓库副本（167 MB，含 ckpt 与靶点结构）
└─ env\                      ← conda-pack 环境（5.0 GB，torch+cu124+RDKit）
```

部署模式的关键设计：exe 启动时若发现**同目录**存在 `repo\` 与 `env\`，自动改用包内资源
（`_EXE_DIR` 探测），因此整个文件夹可拷到任意机器。

**本轮修复的部署缺陷**：旧包硬编码 `repo\structures\` 作为靶点目录，而仓库里结构实际在
`repo\targets\` —— 导致旧部署包启动后**四个靶点全部"未找到"**，只能开界面不能生成。
已改为按候选目录（`structures` → `targets` → `pdb` → `data/structures`）自动探测，
实测 4/4 找到。

### 8.3 脚本体系（38 个 `scripts/`）

| 类别 | 脚本 |
|---|---|
| 数据构建 | `整理配体数据.py`、`build_gpcr_dataset.py`、`merge_gpcr_datasets.py`、`purge_orphans.py`、`kekulize_dataset.py`、`verify_gpcr_datasets.py`、`audit_bondtypes.py`、`audit_conformers.py`、`cif2pdb.py`、`calc_pocket_centers.py` |
| 生成与配置 | `gen_sample_config.py`、`night_phase.ps1`、`overnight_a2a.ps1`、`overnight_phase2.ps1`、`phase4_retrain.ps1`、`dock_3targets.ps1`、`run_sample.bat` |
| 入库与库维护 | `build_library.py`、`norm_docking_summary.py`、`relax_library.py`、`regen_top50.py` |
| 对接 | `docking_pipeline.py` |
| 分析（11 项） | `top50_novelty_3targets.py`、`roc_decoys.py`、`pose_consistency.py`、`admet_screen.py`、`pareto_analysis.py`、`selectivity_analysis.py`、`assay_rules.py`、`selfconsistency.py`、`final_recommendation.py`、`eval_candidates.py`、`test_strain_proxy.py` |
| 验证与判定 | `verify_all_fixes.py`、`verify_hardening.py`、`test_tri_edges_equiv.py`、`e2e_full_run.py`、`judge_promotion.py` |

### 8.4 夜间自动化批跑

`scripts/overnight_a2a.ps1` / `overnight_phase2.ps1`（PowerShell）：

- 分批生成 → 分批入库 → 分批对接 → 分析，全链路无人值守；
- **磁盘守卫**：每批结束清理中间快照（早期曾因 `samples_*.pt` 堆积把磁盘写满）；
- **决策解析**：用锚定正则 `'DECISION:\s*PROMOTED\s*$'` 判定晋级
  （早期 `-match 'PROMOTED'` 会把 `NOT_PROMOTED` 误判为晋级）；
- **等待超时**：14 小时上限、跨日安全（`$end.AddDays(1)`）；
- **数值排序**：checkpoint 按数值而非字典序取最新；
- **编码**：`.ps1` 必须带 UTF-8 BOM，否则 Windows PowerShell 5.1 按 GBK 读取会语法报错。

---

## 第 9 章 迭代进程史

> 以 git 提交为主线（73 次提交，仓库 fork 自 Pocket2Mol 官方）。

### 9.1 时间线

| 日期 | 阶段 | 关键动作 |
|---|---|---|
| 08-09 | **接管与兼容** | fork 官方仓库；PyTorch 2.6 / Windows 兼容补丁；快速测试配置 |
| 08-14 | 资产整理 | 训练文件清单、靶点 PDB、配体数据、已知药物库入库 |
| 08-28 | **方向 1：算法改进** | 实现 `utils/guidance.py` 引导束搜索（配 `guided` 开关，默认关闭）；GPCR 微调数据管线；训练/采样稳定性修复 8 项 |
| 08-29 | 数据治理 | 数据集拓扑孤立原子循环清洗；`transforms` BFS 容错；审计工具 |
| 08-30 | 规模生成 + 微调 | 三靶点生成 4,673 个；四靶点微调 2,000 步（**判定不晋级**）；`purge_orphans`；三靶点对接总控；**外行友好版 GUI + 首个 exe** |
| 09-16 | 部署 | exe 部署模式（自动发现包内 repo/env）；三靶点 Top50 新颖性 |
| 09-19 | **任务 6 / 9** | 实验可测性约束引擎（含已知集回归校准）；模型自洽性 kNN 域外检测 |
| 09-20 | **任务 7 + 分析扩展** | 选择性谱（全库双打分）；构象精修 `relax_mol_geometry`；`tri_edges` 向量化（4×）；`get_init` NaN 防护；条件性缺陷加固 7 项；C11 ADMET；A3 Pareto；C10 姿态置信度；六维送验推荐；C9 ROC 收官 |
| 09-21 | **质量攻坚 + 交付** | 三批全工程 bug 修复（配置生成/静默丢样/编码/数值）；11 份下游分析全部重算；`verify_all_fixes.py` 36/36；PPT 38 页；**GUI v2.0 Pro**；部署包同步 |
| 09-22 | 交付打磨 | PPT 压缩至 20 页；GUI 启动优化（22 s → 1.3 s）；修复点击行不显示 2D 结构；更新自包含部署包并修复其靶点目录缺陷 |

### 9.2 关键决策记录

| 决策 | 内容 | 理由 |
|---|---|---|
| **D1 改进零侵入** | 引导束搜索默认关闭，配置不存在时逐位等同原版 | 保证基线真实、改进可证伪、失败可回退 |
| **D2 微调不晋级** | A/B 判据下微调更差 → 停用，如实汇报 | 科研诚信；避免把噪声当提升 |
| **D3 标注优先于删除** | 可疑分子打标签（ADMET/可测性/tier），不物理删除 | 保留化学空间，让下游按需筛选 |
| **D4 一切结论脚本化** | 每个数字都有生成脚本与 CSV | 可追溯、可复现、可反驳 |
| **D5 修复必须带验证** | 每轮 bug 修复都配验证脚本（36 项 / 4 项 / 端到端 8 阶段） | 防止"修好了 A 坏了 B" |
| **D6 失败要吵** | 静默丢弃改为显式告警/计数/落盘 | 静默丢样是本项目最隐蔽的一类缺陷 |

### 9.3 数据修复对结论的影响（稳健性验证）

09-21 的修复中，对接管线有 3 处静默丢样（中文路径 SDF、显式氢 SMILES、BOM 首行），
累计造成每靶点 1 个候选丢失。修复后**重算了全部 11 份下游分析**：

> 结论：对接均值、富集倍数、Pareto 前沿、选择性分布等**全部结论方向不变**，
> 数值变动在小数点后第 2 位量级。

这一步的价值在于：证明了此前基于"略有缺失的数据"得出的结论是稳健的。

---

## 第 10 章 缺陷审计与修复方法论

### 10.1 缺陷分类与代表案例

| 类别 | 代表缺陷 | 根因 | 修复 |
|---|---|---|---|
| **静默丢样** | 对接管线 3 处：① RDKit C++ 文件 API 打不开含中文路径的 SDF；② 显式氢 SMILES 解析失败；③ 文件首行 BOM 污染 | 异常被 `except` 吞掉，返回空 SMILES 后照常写入统计 | 改用 `open()` + `MolFromMolBlock`；SMILES 规范化；`utf-8-sig` + `lstrip('\ufeff')`；并写 `norm_docking_summary.py` 修复 7 行历史污染数据 |
| **配置生成** | GUI/脚本"在 YAML 模板末尾追加缩进行" | 模板新增字段（如 `relax_output`）后，追加位置落到错误层级 → `ScannerError` | 新增结构化生成器 `gen_sample_config.py`（PyYAML 读改写 + 回读自检） |
| **编码** | ① 带 BOM 的配置文件 → `config['\ufeffmodel']` 崩溃；② `.ps1` 无 BOM 被 GBK 误读；③ 新版 RDKit 移除 `rdkit.six` 而 `sascorer.py` 依赖它 | Windows 中文环境的多编码混用 | `load_config` 统一 `utf-8-sig`；脚本显式保存为 UTF-8 BOM；内联 `rdkit.six` shim |
| **数值健壮性** | ① `get_init` 全 NaN 概率 → 死循环；② `sample_position` 的 `assert` 在 `-O` 下被剥离 → 返回未定义变量；③ 非有限 loss | 边界条件未防护 | 重试上限 + 空列表返回；显式 `raise NotImplementedError`；非有限 loss 跳过 |
| **条件性缺陷** | 7 项：芳香键契约缺失、索引与 split 不一致、PCA 轴对齐判据不可靠等 | 只在特定输入下暴露 | 逐项加固（`verify_hardening.py` 4/4 PASS） |
| **部署** | ① 打包用系统 Python 3.14 的 Tk 9 → 缺 DLL；② 部署模式硬编码 `structures\` 而结构在 `targets\` | 打包环境与运行环境不一致；路径假设未验证 | 用 conda Tk 8.6 + 显式注入 tcl/tk DLL；靶点目录自动探测（4/4 验证） |
| **GUI** | ① `TreeviewSelect` 回调签名少 event → 点击行永不渲染结构；② 预览 Label 字符单位裁剪图片；③ 自绘控件固定像素宽 → 高 DPI 裁字 | Tkinter 细节与 DPI 假设 | 逐项修复并实拍验证（见 §8.1） |
| **判据逻辑** | `-match 'PROMOTED'` 把 `NOT_PROMOTED` 判成晋级 | 子串匹配 | 锚定正则 `'DECISION:\s*PROMOTED\s*$'` |

### 10.2 审计方法

1. **静态扫描**：全仓库检索 `except:`（裸异常）、`assert`、硬编码路径、`stdout=PIPE`（沙箱禁用）；
2. **动态取证**：扫描 237 份历史运行日志，统计真实异常（仅 1 类：CUDA CUBLAS，已修复）；
3. **等价性测试**：对性能优化（`tri_edges` 向量化）写逐位等价测试，确保优化不改变数值结果；
4. **端到端验证**：`e2e_full_run.py` 串起 8 个阶段（环境 → 配置 → 采样 → 入库 → 对接 → 新颖性 → 训练冒烟 → 下游分析链），
   全部通过才算健康（实测 8/8）；
5. **回归清单**：`verify_all_fixes.py` 将 36 条断言固化为可重复执行的检查（36/36 PASS）。

### 10.3 一条重要的方法论教训

**`sanitize=True` 会掩盖真实键值**：早期怀疑数据集的键类型被芳香标志污染，
但直接读 LMDB 得到的全是 kekulize 后的单/双/三键（30,304 条键**零芳香标志**）。
原因是解析器的 `sanitize=True` 会把芳香键"修正"成 kekulize 形式 ——
**用带 sanitize 的接口做数据审计，会得到错误的安心**。
改用原始键值读取（`audit_bondtypes.py`）后才确认：训练特征无污染。

---

## 第 11 章 局限、未采纳路线与展望

### 11.1 已知局限（必须如实说明）

| 局限 | 具体表现 | 影响 |
|---|---|---|
| **刚性受体对接** | Vina 固定受体构象，不含诱导契合 | 亲和力估计有系统偏差 |
| **打分非自由能** | Vina 打分是经验性的，且姿态置信度显示 64% 分子存在近简并模式 | 排序可用，绝对值不可当 ΔG |
| **ADMET 是规则筛查** | 基于理化与结构片段的警示分层，不是预测模型 | 只能筛风险，不能给 ADMET 端点预测 |
| **选择性基于对接分差** | A2A vs A1 的百分位差，不是实验选择性 | 提示性结论 |
| **无湿实验验证** | 全部为计算结论 | 最终成药性未经实验确认 |
| **微调未成功** | 479 对数据不足；CrossDocked LMDB 为空壳（0 条目），官方规模复训不可行 | 模型本体未针对 GPCR 适配 |
| **两处待确认** | hERG `!$(Na)` SMARTS 语义；mmCIF 引号字段解析 | 影响面小，已标注 |

### 11.2 评估过但未采纳的路线

| 路线 | 为什么没做 |
|---|---|
| A1 对接引导搜索（把 Vina 分数放进束排序） | 每步都调 Vina，成本不可接受（比生成本身贵 2–3 个数量级） |
| A2 DPP 行列式点过程多样性 | 多样性目标已由 Tanimoto 惩罚覆盖，收益边际 |
| A5 ε-greedy 束搜索 | λ 扫描已找到合适强度，探索性问题非当前瓶颈 |
| B6/B7 强化学习微调、LoRA | 微调本身未晋级，先解决数据问题再谈算法 |
| B8 实验构象微调数据 | 依赖湿实验反馈，暂不具备条件 |
| C10 全库分子动力学 | 2908 个分子做 MD 需数千 GPU 小时，与赛程不匹配 |
| D12 全 GPCR 家族扩展 | 当前四靶点已足够验证通用性（同一管线零失败），扩展属于工程复制 |

### 11.3 后续工作优先级

1. **实验闭环**：Top200 清单送验 → 用实验构象与活性数据重建微调集（把 §4.4 的失败路线变成可行路线）；
2. **亲和力精算**：对 Top50 做 MM-GBSA / 短程 MD 复核，替代当前 Vina 单点打分；
3. **主动学习循环**：湿实验反馈 → 引导分数权重自动演化 → 重新生成；
4. **选择性扩展**：A1/A2A 之外增加 A2B、A3 亚型，做亚型选择性谱；
5. **工程收敛**：把 38 个脚本编排为一条带断点续跑与磁盘守卫的流水线（当前夜间脚本已具备雏形）。

---

## 第 12 章 附录

### 12.1 四个靶点参数

| 靶点 | PDB | 口袋中心 (x,y,z) | 适应症 |
|---|---|---|---|
| A2A 腺苷受体 | 4EIY | −0.4, 8.5, 17.1 | 帕金森病 · 肿瘤免疫 · 炎症 |
| β2 肾上腺素受体 | 2RH1 | −29.5, 9.2, 6.9 | 哮喘 · 心血管疾病 |
| D3 多巴胺受体 | 3PBL | 0.085, −14.828, 10.432 | 精神分裂症 · 成瘾 |
| 5-HT2B 血清素受体 | 4IB4 | 22.448, 18.284, 11.726 | 偏头痛 · 肺动脉高压 |
| （选择性对照）A1 腺苷受体 | 5UEN | 55.969, 58.897, 143.624 | A2A 选择性评估用 |

中心值来源：4EIY/2RH1 的参考配体质心；3PBL/4IB4 用同一方法计算（方法一致性已验证）。

### 12.2 配置速查

**采样（`configs/sample_for_pdb_guided_l3.yml`）**

```yaml
model: { checkpoint: ./ckpt/pretrained_Pocket2Mol.pt }
sample:
  seed: 2024 ; num_samples: 50 ; beam_size: 100 ; max_steps: 50
  threshold: { focal_threshold: 0.5, pos_threshold: 0.25, element_threshold: 0.3,
               hasatom_threshold: 0.6, bond_threshold: 0.4 }
  guided: { enabled: true, qed_w: 1.0, sa_w: 1.0, lam: 3 }
  relax_output: true
```

**训练（`configs/train_multitarget_v2.yml`）** 关键项见 §3.3 / §3.6 / §3.7。

### 12.3 复现命令

```powershell
# 1) 生成（A2A，正式参数）
D:\Miniconda3\envs\Pocket2Mol\python.exe sample_for_pdb.py `
  --pdb_path "D:\MMModel\靶点结构\4EIY_A2A受体.pdb" --center " -0.4,8.5,17.1" `
  --config configs\sample_for_pdb_guided_l3.yml --outdir outputs\my_run

# 2) 过滤入库
D:\Miniconda3\envs\Pocket2Mol\python.exe scripts\build_library.py `
  --runs outputs\my_run --library outputs\my_run\library

# 3) 对接
D:\Miniconda3\envs\Pocket2Mol\python.exe scripts\docking_pipeline.py `
  --receptor "D:\MMModel\靶点结构\4EIY_A2A受体.pdb" --center= -0.4,8.5,17.1 `
  --smiles-file library_smiles.txt --out outputs\docking_my

# 4) 全链路自检（8 阶段）
D:\Miniconda3\envs\Pocket2Mol\python.exe scripts\e2e_full_run.py

# 5) 修复回归验证（36 项）
D:\Miniconda3\envs\Pocket2Mol\python.exe scripts\verify_all_fixes.py

# 6) GUI
D:\Miniconda3\envs\Pocket2Mol\python.exe gui\app_pro.py            # 源码运行
gui\dist\7-eonmol.exe                                        # 打包版
```

### 12.4 结果文件对照

| 文件（`D:\MMModel\化合物库\`） | 内容 | 行数 |
|---|---|---|
| `compounds_tagged.csv` | 全库分子 + 全部理化指标 + 靶点标签 | 2,931 |
| `compounds.csv` / `library.db` | 主库与结构化索引 | 2,931 |
| `rejected.csv` | 被过滤分子及原因 | — |
| `top_candidates.csv` | 对接后候选 | — |
| `top50_novelty.csv` | Top50 及新颖性分类 | 50 |
| `selectivity_full.csv` | 全库 A2A/A1 双打分与 SI | 2,931 |
| `admet_screen.csv` | 全库 ADMET 规则分层 | 2,931 |
| `pareto_front.csv` | 五目标非支配解 | 91 |
| `pose_consistency.csv` | 对接模式一致性 | 2,908 |
| `roc_decoys.csv` | ROC 曲线点（100 阈值） | 100 |
| `final_recommendation.csv` | 六维推荐清单 | 200 |
| `top50_assayranked.csv` | 可测性重排后的 Top50 | 50 |
| `strain_energy.csv` | 构象应变能 | — |

### 12.5 术语表

| 术语 | 含义 |
|---|---|
| **口袋（pocket）** | 蛋白上可结合小分子的空腔，由中心坐标 + 盒子边长定义 |
| **束搜索（beam search）** | 每步保留 top-k 候选的启发式搜索，k = `beam_size` |
| **focal（前沿原子）** | 下一步生长要挂靠的已有原子 |
| **MDN** | 混合密度网络，用多个高斯分量建模多模态连续分布 |
| **QED** | 定量估计类药性（0–1，越大越像药） |
| **SA score** | 合成可及性评分（1 易合成 – 10 极难） |
| **Murcko 骨架** | 去掉侧链后的分子核心结构，用于衡量化学多样性 |
| **Tanimoto** | 分子指纹相似度（0–1） |
| **MMR** | 最大边际相关，兼顾相关性与多样性的选择策略 |
| **Pareto 非支配解** | 在所有目标上都不劣于其他解、且至少一个目标更优的解 |
| **SI（选择性指数）** | 本项目中为 A2A 与 A1 对接分百分位之差 |
| **富集倍数** | 在给定阈值下，命中已知活性分子的比例相对随机的倍数 |
| **conda-pack** | 把 conda 环境打包成可迁移归档的工具 |
| **onedir / onefile** | PyInstaller 的两种产物形态：目录版（免解压）/ 单文件版（启动时解压） |

---

**文档结束** · 本文与仓库代码同步于提交 `ace0596` 之后的工作区状态；
如需更新，请同步修改 `docs/` 下的说明与 `gui/README.md`。

