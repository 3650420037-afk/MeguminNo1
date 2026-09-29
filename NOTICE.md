# NOTICE — 来源、许可与第三方组件

本软件（7-eonmol · GPCR 靶向分子生成器）由以下部分构成。分发或商用前请阅读本文。

---

## 一、本软件构成与来源

| 部分 | 归属 | 许可 |
|---|---|---|
| 基础模型与框架 Pocket2Mol | Xingang Peng 等（北京大学），ICML 2022 | **MIT License**，Copyright (c) 2022 Xingang Peng |
| 本项目的改进与新增（引导束搜索、构象精修、过滤入库、对接管线、桌面界面与打包脚本） | 本项目作者 | 归本项目作者所有 |
| 依赖库（PyTorch / RDKit / PyG 等） | 各自作者 | 见第三节 |

**Pocket2Mol 原始许可（MIT）要求：在软件的所有副本或实质部分中，必须保留上述版权声明与本许可声明。**
因此分发本软件时，`LICENSE` 文件（Pocket2Mol 的 MIT 原文）必须随附，不得删除。

### 引用要求（学术场合）
使用本项目或其改进的引导束搜索方法发表工作时，请引用原论文：

> Peng, X., Luo, S., Guan, J., Xie, Q., Huang, J., Ma, J. (2022).
> Pocket2Mol: Efficient Molecular Sampling Based on 3D Protein Pockets. ICML 2022.

---

## 二、本项目相对基础模型做了什么（改进说明）

基础权重与网络结构未替换；改进集中在**采样与后处理环节**：

1. **引导束搜索**（`utils/guidance.py`）：束排序概率乘以化学分数增益
   `P ∝ (exp(ΣlogP) + 1) · exp(λ·[QED + (1 − SA/10) − w_div·maxTanimoto])`，
   实测 QED 中位 0.667 → 0.803（+20%，λ=3，双随机种子复现）。
2. **多样性惩罚**：对与已生成分子相似的结构施加 MMR 式惩罚，防止搜索塌缩。
3. **构象精修与几何审计**（`utils/reconstruct.py` 等）。
4. **全自动下游管线**：过滤入库、AutoDock Vina 对接、新颖性/选择性/ADMET/Pareto 多维评估。
5. **桌面应用**：面向非编程使用者的图形界面与免安装打包。

多靶点微调实验（479 对数据集）已如实记录为**未晋级**；而后续在对接构建的高密度数据集上
微调得到的 `models/7-eonmol_ft_gpcr_v2.pt` 最初判为晋级，但该结论已在 2026-09-27 被推翻
（原 A/B 用的 `diversity_w = 0` 并非出货目标；在出货目标 `0.5` 下复验，QED 两个种子都更低）。
详见 `ModelCard.md` 与 `docs/任务状态与记忆.md`。

**当前默认权重（2026-09-29 起）= `models/7-eonmol_ft_gpcr_v10.pt`**
（官方权重起点 + 自建 `gpcr_mass_v1` 12,244 对微调；十臂同种子成对 vs 官方：原始对接分
4 更好 / 4 打平 / 2 更差、均值 −0.15，生成尺寸与官方一致。**如实限定**：按项目预登记判据
v10 本判 NOT_PROMOTED，本次替换默认权重是用户显式决策，**不构成"全面达标/全面优于官方"**）。
`models/7-eonmol_ft_gpcr_v3.pt` 转为**配体效率/QED 为先的备选**（推荐阈值 −0.5）。
官方预训练权重（`models/pretrained_Pocket2Mol.pt`）始终保留为基线与回退。

---

## 三、第三方组件与许可（随包分发部分）

随包 `env\` 分发的 Python 环境含 61 个依赖包，主要许可如下（均为允许商用的许可）：

| 组件 | 许可 | 商用 |
|---|---|---|
| PyTorch / torchvision | BSD-3-Clause | ✅ |
| RDKit | BSD-3-Clause | ✅ |
| NumPy / SciPy | BSD-3-Clause | ✅ |
| PyTorch Geometric（含 cluster/scatter/sparse） | MIT | ✅ |
| PyYAML | MIT | ✅ |
| lmdb | OpenLDAP Public License 2.8 | ✅ |
| Pillow | MIT-CMU（HPND） | ✅ |
| easydict | LGPL-3.0 | ✅（动态使用，未修改其源码） |
| certifi / tqdm | MPL-2.0（tqdm 另含 MIT） | ✅ |
| PyInstaller（打包工具） | GPLv2-or-later **附例外条款**，明确允许打包并分发商业程序 | ✅ |
| Python 3.10 / Tcl-Tk | PSF-2.0 / TCL | ✅ |

> 说明：MPL-2.0 为文件级 copyleft（修改其源文件才需开源该文件）；
> LGPL-3.0 允许闭源商用（本项目以库形式调用，未修改其源码）。

**未随包分发的外部工具**（若将来一并分发，需另行遵守其许可）：
- Open Babel — **GPL-2.0**（强 copyleft；仅以独立进程调用通常视为聚合，但仍应谨慎评估）
- AutoDock Vina — Apache-2.0

---

## 四、数据来源

| 数据 | 来源 | 说明 |
|---|---|---|
| 靶点蛋白结构（4EIY / 2RH1 / 3PBL / 4IB4） | RCSB PDB | PDB 数据为公共领域（CC0），可自由使用 |
| 生成分子库与全部评估结果 | 本项目产出 | 归本项目作者 |
| 已知活性分子对照集（A2A 等） | ChEMBL | ChEMBL 数据通常为 CC BY-SA 类许可：**用于验证时需署名**；如对外分发该数据集，需遵守其"相同方式共享"条款 |
| 微调数据集（479 对，未晋级） | 由 PDB 复合物构建 | 未随本包分发 |
| 原始训练数据 CrossDocked2020 | 由 PDBbind 构建 | **未随本包分发**；PDBbind 对商业使用有限制，若需再训练请自行核实 |

---

## 五、商用前建议核实的三件事

1. **预训练权重**：`ckpt/README.md` 仅给出官方下载链接，未单独声明权重许可（通常随项目 MIT）。
   若要作为商业产品发布，建议向原作者书面确认权重的商用授权。
2. **训练数据上游**：若产品需要重新训练模型，须自行取得 CrossDocked2020 / PDBbind 的商用授权。
3. **医药用途合规**：本软件输出的是计算候选分子，不构成药品或临床建议；
   涉及药物研发的商业活动另有药监与伦理法规约束，与软件许可无关。

> 本文档为工程视角的许可梳理，**不构成法律意见**。正式商业化前建议由法务/律师复核。

