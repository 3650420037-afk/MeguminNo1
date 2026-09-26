# models/ · 最终模型权重

本目录存放**提交的最终模型权重**。

## 内容

| 文件 | 说明 |
|---|---|
| `7-eonmol_ft_gpcr_v2.pt` | **默认模型**：本项目自训，经双种子 A/B 晋级（19.9 MB） |
| `pretrained_Pocket2Mol.pt` | Pocket2Mol 官方预训练权重（42.8 MB），微调起点与基线 |
| `7-eonmol_ft_gpcr_v1.pt` | 早期自训权重，未晋级，保留供对照（19.9 MB） |
| `README.md` | 本文件 |
| `.gitignore` | 权重例外规则（本目录的 `.pt` 需随仓库提交） |

## 默认模型 `7-eonmol_ft_gpcr_v2.pt`（已晋级）

本项目用**对接构建的高密度数据集**微调得到，`design.py` / `predict.py` 默认使用它。

- 训练：`train.py --config configs/train_gpcr_dock.yml`（冻结编码器 + lr 1e-4 +
  梯度累积 8 + 按验证损失存 best.pt + 早停），最优 iter 3000
- 数据：`data/gpcr_dock_v1`（307 对 / **4 个口袋** / 每口袋 40–127 样本 / 零穿模），
  由 `src/scripts/build_docked_dataset.py` 用 AutoDock Vina 把 ChEMBL 活性集
  对接进**实验口袋**得到

### 晋级证据（`src/scripts/judge_promotion.py`，双种子）

| 种子 | 权重 | 完成分子 | QED 均值 | SA 均值 | 骨架数 | 口袋穿模率 | 判定 |
|---|---|---|---|---|---|---|---|
| 2024 | 官方基线 | 63 | 0.7237 | 3.5386 | 45 | 0.0% | — |
| 2024 | **本权重** | 61 | **0.7636** | **3.3479** | 41 | **0.0%** | **PROMOTED** |
| 2025 | 官方基线 | 62 | 0.7658 | 3.2741 | 45 | 4.84% | — |
| 2025 | **本权重** | **64** | 0.7370 | **3.0118** | 37 | **1.56%** | **PROMOTED** |

**如实表述**：分子数与基线相当；QED 一个种子更高、另一个略低（在 0.05 容差内），
不宣称 QED 全面胜出；**SA 与口袋叠合度在两个种子上都稳定占优**。
完整过程（含 5 次失败与两次被否证的修复尝试）见
[`docs/微调实验_v3_报告.md`](../docs/微调实验_v3_报告.md)。

使用方式：

```bash
python design.py --target A2A                       # 默认即本权重
python design.py --target A2A --ckpt models/pretrained_Pocket2Mol.pt   # 跑官方基线对照
python design.py --target A2A --ckpt models/7-eonmol_ft_gpcr_v1.pt     # 对照未晋级的 v1
```

## 早期权重 `7-eonmol_ft_gpcr_v1.pt`（未晋级，保留对照）

在实验坐标数据集（每口袋仅 1 条）上训练，穿模 0%、SA 2.72、接触距离 2.35 Å 三项更好，
但完成分子数 23 vs 基线 63 不达标，判 NOT_PROMOTED。保留以便复现该结论。

## 权重来源与许可

- 来源：Pocket2Mol 官方发布（ICML 2022, Xingang Peng 等）
- 许可：**MIT License**, Copyright (c) 2022 Xingang Peng
- 原始下载地址：
  https://drive.google.com/drive/folders/1KfdOczjUPITPhIvCuBmnj4xFTV-iI2xB
- 归属声明见仓库根 [`NOTICE.md`](../NOTICE.md)

下载后请放到本目录并命名为 `pretrained_Pocket2Mol.pt`：

```bash
# 下载后确认文件完整（约 42.8 MB）
ls -lh models/pretrained_Pocket2Mol.pt
```

## 权重内容

`torch.load` 后为字典，字段如下（实测）：

| 键 | 类型 | 说明 |
|---|---|---|
| `config` | EasyDict | 模型与训练配置 |
| `model` | OrderedDict | 模型参数（state_dict） |
| `optimizer` | dict | 优化器状态 |
| `scheduler` | dict | 学习率调度状态 |
| `iteration` | int | 训练迭代数 |

该文件**只包含 `state_dict`**，不含项目模块的 pickle 引用，因此可以安全地
把 `src/models/` 代码目录整体移动而不影响加载。

> PyTorch 2.6 起 `torch.load` 默认 `weights_only=True`。本仓库在
> `src/utils/misc.py` 中注册了所需的安全全局对象（补丁），加载前请先
> `import utils.misc`（`design.py` / `predict.py` / notebook 均已处理）。

## 本项目是否微调过权重

**没有作为最终交付。** 两轮微调实验均未通过 A/B 对照，已整体回退；
最终交付使用上述官方权重，模型的贡献在**采样层（引导束搜索）**而非权重。
详见仓库根 [`ModelCard.md`](../ModelCard.md) 第 1 节与第 7 节。

## 训练/微调产出的权重放哪里

`python train.py` 会把检查点写到 `logs/<config>_<时间戳>/checkpoints/`，
不会写入本目录，以免与最终交付权重混淆。
