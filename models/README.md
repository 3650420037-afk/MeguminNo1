# models/ · 最终模型权重

本目录存放**提交的最终模型权重**。

## 内容

| 文件 | 说明 |
|---|---|
| `pretrained_Pocket2Mol.pt` | Pocket2Mol 官方预训练权重，**42.8 MB**，本项目最终交付使用 |
| `README.md` | 本文件 |
| `.gitignore` | 权重例外规则（本目录的 `.pt` 需随仓库提交） |

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
