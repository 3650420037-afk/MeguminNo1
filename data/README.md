# 数据说明 · 来源、许可与划分

> 本目录**提交**小体积的说明、划分记录与参考数据；**不提交** GB 级 LMDB 训练张量与
> 原始结构下载。大数据集按下方「六、大数据集如何获取」重建。

---

## 一、目录内容

| 路径 | 提交 | 内容 |
|---|---|---|
| `README.md` | ✅ | 本文件：来源、许可、清洗、划分、泄漏防控 |
| `targets/` | ✅ | 5 个 GPCR 靶点的受体坐标与共晶配体定义（RCSB PDB，CC0） |
| `example/` | ✅ | 最小可跑示例：`4yhj.pdb` + 参考配体 `4yhj_reference_ligand.sdf` |
| `ligands/` | ✅ | 四靶点活性配体清单（由 ChEMBL 整理，CSV） |
| `known_drugs/` | ✅ | 评估用参考集：各靶点已知活性分子、诱饵集、选择性集合 |
| `druglib/` | ✅ | 已知药物小样（`ChEMBL_已知药物_100.csv`） |
| `ligand_data/` | ✅ | ChEMBL 原始拉取记录（JSON），用于溯源 |
| `gpcr_a2a/` `gpcr_b2ar/` `gpcr_d3/` `gpcr_5ht2b/` | ✅ | 各靶点数据集的 `index.pkl`、`split_by_name.pt`、`build_report.txt`、口袋 PDB |
| `gpcr_multitarget/` `gpcr_multitarget_v2/` | ✅ | 合并数据集及其划分记录 |
| `*.lmdb` | ❌ | 训练用的 LMDB 张量（GB 级），本地生成 |
| `gpcr_v3*/` `gpcr_pdb_raw/` `_pdb_cache/` | ❌ | 数据扩充实验的中间产物（该实验已回退，见第七节） |

---

## 二、数据来源与许可

| 数据 | 来源 | 版本/获取时间 | 许可 | 用途 |
|---|---|---|---|---|
| CrossDocked2020 处理后数据 | [3D-Generative-SBDD](https://github.com/luost26/3D-Generative-SBDD) 提供的 `crossdocked_pocket10.tar.gz` | 2026-08 获取 | 随 CrossDocked2020 原始许可（学术研究可用） | 原版 Pocket2Mol 预训练/训练 |
| 受体-配体复合物结构 | [RCSB PDB](https://www.rcsb.org/) | 2026-08 下载 | **CC0**（公有领域） | 靶点口袋定义、共晶配体参照 |
| 活性分子数据 | [ChEMBL](https://www.ebi.ac.uk/chembl/) | 2026-08 拉取 | **CC BY-SA 3.0**（需署名） | 构建靶点特异性数据集、已知活性对照集 |
| 开源预训练权重 | Pocket2Mol 官方 `pretrained_Pocket2Mol.pt` | 随仓库分发 | **MIT**（(c) 2022 Xingang Peng） | 生成模型的初始与最终权重 |

署名与第三方归属见仓库根 `NOTICE.md`。

### 使用的 5 个靶点结构

口袋中心 = **共晶配体重原子的实验质心**，由 `src/scripts/build_target_registry.py`
从结构文件自动提取（不手工录入），登记于 `configs/targets.json`。

| 靶点名 | PDB | 受体 | 共晶配体 | 口袋中心 (Å) | 沉积配体完整性 |
|---|---|---|---|---|---|
| A2A | 4EIY | 腺苷 A2A 受体（A2AAR-BRIL） | ZMA (ZM241385) | −0.42, 8.53, 17.13 | 完整（25/25 重原子） |
| B2AR | 2RH1 | β2 肾上腺素受体 | CAU (卡拉洛尔) | −29.52, 9.23, 6.94 | 完整（22/22） |
| D3 | 3PBL | 多巴胺 D3 受体 | ETQ (依替必利) | 0.09, −14.83, 10.43 | 完整（23/23） |
| 5HT2B | 4IB4 | 5-HT2B 受体（5-HT2B-BRIL） | ERM (麦角胺) | 22.45, 18.28, 11.73 | 完整（43/43） |
| A1AR | 5UEN | 腺苷 **A1** 受体（A1AR-bRIL） | DU1 (DU172) | 55.97, 58.90, 143.62 | **不完整**（沉积 35，CCD 定义 36，缺磺酰氟上的 F） |

> A1AR 的共晶配体在沉积结构中漏建了 1 个氟原子，因此该条目无法同时满足
> "使用实验坐标"与"使用完整分子"。建库脚本对这类条目**默认丢弃并在报告中留痕**
> （详见 `docs/数据建库穿模修复_验证报告.md`）。

---

## 三、自建数据集与划分

各靶点数据集由「受体结构 + ChEMBL 活性分子」构成：以共晶配体质心为口袋中心，
按 12 Å 提取口袋残基，配体取**实验坐标**。划分按 `split_by_name.pt` 记录。

| 数据集 | 来源 | 候选 SMILES | 有效样本 | 训练 | 验证 | 跳过 |
|---|---|---|---|---|---|---|
| `gpcr_a2a` | 4EIY + A2A 活性集 | 175 | 141 | 112 | 29 | 34 |
| `gpcr_b2ar` | 2RH1 + B2AR 活性集 | 191 | 179 | 143 | 36 | 12 |
| `gpcr_d3` | 3PBL + D3 活性集 | 99 | 95 | 76 | 19 | 4 |
| `gpcr_5ht2b` | 4IB4 + 5HT2B 活性集 | 81 | 73 | 58 | 15 | 8 |
| `gpcr_multitarget` | a2a + b2ar | — | 320 | 256 | 64 | — |
| `gpcr_multitarget_v2` | a2a(136) + b2ar(179) + d3(91) + 5ht2b(73) | — | 479 | 383 | 96 | — |

数字取自各数据集目录内的 `build_report.txt`（构建时的原始记录，未事后修改）。

划分方式：**按样本随机划分并固定种子**（`seed=2021`），训练/验证互斥；
各靶点数据集独立划分，合并数据集（`gpcr_multitarget*`）在合并后重新划分，
因此合并集的验证集与单靶点验证集**不保证互斥**，评测一律以 `gpcr_multitarget_v2`
的划分为准。

---

## 四、预处理与清洗记录

1. **配体坐标取实验坐标**（而非用 SMILES 重新生成构象后平移质心）。
   早期实现只平移质心，扭转角与实验构象无关，导致配体插进蛋白壁：
   全量诊断 16.0 % 样本最近重原子距离 < 0.5 Å（实验位姿为 2.5–3.3 Å）。
   修复后同口径下可用样本穿模为 0。详见 `docs/数据建库穿模修复_验证报告.md`。
2. **元素过滤**：词表外元素直接跳过。实测跳过的为 Br（原子序数 35）、I（53）、
   Na（11），占四靶点候选的 4 %–19 %（见各 `build_report.txt` 的 `skipped_*` 行）。
3. **掩码配置缺陷修复**：原 `max_ratio: 1.1` 且 `min_num_unmasked: 0` 会让配体以
   约 9 % 概率被**全部掩码**，使 context 为空并产生 `NaN` 损失。已在 6 个配置中把
   `min_num_unmasked` 改为 1，并在 `src/models/maskfill.py` 的 `get_loss` 中加边界守卫。
4. **建库可复现**：原 `base_index()` 使用内置 `hash(str)`，其进程随机盐会让同一份
   输入每次建出不同数据集；已改为 `zlib.crc32`。

---

## 五、数据泄漏防控

- **训练/验证按样本互斥**，划分文件随数据集提交，可复核。
- **已知活性对照集只用于评估**（新颖性、富集倍数、ROC），**不进入任何训练集**。
- **靶点结构来自实验晶体结构**，不使用生成模型产出的结构做训练标签。
- **去重**：入库按规范化 SMILES 去重；同一分子不重复计入统计。
- 本项目为**生成式**任务（从口袋生成新分子），与数据集内的活性分子不构成
  "记忆-检索"式泄漏；但为稳妥，`predict.py` 会报告每个候选与已知活性的
  最大 Tanimoto 相似度并在 `备注` 中标注 ≥ 0.8 的高相似候选。

---

## 六、大数据集如何获取

`*.lmdb`（GB 级）不在仓库内，按下列步骤重建：

```bash
# 1) 原版预训练数据(CrossDocked2020 处理版)
#    按本文件顶部说明从 3D-Generative-SBDD 的发布链接下载
#    crossdocked_pocket10.tar.gz 与 split_by_name.pt, 解压到 data/

# 2) 靶点结构自带于 data/targets/, 无需下载; 如需扩充可从 RCSB 重新获取
#    只保留 PF00001(人源 GPCR)且有类药配体的条目
python src/scripts/fetch_gpcr_structures.py --limit 800 --workers 8
#    再用抓下来的结构重建"实验坐标 + 口袋提取"数据集
python src/scripts/build_gpcr_v3.py --src data/gpcr_pdb_raw --out data/gpcr_v3 --workers 6

# 3) 合并为多靶点训练集(在合并后重新划分, 默认验证集比例 0.2)
python src/scripts/merge_gpcr_datasets.py \
    --inputs data/gpcr_a2a data/gpcr_b2ar data/gpcr_d3 data/gpcr_5ht2b \
    --output data/gpcr_multitarget_v2 --val-ratio 0.2

# 4) 训练脚本读到 .lmdb 不存在时会自动由 pocket/ligand 文件生成
python train.py --config configs/train_multitarget_v2.yml
```

数据质量复核（构建后可运行）：

```bash
python src/scripts/diagnose_dataset.py --data data/gpcr_multitarget_v2 --mask-samples 0
```

---

## 七、关于历史记录中的旧路径

各 `build_report.txt`、`merge_log.txt` 等**构建时生成的记录**里出现的
`D:\MMModel\配体数据\...`、`D:\MMModel\靶点结构\...`、`D:\MMModel\Pocket2Mol\data\...`
是**重构前的目录**，对应现在的 `data/ligands/`、`data/targets/`、`data/`。
这些记录为保持原始性未作修改。

此外，`data/` 根目录残留的 `gpcr_v3_clean_*`、`gpcr_v3_merged_*` 文件属于
**已回退的数据扩充实验**（该实验的微调未通过 A/B 对照，已整体回退，
详见 `docs/7-eonmol_技术全书.md`），保留仅为可追溯，不参与当前训练与评测。
