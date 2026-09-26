# results/ · 最终结果与候选清单

本目录包含**两代药物库**：当前交付库（由本项目**自训模型**生成）与历史库
（由**官方权重**生成，完整归档保留，未删除）。

---

## 一、当前交付库（自训模型生成）

**生成模型**：`models/7-eonmol_ft_gpcr_v2.pt`（本项目自训，经种子 2024/2025 双轮 A/B 晋级）
**靶点**：A2A 腺苷受体（4EIY）｜**采样**：束宽 50 / 上限 60 步 / **12 个批次**（种子 2024–2042）

| 阶段 | 数量 |
|---|---|
| 生成完成的分子 | 1,857 |
| 过滤入库（PAINS+BRENK、MW[250,500]、LogP[1,5]、QED≥0.4、SA≤6、元素、环数） | **521** |
| 对接成功（AutoDock Vina） | 521 / 521 |
| 最终候选清单（按 Vina 结合分排序取前 300） | **300** |

### 关键指标

| 指标 | 全库（521） | 候选清单（300） |
|---|---|---|
| Vina 结合分 中位 | **−9.70 kcal/mol** | **−10.53** |
| 强于 −10 kcal/mol | 220（**42.2%**） | 220（**73.3%**） |
| 最强结合分 | **−14.45** | −14.45 |
| QED 中位 | 0.757 | 0.744 |
| SA 中位（越低越易合成） | 3.224 | 3.290 |
| MW 中位 | 390.5 | 400.5 |
| 与已知活性最大 Tanimoto 中位（越低越新） | 0.281 | 0.281 |
| 唯一 Murcko 骨架 | **373 / 521** | — |

### 文件

| 文件 | 说明 |
|---|---|
| `results.csv` | **最终候选清单**（300 个，34 列，含对接分/QED/SA/新颖性/理化性质/模型版本/种子） |
| `structures/` | 与清单「结构文件」列一一对应的 SDF（300 个） |
| `run_manifest.json` | 参数、种子、权重路径、环境、耗时（可溯源） |
| `compounds.csv` | 入库全量分子与描述符（521 个） |
| `rejected.csv` | 被过滤分子的淘汰原因（1,336 个） |
| `docking_summary.csv` | 全库对接结果（smiles / vina_kcal / 受体 / 盒子） |
| `top_candidates.csv` | 按对接分排序的全库榜单 |
| `sdf/` | 入库分子的三维构象 |
| `library.db` | SQLite 库 |

**复现命令**

```bash
# 1) 生成（12 批次，种子 2024–2042）
python design.py --target A2A --num-samples 200 --beam 50 --max-steps 60 --seed <SEED> \
    --outdir outputs/library_v2
# 2) 入库 + 打分 + 出清单
python predict.py --skip-design --target A2A --top 300 \
    --runs outputs/library_v2/<session>/sample_for_pdb_* \
    --outdir outputs/predict_v2 --results-dir outputs/results_v2
# 3) 对接（AutoDock Vina）
python src/scripts/docking_pipeline.py --receptor data/targets/4EIY_A2A受体.pdb \
    --center=-0.42,8.53,17.13 --smiles-file <SMILES> --out <OUT> --exhaustiveness 4
# 4) 带对接分重出清单
python predict.py --library <库目录> --target A2A --top 300 --rank-by vina
```

---

## 二、历史库（官方权重生成，已归档）

路径：**`official_weight_library/`**（9,179 文件，完整保留）

| 项 | 内容 |
|---|---|
| 生成模型 | `models/pretrained_Pocket2Mol.pt`（ICML 2022 官方权重，未微调） |
| 来源 | 73 个历史批次（`batch_*`） |
| 规模 | 原始 4,231 → 入库 **2,931** |
| Vina | 均值 **−9.68**，1,137 个强于 −10（38.8%），最强 −13.51 |
| 附带分析 | 对接、ADMET 分级、选择性、Pareto 前沿、ROC 富集、姿态一致性、六维 Tier 推荐（`final_recommendation.csv`，Tier1 166 个） |

**为什么保留**：这一整套下游分析（对接/ADMET/选择性/Pareto/ROC/六维推荐）是在该库上完成的，
与当前交付库不可互换；且它是本项目的完整评估证据链，删除会丢失可追溯性。

**两库对比（同一靶点、同一对接流程）**

| | 历史库（官方权重） | **当前库（自训模型）** |
|---|---|---|
| 入库分子 | 2,931 | 521 |
| Vina 中位 / 均值 | 均值 −9.68 | 中位 **−9.70** |
| 强于 −10 的比例 | 38.8% | **42.2%**（候选清单 73.3%） |
| 最强结合分 | −13.51 | **−14.45** |
| 下游六维分析 | 完整 | 未做（需另行运行 `src/scripts/` 下的分析链） |

> **如实说明**：当前库的规模远小于历史库（521 vs 2,931），因为历史库来自 73 个批次而当前
> 为 12 个批次；但**单位分子的对接表现更好**（强于 −10 的比例 42.2% vs 38.8%，最强分更低）。
> 当前库尚未跑 ADMET/选择性/Pareto/六维推荐，因此 `results.csv` 中这些列为空
> —— 这是**没有做**，不是"结果为零"。
