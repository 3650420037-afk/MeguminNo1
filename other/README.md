# other/ · 非提交内容

> 本目录**内容不随仓库提交**（仅本说明与迁移记录 `_restructure_log.json` 除外，
> 以便评测方了解哪些内容被排除、为什么），也不属于大赛提交物。
> 放在这里是为了：① 保持提交目录结构整洁；② 不删除证据，便于随时回溯。
>
> 因此：**克隆仓库后不会看到下表中的大体积目录**，它们只存在于作者的本地工作副本中；
> 需要用到的数据请按 `data/README.md` 第六节自行重建。

---

## 目录内容

| 路径 | 体积 | 内容 | 为什么放这里 |
|---|---|---|---|
| `outputs_history/` | ~180 GB | 历代采样输出（四靶点 9 万余个分子文件与过程快照） | 运行产物，非提交内容；当前结果已汇总进 `results/` |
| `backup_original_20260826/` | ~50 MB | 初始克隆的原始 Pocket2Mol 仓库快照 | 历史备份，仅用于比对原版 |
| `backup_before_local_finetune_20260826/` | ~40 MB | 本地微调实验前的快照 | 同上 |
| `target_structures_raw/` | ~10 MB | 靶点结构原始下载（10 个文件） | 与提交的 `data/targets/`（已补全为 10 个文件）重复，保留原始副本 |
| `ligand_data_duplicates/` | ~40 KB | 与 `data/ligands/` **逐字节相同**的 4 个 CSV | 重复文件 |
| `gui_build_artifacts/` | ~170 MB | PyInstaller 的 `build/`、`dist/`（含打包好的 exe） | 构建产物，可用 `src/gui/build_exe.ps1` 重现 |
| `训练文件清单.txt` | 9 KB | 与 `docs/训练文件清单.txt` 重复 | 重复文件 |
| `run_sample.bat` | 1 KB | 根目录散落的旧启动脚本 | 已被 `design.py` / `predict.py` 取代 |
| `整理配体数据.py` | 3.5 KB | 与 `src/scripts/整理配体数据.py` 重复 | 重复文件 |
| `训练数据/` | 0 | 空目录 | 无内容 |
| `_restructure_log.json` | — | 本次重构的完整迁移记录（52 项移动） | 变更留痕 |

---

## 重构前后对照

本次按大赛附件5《代码提交要求》的目录模板重组，**仓库根由 `D:\MMModel\Pocket2Mol`
上移为 `D:\MMModel`**。主要对应关系：

| 重构前 | 重构后 |
|---|---|
| `Pocket2Mol/models/`（架构代码） | `src/models/` |
| `Pocket2Mol/utils/` | `src/utils/` |
| `Pocket2Mol/scripts/` | `src/scripts/` |
| `Pocket2Mol/evaluation/` | `src/evaluation/` |
| `Pocket2Mol/gui/` | `src/gui/` |
| `Pocket2Mol/sample.py`、`sample_for_pdb.py` | `src/`（并新增根目录入口 `design.py` / `predict.py`） |
| `Pocket2Mol/ckpt/`（权重） | `models/` |
| `Pocket2Mol/configs/` | `configs/` |
| `Pocket2Mol/data/` | `data/` |
| `Pocket2Mol/targets/`、`ligands/`、`druglib/`、`example/` | `data/targets/`、`data/ligands/`、`data/druglib/`、`data/example/` |
| `Pocket2Mol/logs/` | `logs/` |
| `Pocket2Mol/outputs/` | `outputs/`（运行时）；历史产物 → `other/outputs_history/` |
| `Pocket2Mol/docs/`、`assets/` | `docs/`、`docs/assets/` |
| `Pocket2Mol/train.py` | `train.py`（适配新路径） |
| `D:\MMModel\化合物库\` | `results/` |
| `D:\MMModel\已知药物库\` | `data/known_drugs/` |
| `D:\MMModel\配体数据\*.json` | `data/ligand_data/`（其余 CSV 为重复） |
| `D:\MMModel\靶点结构\` | `data/targets/`（原始副本留在 `other/target_structures_raw/`） |

**未移动**：`tools/`（AutoDock Vina 与 Open Babel 二进制）仍留在仓库根，
它是分子对接这一步的真实依赖，不属于"多余文件"；已在 `.gitignore` 中排除，
并在 README 环境说明中标注。

完整的逐项迁移记录见 `_restructure_log.json`。

---

## 想彻底清理时

确认不再需要回溯后，可直接删除本目录下任意子目录；它们不参与任何构建或运行，
删除后执行 `python predict.py` 仍可正常工作。
唯一例外是 `target_structures_raw/`：若你还需要原始下载副本请先确认
`data/targets/` 完整。
