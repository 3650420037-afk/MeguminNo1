# -*- coding: utf-8 -*-
"""§5.63 追加：P29 时间线补写 + 本代 ADMET 回填 + 最终备份。"""
import io, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
DOC = r"D:\MMModel\docs\任务状态与记忆.md"
text = """
**补充（同日 16:05）：又发现两处并已处理**

1. **P29「项目历程」原是半成品**（只有标题 + 一条时间轴线图片，正文只有 46 个字符）——
   已按仓库真实历元补 7 条时间线（`fill_p29_timeline.py`，脚本可重跑）：
   08-26 基线打通 → 08-29 官方权重库建成（8,904→5,864）→ 09-21 库级评估体系（ROC 5.80×/ADMET/选择性/Pareto/六维）→
   09-27 高密度数据集 12,244 对 + 架构改造（终止模块 4.2×、编码器解冻）→ 09-29 v10 晋级默认权重（十臂 4/4/2，均值 −0.15）→
   09-29 引导束搜索定型（λ=3、QED +20%）→ 09-30 采样侧尺寸控制 + 药库重出（1,151→707）+ B2AR 入口禁用。
   （其余各页文本量 200–350 字符，无第二处半成品；用"逐页文本量统计"筛出来的。）
2. **候选表的 ADMET 列本代是空的**（`results/results.csv` 只重跑了对接），而附件3 §五写着"均含…ADMET 分级" ——
   已用本项目自己的规则脚本 `src/scripts/admet_screen.py --library <dir>` **对本代三靶点重跑**并回填
   （`fill_admet_columns.py`，按 SMILES 精确匹配 300/106/206 全部命中）：
   A2A 库 A 102 / B 237 / C 56（候选 300：A 77 / B 177 / C 46）、D3 候选 A 21 / B 70 / C 15、5HT2B 候选 A 27 / B 148 / C 31。
   上一代 `admet_screen.csv` 四个文件已备份到 `other/library_backups/20260930_admet_prev/`（不删只归档）。
   **仍未重跑（列保持空白，已知并留痕）**：选择性 SI 百分位、六维 Tier、库级 ROC/Pareto、姿态置信度 ——
   只有历史库（`results/official_weight_library/`，09-21 完成）有；已在 `results/README.md` 新增
   「本代库的重跑分析范围（如实说明）」与根 README 的成药性行注明。
3. **最终备份**：修订后的 PPT 存 `outputs/ppt_backup/gpcr_(1)_修订后_20260930_1600.pptx`；
   终版逐页渲染在 `outputs/ppt_png/final/`（32 页）。
"""
io.open(DOC, "a", encoding="utf-8", newline="\n").write(text)
print("追加 %d 字符" % len(text))
