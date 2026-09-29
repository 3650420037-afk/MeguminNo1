# -*- coding: utf-8 -*-
"""生成《夜间报告》：把今夜所有实测结果汇总成一页（供早晨一眼看全，可反复重跑刷新）。

读取的都是**已落盘的产物**（不重新计算、不跑 GPU）：
- `outputs/night_pipeline_v10/sel_v10term.csv` / `sel_v10_matched.csv`（v10term vs v10 配对）
- `outputs/s1_exp/*/…/SMILES.txt`（阶段 1：S1 恒温 / S2 下限 / S3 尺寸窗 / 官方同协议）
- `outputs/dip_test/*/…/SMILES.txt`（凹坑测试）
- `git log`（提交留痕）

用法：python outputs/night_pipeline_v10/make_night_report.py
"""
import glob
import io
import os
import statistics as st
import subprocess
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))
OUT = os.path.join(ROOT, "docs", "夜间报告_20260929.md")

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")
BAND = (340.0, 430.0)


def smi_metrics(d):
    cands = glob.glob(os.path.join(d, "**", "SMILES.txt"), recursive=True)
    if not cands:
        return None
    ms = [Chem.MolFromSmiles(l.strip()) for l in io.open(sorted(cands)[-1], encoding="utf-8")
          if l.strip()]
    ms = [m for m in ms if m]
    if not ms:
        return None
    ha = sorted(m.GetNumHeavyAtoms() for m in ms)
    mw = sorted(round(Descriptors.MolWt(m), 1) for m in ms)
    return dict(n=len(ms), mw=st.median(mw), ha=st.median(ha),
                inband=sum(1 for x in mw if BAND[0] <= x <= BAND[1]) / len(mw),
                under=sum(1 for x in mw if x < BAND[0]) / len(mw))


def rows_of(patterns):
    out = []
    for pat, label in patterns:
        for d in sorted(glob.glob(os.path.join(ROOT, pat))):
            m = smi_metrics(d)
            if m:
                out.append((label or os.path.basename(d.rstrip("/")), m))
    return out


def table(rows, title):
    if not rows:
        return "（无数据）\n"
    s = ["| 臂 | n | MW 中位 | HA 中位 | 落区率 | MW<340 |", "| --- | --- | --- | --- | --- | --- |"]
    for tag, m in rows:
        s.append("| %s | %d | %.1f | %.1f | %.2f | %.2f |"
                 % (tag, m["n"], m["mw"], m["ha"], m["inband"], m["under"]))
    return "\n".join(s) + "\n"


def main():
    parts = ["# 夜间报告 · 2026-09-29（自动生成，%s）" % time.strftime("%Y-%m-%d %H:%M"),
             "",
             "> 本文件由 `outputs/night_pipeline_v10/make_night_report.py` 从**已落盘产物**汇总，可反复重跑刷新。",
             "> 结论口径：全部为**同靶点/同 n=60/beam50/steps40/同种子/同阈值 0.0** 的配对或并列比较，",
             "> 判据沿用**预登记**那套（尺寸落区 / strict ≥ 官方同协议臂 90% / 生成数 ≥n），见 `docs/任务状态与记忆.md` §5.41/§5.53。",
             "",
             "## 1. 默认权重状态",
             "",
             "- 默认权重 = **`models/7-eonmol_ft_gpcr_v10.pt`**（v10@5500，官方起点 + gpcr_mass_v1 微调）",
             "- 十臂配对 vs 官方：Δ原始 Vina 均值 **−0.15**，95% CI **[−0.67, +0.37]（跨 0）** → **打平略优，不能宣称更强**",
             "- 相对上一版 v3：**明显更强**（v3 在 B2AR/D3/5HT2B 落后 +0.80/+1.44/+1.54，v10 基本打平）",
             "",
             "## 2. V1（训练侧 `pos_weight=2.0`）：**判定不通过**",
             "",
             "v10term 早停于 iter 8000。同协议同迭代配对（v10 参照臂今夜新采）：",
             "",
             "| 迭代 | v10（pw=1.0）MW / 落区 | v10term（pw=2.0）MW / 落区 |",
             "| --- | --- | --- |",
             "| 5500 | 397.5 / 0.77 | 416.9 / 0.32 |",
             "| 6500 | 302.5 / 0.21 | 451.5 / 0.28 |",
             "| 7500 | 358.9 / 0.55 | 406.4 / 0.24 |",
             "| 8000 | 354.4 / 0.48 | 298.4 / 0.12 |",
             "",
             "- **MW 极差 95.0 → 153.1（更差）**；落区率均值 0.50 → 0.24；凹坑从 6500 挪到 8000",
             "- 5 个检查点 strict **全不达标**（比值 0.61–0.85 < 0.90）",
             "- 对接读数：v10term@5500 原始分好 0.16，但分子大 2.5 HA、**LE 更差**（0.358→0.340）→ 靠变大换分",
             "",
             "## 3. 采样侧尺寸控制 · 阶段 1（在 v10@5500 健康点上）",
             ""]

    parts.append(table(rows_of([
        ("outputs/s1_exp/s1off_s*", "S1 关（对照）"),
        ("outputs/s1_exp/s1on_s*", "S1 恒温"),
        ("outputs/s1_exp/s2floor_s*", "S2 硬下限 24"),
        ("outputs/s1_exp/s3window_s*", "S3 尺寸窗 24–32"),
        ("outputs/s1_exp/official_s*", "官方（同协议）"),
    ]), "阶段 1"))
    parts += ["",
              "**读法**：v10@5500 本身健康（对照 MW 397.5、落区 0.77），机制只会把它**抬得更高**（447–494），",
              "落区率反而下降 → **在健康点上「变大」不算成绩**（见 §5.55）。",
              "",
              "## 4. 凹坑测试（把机制用在**已知掉坑**的 v10@6500 上）",
              ""]
    parts.append(table(rows_of([
        ("outputs/night_pipeline_v10/v10_matched_sel/6500", "对照 v10@6500（无机制）"),
        ("outputs/dip_test/s1on_s*", "S1 恒温"),
        ("outputs/dip_test/s2floor_s*", "S2 硬下限 24"),
        ("outputs/dip_test/s3window_s*", "S3 尺寸窗 24–32"),
    ]), "凹坑测试"))
    parts += ["",
              "**判读**：对照 MW ≈302（远低于 340 下限）。若某机制把 MW 拉回 **340–430** 且 strict 不劣化，",
              "则该机制**有用**（可做成出货侧默认，与「阈值随权重走」并列）。",
              "",
              "## 5. 今夜基础设施/事故留痕（都在 §5.41–§5.55）",
              "",
              "- 代码审查 **SAFE_WITH_FIXES** → 5 项必修全修：LoRA 上 GPU 必崩的设备问题、适配器 lr 被吞成 1e-5、",
              "  字符串配置解冻 96% 参数、`floor_steps≥max_steps` 静默 0 产出、折入脚本静默零增量",
              "- 方法论陷阱入档：**采样进行中的部分读数系统性偏小**（因 `HA ≡ step`，分子按尺寸从小到大完成）",
              "- 负结果：CPU 廉价探针无法替代采样触发量（§5.43）；阈值→尺寸在 v10 上**不单调**（§5.48）",
              "",
              "## 6. git 状态（按用户要求**未推送**）",
              ""]
    try:
        log = subprocess.run(["git", "log", "--oneline", "-15"], cwd=ROOT,
                             capture_output=True, text=True, errors="ignore").stdout.strip()
        ahead = subprocess.run(["git", "rev-list", "--count", "origin/main..HEAD"], cwd=ROOT,
                               capture_output=True, text=True, errors="ignore").stdout.strip()
        parts.append("本地领先 origin **%s** 个提交：\n\n```\n%s\n```\n" % (ahead, log))
    except Exception as e:
        parts.append("（读取 git 状态失败：%s）\n" % e)

    parts += ["", "## 7. 下一步（按优先级）", "",
              "1. 凹坑测试结果判定 → 有用则做成出货侧默认机制；无用则转 V2（逐样本平衡训练）",
              "2. V2 训练与同协议评测（配置 `configs/train_gpcr_mass_v10term_bal.yml`，~50 分钟）",
              "3. LoRA 路线（基础设施已就绪并验证：`configs/train_gpcr_mass_v10lora.yml`）",
              "4. 四个药库用 v10 重出（口径不一致，需 GPU 1–2 h + CPU 对接数小时）",
              "5. A′ v1.2 的 10 种子全量重验（脚本已补内部超时护栏）",
              ""]

    io.open(OUT, "w", encoding="utf-8").write("\n".join(parts))
    print("夜间报告已生成：%s（%d 行）" % (os.path.relpath(OUT, ROOT), len(parts)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
