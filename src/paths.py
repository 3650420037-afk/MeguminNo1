# -*- coding: utf-8 -*-
"""项目路径的唯一来源 —— 禁止在其它模块里硬编码盘符路径(附件5 第六节)。

用法(脚本位于 src/<sub>/x.py 时):

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from paths import DATA, RESULTS, OUTPUTS   # 用到什么引什么

所有路径都由本文件所在的 src/ 反推, 因此整棵目录可以整体搬移而无需改代码。
需要覆盖默认位置时用环境变量, 不要改源码。
"""
import json
import os
import sys

__all__ = [
    "ROOT", "SRC", "DATA", "MODELS", "CONFIGS", "RESULTS", "LOGS", "OUTPUTS",
    "DOCS", "TOOLS", "OTHER", "NOTEBOOKS",
    "TARGETS", "EXAMPLE", "LIGANDS", "KNOWN_DRUGS", "LIGAND_DATA", "DRUGLIB",
    "PRETRAINED", "VINA", "OBABEL", "PYTHON", "FINETUNED", "DEFAULT_CKPT",
    "src_on_path", "ensure_dir", "require", "missing_hint",
]

# ---------------------------------------------------------------- 根目录
# src/paths.py -> src/ -> 仓库根
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")

# ---------------------------------------------------------------- 模板要求的顶层目录
DATA = os.environ.get("EONMOL_DATA") or os.path.join(ROOT, "data")
MODELS = os.environ.get("EONMOL_MODELS") or os.path.join(ROOT, "models")
CONFIGS = os.path.join(ROOT, "configs")
RESULTS = os.path.join(ROOT, "results")
LOGS = os.path.join(ROOT, "logs")
OUTPUTS = os.path.join(ROOT, "outputs")
DOCS = os.path.join(ROOT, "docs")
NOTEBOOKS = os.path.join(ROOT, "notebooks")

# 第三方二进制与本地非提交内容
TOOLS = os.environ.get("EONMOL_TOOLS") or os.path.join(ROOT, "tools")
OTHER = os.path.join(ROOT, "other")

# ---------------------------------------------------------------- data/ 子目录
TARGETS = os.path.join(DATA, "targets")
EXAMPLE = os.path.join(DATA, "example")
LIGANDS = os.path.join(DATA, "ligands")
KNOWN_DRUGS = os.path.join(DATA, "known_drugs")
LIGAND_DATA = os.path.join(DATA, "ligand_data")
DRUGLIB = os.path.join(DATA, "druglib")

# ---------------------------------------------------------------- 权重
PRETRAINED = os.path.join(MODELS, "pretrained_Pocket2Mol.pt")

# 交付权重 v3 = 架构实验 v8 的 iter 800(解冻 encoder + 扩容 frontier_pred)。
# 晋级依据(全部在**出货目标** diversity_w=0.5 下、A2A、同种子成对、双种子):
#   - 结合强度(第一优化轴): **与官方基线打平**(Vina 中位差 -0.04 / -0.03, 远小于
#     同臂运行间噪声 ~0.14), 且 <=-10 占比两个种子都更高(36.5% vs 35.8%; 60.0% vs 55.7%);
#   - QED 中位两个种子都更高(+0.065 / +0.082)、SA 两个种子都更低(更好)、
#     骨架数打平或更多(47/47 与 55/39)、穿模率两个种子均为 0%;
#   - 完成数相当(64/63 与 60/61) -> 停止策略健康。
# 因此它是目前唯一在所有已测轴向上都不劣于官方权重、且多项更优的模型。
# ⚠ 局限: 对接验证仅 A2A 一个靶点、每臂 50-60 个分子; 不得宣称"结合强度更强"。
# ⚠ 架构随权重走: v3 的 frontier 容量为 256/64/2, 记录在权重文件自己的 config 里,
#   sample_for_pdb.py 正是按 ckpt['config'].model 构建模型, 故可直接加载。
#
# 历史权重(保留供对照, 勿删):
#   - 7-eonmol_ft_gpcr_v2.pt: 原判"晋级", 但在出货目标下**双种子结合强度与 QED 都更差**,
#     已判 NOT_PROMOTED(见 ModelCard.md); 保留供复现该结论。
#   - 7-eonmol_ft_gpcr_v1.pt: 早期权重, 未晋级。
FINETUNED = os.path.join(MODELS, "7-eonmol_ft_gpcr_v3.pt")
DEFAULT_CKPT = FINETUNED if os.path.exists(FINETUNED) else PRETRAINED

# ---------------------------------------------------------------- 外部程序
# 默认取仓库内 tools/; 可用环境变量指向别处(例如系统安装的 vina)
VINA = os.environ.get("EONMOL_VINA") or os.path.join(TOOLS, "vina", "vina.exe")
OBABEL = os.environ.get("EONMOL_OBABEL") or os.path.join(
    TOOLS, "openbabel", "Library", "bin", "obabel.exe")
OBABEL_DATA = os.environ.get("EONMOL_OBABEL_DATA") or os.path.join(
    TOOLS, "openbabel", "Library", "share", "openbabel")
# 当前解释器; GUI 需要以子进程方式跑采样时使用
PYTHON = os.environ.get("EONMOL_PYTHON") or sys.executable


def src_on_path():
    """把 src/ 放进 sys.path, 使 `import models / utils / scripts` 保持可用。"""
    if SRC not in sys.path:
        sys.path.insert(0, SRC)
    return SRC


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def require(path, what):
    """路径不存在时给出可操作的提示, 而不是抛出难懂的 IOError。"""
    if not os.path.exists(path):
        raise SystemExit(
            "[paths] 缺少%s:\n    %s\n请按 README 的“数据与权重获取”一节准备, "
            "或用环境变量覆盖默认位置。" % (what, path))
    return path


def missing_hint(paths_with_names):
    """返回缺失项的清单文本, 供入口脚本一次性提示。"""
    lines = []
    for name, p in paths_with_names:
        if not os.path.exists(p):
            lines.append("  - %s: %s" % (name, p))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 推荐采样阈值随权重走（2026-09-28）
# ---------------------------------------------------------------------------
# 动机（实测，见 docs/任务状态与记忆.md §5.16 / §5.18 / §5.21）：
#   * v3 在 frontier_threshold=-0.5 下，A2A 9 个种子中 8 个原始 Vina 更好
#     （符号检验单侧 p=0.0195，中位 -0.765），尺寸 20->25 重原子、strict 63.9%->78.2%；
#   * 但同一 -0.5 会伤害 v10（分子本来就到位）—— 即"最优阈值与权重配套"。
# 因此阈值不做全局硬编码，而是**跟着权重走**：权重自带推荐值时用它，否则查
# models/threshold_recommendations.json，再否则回退 0.0（= 原版行为，向后兼容）。
THRESHOLD_TABLE = os.path.join(MODELS, "threshold_recommendations.json")


def recommended_frontier_threshold(ckpt_path, default=0.0, peek_ckpt=True):
    """返回该检查点的推荐 frontier_threshold，以及来源说明。

    解析顺序：① 权重文件 config 里的 frontier_threshold（未来晋级的权重会内置）
              ② models/threshold_recommendations.json 侧车表（按文件名）
              ③ default（0.0，与官方/原版行为一致）
    peek_ckpt 为 True 时会尝试 torch.load（只读 config，失败则跳过）——15–45 MB 权重，
    仅在 design.py / 评测脚本启动时调用一次。
    """
    name = os.path.basename(ckpt_path)
    if peek_ckpt:
        try:
            import torch
            ck = torch.load(ckpt_path, map_location="cpu", weights_only=False)
            val = (ck.get("config") or {}).get("frontier_threshold")
            if val is not None:
                return float(val), "权重内置(frontier_threshold)"
        except Exception:
            pass
    try:
        with open(THRESHOLD_TABLE, encoding="utf-8") as f:
            tab = json.load(f)
        if name in tab:
            return float(tab[name]), "models/threshold_recommendations.json"
    except Exception:
        pass
    return float(default), "默认(原版行为)"
