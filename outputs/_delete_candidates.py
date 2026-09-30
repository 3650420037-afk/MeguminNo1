# -*- coding: utf-8 -*-
"""列出**待删候选**（只读，不删除）：按要求分四类，给出体积/文件数与"是否证据"标注。

用法：python outputs/_delete_candidates.py
"""
import glob
import io
import os
import subprocess
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def size_of(paths):
    """返回 (apparent_bytes, files)。apparent = 逻辑大小；lmdb 可能是稀疏文件。"""
    tot = n = 0
    for p in paths:
        if os.path.isfile(p):
            tot += os.path.getsize(p); n += 1
        elif os.path.isdir(p):
            for dp, _dn, fns in os.walk(p):
                for f in fns:
                    try:
                        tot += os.path.getsize(os.path.join(dp, f)); n += 1
                    except OSError:
                        pass
    return tot, n


def allocated(paths):
    """用 du 取**实际占用**（Windows 上稀疏文件会小于逻辑大小）。"""
    args = ["du", "-sc", "--block-size=1"] + [p for p in paths if os.path.exists(p)]
    if len(args) <= 4:
        return None
    try:
        out = subprocess.run(args, capture_output=True, text=True, errors="ignore").stdout
        return int(out.strip().split("\n")[-1].split("\t")[0])
    except Exception:
        return None


def gb(b):
    return "%.2f GB" % (b / 1e9) if b else "-"


def show(title, paths, note):
    paths = [p for p in paths if os.path.exists(p)]
    if not paths:
        print("## %s\n（无）\n" % title)
        return 0
    ap, n = size_of(paths)
    al = allocated(paths)
    print("## %s" % title)
    print("- 条目数 **%d**，文件 **%d**，逻辑大小 **%s**，实际占用 **%s**" % (len(paths), n, gb(ap), gb(al) if al else "?"))
    print("- 性质：%s" % note)
    for p in sorted(paths)[:12]:
        one, _ = size_of([p])
        print("  - `%s` — %s" % (os.path.relpath(p, ROOT), gb(one)))
    if len(paths) > 12:
        print("  - …另有 %d 项" % (len(paths) - 12))
    print()
    return ap


def main():
    print("# 待删候选清单（**未删除任何文件**）\n")
    total = 0

    # A. B2AR 相关（用户明确要求清除）
    a = ["results/B2AR"] + \
        glob.glob(os.path.join(ROOT, "data", "gpcr_b2ar*")) + \
        glob.glob(os.path.join(ROOT, "other", "library_backups", "*", "B2AR")) + \
        glob.glob(os.path.join(ROOT, "outputs", "*B2AR*")) + \
        glob.glob(os.path.join(ROOT, "outputs", "*b2ar*"))
    total += show("A. B2AR 相关（生成入口已禁用 → 库与中间产物可清）", a,
                  "B2AR 药库/中间产物；**回归证据（CSV/日志）建议保留**，见 D")

    # B. 历史药库
    b = ["results/official_weight_library"] + \
        glob.glob(os.path.join(ROOT, "other", "library_backups", "*")) + \
        glob.glob(os.path.join(ROOT, "other", "outputs_history", "*official*")) + \
        glob.glob(os.path.join(ROOT, "other", "outputs_history", "*library*"))
    total += show("B. 历史药库（官方权重库 5,864 分子 + v2 时代四靶点库备份）", b,
                  "**交付历史与对照证据**；删除后 `official_weight_library` 的对账数字无法本地复算")

    # C. 纯中间文件
    c = glob.glob(os.path.join(ROOT, "outputs", "**", "samples_*.pt"), recursive=True) + \
        glob.glob(os.path.join(ROOT, "outputs", "mass_dock")) + \
        glob.glob(os.path.join(ROOT, "outputs", "mass_dock_build")) + \
        glob.glob(os.path.join(ROOT, "logs", "*", "checkpoints")) + \
        glob.glob(os.path.join(ROOT, "**", "__pycache__"), recursive=True)
    total += show("C. 纯中间文件（采样快照 / 对接中间 / 训练检查点 / 缓存）", c,
                  "**可重建或已无用**：`samples_*.pt` 是无下游读取的束状态快照；对接 pdbqt 可重跑；"
                  "`__pycache__` 自动重建；训练检查点除已抽取的 v10@5500 外均非交付物（但为 §5.33 检查点分析的证据）")

    # D. 需你判断的巨型归档/缓存
    d = [os.path.join(ROOT, "data", "%s_processed.lmdb" % n) for n in
         ("gpcr_b2ar", "gpcr_dock_v1", "gpcr_ft_v3", "gpcr_ft_v3c", "gpcr_mass_v1",
          "gpcr_multi_v1", "gpcr_multitarget", "gpcr_multitarget_v2", "gpcr_v3_clean", "gpcr_v3_merged")] + \
        [os.path.join(ROOT, "other", "outputs_history")]
    total += show("D. 巨型缓存/归档（收益最大，但需你确认）", d,
                  "`*_processed.lmdb` = **训练时自动生成的数据集缓存**（lmdb map_size 固定 10 GB/个），"
                  "删除后下次用该数据集会自动重建（需原始数据在库，均在）；`other/outputs_history` = 历史 outputs 归档")

    print("本清单涉及逻辑总量约 **%s**（实际占用更少，lmdb 多为稀疏分配）。" % gb(total))
    print("\n**等你确认后我再删除**；建议顺序：C（纯中间）→ A（B2AR）→ D（缓存/归档）→ B（历史药库，最后）。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
