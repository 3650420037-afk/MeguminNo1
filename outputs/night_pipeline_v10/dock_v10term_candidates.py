# -*- coding: utf-8 -*-
"""扫描后贯通：从 sel_v10term.csv 里挑**通过预登记判据**的候选，复用**已采分子**做对接 A/B。

为什么复用：`select_ckpt_by_generation.py --keep-dir` 已把每个检查点的采样产物（含 SMILES.txt）
留在盘上 → 对接级 A/B **不需要重新采样**（采样是 GPU 瓶颈，对接是 CPU）。
工具链：`arm_gen_metrics.py`（尺寸/strict 口径）→ `dock_smiles_ab.py`（同种子成对 Vina 比较）。

判据（预登记，§5.41）：① MW 中位 ≥340 且重原子中位 ≥25；② strict 通过率 ≥ 官方基线的 90%；
③ 生成数 ≥ 请求值（并注意 §5.24 的"基线超产"假阳性）。
只有 ①② 通过的候选才进入对接级；③ 单独记录（因为它有已知假阳性，不作为淘汰依据）。

用法：
    python outputs/night_pipeline_v10/dock_v10term_candidates.py            # 自动跑完整个判定
"""
import argparse
import csv
import glob
import io
import os
import subprocess
import sys

# 无人值守运行时 stdout 常被重定向为 GBK 管道：非 GBK 字符（⚠/… 等）会直接崩脚本。
# 这是 A′ v1.0 采样器真实踩过的坑（当时两个臂全崩、CSV 0 行），此处统一加固。
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PY = sys.executable
NIGHT = os.path.join(ROOT, "outputs", "night_pipeline_v10")
SEL = os.path.join(NIGHT, "sel_v10term.csv")
ARM_CSV = os.path.join(NIGHT, "arm_gen_v10term.csv")
BASELINE_ARM = os.path.join("outputs", "ab_v10_s60", "a", "pretrained_Pocket2Mol.pt_s2024")
MW_MIN, HA_MIN, RATIO_MIN = 340.0, 25.0, 0.90


def find_smiles(path):
    if os.path.isfile(path):
        return path
    cands = glob.glob(os.path.join(ROOT, path, "**", "SMILES.txt"), recursive=True)
    return sorted(cands)[-1] if cands else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sel", default=SEL)
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--parallel", type=int, default=6)
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--force", action="store_true", help="即使候选未过 ①② 也做对接（人工覆盖用）")
    args = ap.parse_args()

    if not os.path.exists(ARM_CSV):
        print("缺少 %s —— 先跑 decide_v10term.py" % ARM_CSV)
        return 1
    rows = {r["arm"]: r for r in csv.DictReader(io.open(ARM_CSV, encoding="utf-8-sig"))}
    base = rows.get("official")
    if not base:
        print("arm CSV 里没有 official 行")
        return 1
    b_strict = float(base.get("strict_rate") or 0)

    sel = list(csv.DictReader(io.open(args.sel, encoding="utf-8-sig"))) if os.path.exists(args.sel) else []
    run_dir_by_name = {os.path.splitext(os.path.basename(r["ckpt"]))[0]: r.get("run_dir") for r in sel}

    chosen = []
    for name, r in rows.items():
        if name == "official":
            continue
        mw, ha = float(r.get("mw_med") or 0), float(r.get("ha_med") or 0)
        sr = float(r.get("strict_rate") or 0)
        size_ok = mw >= MW_MIN and ha >= HA_MIN
        strict_ok = bool(b_strict) and (sr / b_strict) >= RATIO_MIN
        if (size_ok and strict_ok) or args.force:
            chosen.append((name, mw, ha, sr))
    if not chosen:
        print("没有候选通过 ①②（尺寸 + strict≥基线90%）。")
        print("→ 按预登记升级阶梯：V1（常数 pos_weight=2.0）不解决尺寸/口径问题，应转 V2"
              "（configs/train_gpcr_mass_v10term_bal.yml）。")
        print("  若仍要对某个具体检查点看对接分，可加 --force。")
        return 0

    print("进入对接级的候选：%s" % ", ".join(c[0] for c in chosen))
    arms = ["official=%s" % find_smiles(BASELINE_ARM)]
    for name, *_ in chosen:
        d = run_dir_by_name.get(name)
        smi = find_smiles(d) if d else None
        if not smi:
            print("跳过 %s：找不到 SMILES.txt（run_dir=%s）" % (name, d))
            continue
        arms.append("%s=%s" % (name, smi))
    if len(arms) < 2:
        print("可用臂不足 2 个（候选分子缺失）")
        return 1

    out_csv = os.path.join(NIGHT, "dock_v10term.csv")
    cmd = [PY, os.path.join(ROOT, "src", "scripts", "dock_smiles_ab.py"),
           "--arms"] + arms + ["--target", args.target,
                               "--parallel", str(args.parallel),
                               "--exhaustiveness", str(args.exhaustiveness)]
    print("对接命令：\n  %s" % " ".join(cmd))
    with io.open(os.path.join(NIGHT, "dock_v10term.log"), "w", encoding="utf-8") as fh:
        r = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    print("对接返回码 %d，日志 %s" % (r.returncode, os.path.join(NIGHT, "dock_v10term.log")))
    # dock_smiles_ab 自己写 CSV，这里只报告位置供下一轮读取
    for cand in glob.glob(os.path.join(ROOT, "outputs", "dock_smiles_ab", "*.csv")):
        print("候选产物：%s" % os.path.relpath(cand, ROOT))
    print("期望产物（若工具写在自己的 work 目录）：%s" % out_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
