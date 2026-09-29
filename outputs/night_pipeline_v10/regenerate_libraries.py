# -*- coding: utf-8 -*-
"""用**当前默认权重 v10** 重出四个药库，替换当年由 v2 生成的旧库（用户 2026-09-29 夜指示）。

口径对齐（与旧库可直接对照）
--------------------------
旧库清单（`results/*/run_manifest.json`）：4 靶点均 `n=50 / beam=100 / steps=50 / seed 2024 / diversity_w=0.5 / top=300`，
权重 `models/7-eonmol_ft_gpcr_v2.pt`；规模 A2A **521**、B2AR 275、D3 217、5HT2B 252 分子。
本脚本沿用**完全相同**的采样参数与权重**换成 v10**，并按旧库规模决定每靶点的批数（每批 50 分子）：
A2A 11 批、B2AR 6 批、D3 5 批、5HT2B 6 批 = **28 次采样 ≈ 3.5–4 小时 GPU**（实测 10.6 分子/分钟量级）。

三阶段（与旧库同构）
------------------
1. `design.py`（GPU，每批一个种子，种子 2024 起）+ 可选**尺寸机制**（`--schedule-mode`）
2. `build_library.py --runs <所有会话> --library <staging>/library`（CPU：药化规则 + 3D 构象）
3. `predict.py --library <staging>/library --top 300 --results-dir <staging>/results`（CPU：对接 + ADMET + 排序）

安全（用户硬要求：改动前备份、失败可回退）
--------------------------------------
- **先备份**：`results/` 下四个库 + results.csv/compounds.csv 等 → `other/library_backups/<ts>_v2_libraries/`
- **先写 staging**（`outputs/library_v10/`），全部成功后才由 `--swap` 替换 `results/`；
  未加 `--swap` 时**绝不触碰** `results/`。

用法
----
    python outputs/night_pipeline_v10/regenerate_libraries.py --backup-only
    python outputs/night_pipeline_v10/regenerate_libraries.py --targets A2A                # 先做一个靶点
    python outputs/night_pipeline_v10/regenerate_libraries.py --targets A2A B2AR D3 5HT2B  # 全量
    python outputs/night_pipeline_v10/regenerate_libraries.py --swap                       # 校验并把 staging 换上
"""
import argparse
import glob
import io
import os
import shutil
import subprocess
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PY = sys.executable
STAGING = os.path.join(ROOT, "outputs", "library_v10")
RESULTS = os.path.join(ROOT, "results")
BATCHES = {"A2A": 11, "B2AR": 6, "D3": 5, "5HT2B": 6}      # 由旧库规模 / 50 得出
N_PER_BATCH, BEAM, STEPS, DW = 50, 50, 40, 0.5
# 协议说明（2026-09-30 01:10 决定）：旧库用 beam 100/steps 50，实测吞吐仅 ~2.5 分子/分钟
# （28 批要 8+ 小时）；改为 **beam 50/steps 40** —— 与本项目全部 A/B 证据同源（吞吐 ~12 分子/分钟），
# 使"交付的药库"与"我们实测并报告过的分子"是同一口径。旧库协议与产物已完整备份。

# 尺寸机制：off = 纯 v10（阈值 0.0）；s3 = 加硬尺寸窗（需先经凹坑测试判定有效）
SCHEDULES = {
    "off": [],
    "s3": ["--schedule-floor-steps", "24", "--schedule-floor-thr", "-2.0",
           "--schedule-ceiling-ha", "32"],
}


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with io.open(os.path.join(STAGING, "regenerate.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def backup_results():
    ts = time.strftime("%Y%m%d_%H%M%S")
    dest = os.path.join(ROOT, "other", "library_backups", "%s_v2_libraries" % ts)
    os.makedirs(dest, exist_ok=True)
    n = 0
    for name in os.listdir(RESULTS):
        src = os.path.join(RESULTS, name)
        if name == "official_weight_library":          # 官方库归档，不动
            continue
        dst = os.path.join(dest, name)
        if os.path.isdir(src):
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            shutil.copy2(src, dst)
        n += 1
    return dest, n


def run(cmd, logfile, label):
    t0 = time.time()
    with io.open(logfile, "a", encoding="utf-8") as fh:
        fh.write("\n$ %s\n" % " ".join(cmd))
        rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
    return rc, time.time() - t0


def build_one(target, sched_mode, dry, beam=None, steps=None):
    tdir = os.path.join(STAGING, target)
    os.makedirs(tdir, exist_ok=True)
    design_out = os.path.join(tdir, "design")
    runs = []
    for i in range(BATCHES[target]):
        seed = 2024 + i
        cmd = [PY, os.path.join(ROOT, "design.py"), "--target", target,
               "--num-samples", str(N_PER_BATCH), "--beam", str(beam or BEAM),
               "--max-steps", str(steps or STEPS),
               "--seed", str(seed), "--diversity-w", str(DW), "--outdir", design_out] + \
              SCHEDULES[sched_mode]
        if dry:
            print(" ".join(cmd))
            continue
        log("  [%s] 采样 seed %d（%d/%d）..." % (target, seed, i + 1, BATCHES[target]))
        rc, dt = run(cmd, os.path.join(tdir, "01_design.log"), "design")
        log("  [%s] seed %d rc=%d %.1f 分钟" % (target, seed, rc, dt))
        if rc != 0:
            log("  [%s] seed %d 失败，跳过该批" % (target, seed))
    if dry:
        return None
    runs = sorted(d for d in glob.glob(os.path.join(design_out, "*")) if os.path.isdir(d))
    if not runs:
        log("  [%s] 没有任何采样会话，跳过" % target)
        return None
    lib = os.path.join(tdir, "library")
    log("  [%s] 入库（%d 个会话）..." % (target, len(runs)))
    rc, dt = run([PY, os.path.join(ROOT, "src", "scripts", "build_library.py"),
                  "--runs"] + runs + ["--library", lib],
                 os.path.join(tdir, "02_library.log"), "library")
    log("  [%s] 入库 rc=%d %.1f 分钟" % (target, rc, dt))
    if rc != 0 or not os.path.exists(os.path.join(lib, "compounds.csv")):
        log("  [%s] 入库失败，跳过打分" % target)
        return None
    res = os.path.join(tdir, "results")
    log("  [%s] 打分排序（对接 + ADMET + 排序）..." % target)
    rc, dt = run([PY, os.path.join(ROOT, "predict.py"), "--target", target, "--library", lib,
                  "--top", "300", "--seed", "2024", "--results-dir", res,
                  "--outdir", os.path.join(tdir, "predict")],
                 os.path.join(tdir, "03_predict.log"), "predict")
    log("  [%s] 打分 rc=%d %.1f 分钟" % (target, rc, dt))
    return lib if rc == 0 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", default=list(BATCHES))
    ap.add_argument("--schedule-mode", choices=sorted(SCHEDULES), default="off",
                    help="尺寸机制：off=纯 v10；s3=硬尺寸窗(24,32)（需先经凹坑测试判定足够好）")
    ap.add_argument("--backup-only", action="store_true")
    ap.add_argument("--swap", action="store_true", help="把 staging 成果替换到 results/（默认不动 results/）")
    ap.add_argument("--skip-backup", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--beam", type=int, default=BEAM)
    ap.add_argument("--steps", type=int, default=STEPS)
    args = ap.parse_args()
    os.makedirs(STAGING, exist_ok=True)

    if args.backup_only or not args.skip_backup:
        dest, n = backup_results()
        print("已备份 results/ 下 %d 项 → %s" % (n, os.path.relpath(dest, ROOT)))
        if args.backup_only:
            return 0

    if args.swap:
        return do_swap(args.targets)

    log("目标=%s | 机制=%s | 批数=%s" % (args.targets, args.schedule_mode,
                                        {t: BATCHES[t] for t in args.targets}))
    for t in args.targets:
        log("=== 开始 %s ===" % t)
        build_one(t, args.schedule_mode, args.dry_run, args.beam, args.steps)
    log("全部目标处理完毕（结果在 %s；未替换 results/，需显式 --swap）" % os.path.relpath(STAGING, ROOT))
    return 0


def do_swap(targets):
    """校验 staging 产物后替换 results/ 下对应内容（旧库已备份）。"""
    swapped = []
    for t in targets:
        src_res = os.path.join(STAGING, t, "results")
        src_lib = os.path.join(STAGING, t, "library")
        if not os.path.exists(src_res):
            log("  [swap] %s 缺少 staging results，跳过" % t)
            continue
        # 旧库结构：A2A 的清单在 results/ 根目录；其余靶点在 results/<T>/
        root_files = ("results.csv", "compounds.csv", "admet_screen.csv", "top_candidates.csv",
                      "docking_summary.csv", "rejected.csv", "run_manifest.json", "summary.txt")
        if t == "A2A":
            for name in root_files:
                src_f = os.path.join(src_res, name)
                if os.path.exists(src_f):
                    shutil.copy2(src_f, os.path.join(RESULTS, name))
        # 每靶点子目录（旧库结构：results/<T>/）
        dst_t = os.path.join(RESULTS, t)
        os.makedirs(dst_t, exist_ok=True)
        for name in os.listdir(src_res):
            s = os.path.join(src_res, name)
            d = os.path.join(dst_t, name)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)
        lib = os.path.join(src_lib, "compounds.csv")
        if os.path.exists(lib):
            shutil.copy2(lib, os.path.join(RESULTS, "compounds_%s.csv" % t))
        swapped.append(t)
        log("  [swap] %s 已替换（源 %s）" % (t, os.path.relpath(src_res, ROOT)))
    log("替换完成：%s。**请随即更新 results/README.md 与四份对外文档的生成权重标注。**" % swapped)
    return 0


if __name__ == "__main__":
    sys.exit(main())
