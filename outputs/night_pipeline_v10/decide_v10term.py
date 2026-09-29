# -*- coding: utf-8 -*-
"""v10term 生成级扫描结束后的**决策**步骤：套预登记判据选点，并打印对接 A/B 命令。

预登记判据（§5.41 已写死，不得临时发明）：
  ① 尺寸：MW 中位 ≥340 且重原子中位 ≥25（`arm_gen_metrics.py` 的 ok_size）
  ② strict 通过率 ≥ 官方基线的 90%
  ③ 生成数：≥ 请求值（n=60）——注意 §5.24 的"基线超产"假阳性
  ④ 通过 ①②③ 的候选里，按 strict 通过率降序取前 2 个做对接级 A/B（同种子成对 vs 官方）

基线 arm 复用既有同协议采样（官方权重，n=60/beam50/steps40/seed2024，阈值 0.0），
**不重跑采样**；扫描产物由 run_v10term_eval.py 留在
outputs/night_pipeline_v10/v10term_sel/<iter> 下。

用法：
    python outputs/night_pipeline_v10/decide_v10term.py [--baseline <arm 目录>]
"""
import argparse
import csv
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
SEL_CSV = os.path.join(NIGHT, "sel_v10term.csv")
KEEP = os.path.join(NIGHT, "v10term_sel")
MW_MIN, HA_MIN, STRICT_RATIO_MIN = 340.0, 25.0, 0.90


def log(msg, fh=None):
    print(msg, flush=True)
    if fh:
        fh.write(msg + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline",
                    default=os.path.join("outputs", "ab_v10_s60", "a", "pretrained_Pocket2Mol.pt_s2024"),
                    help="官方基线的 arm 目录。⚠ 默认这个是 **n=50/steps60** 的旧协议臂；"
                         "扫描用的是 n=60/steps40 —— 报告里会显式标注该口径差，"
                         "同协议官方臂由 outputs/night_pipeline_v10/run_s1_experiment.py 产出。")
    ap.add_argument("--baseline-protocol", default="n=50/steps60（旧臂；与本次扫描 n=60/steps40 不完全同协议）",
                    help="写进报告的口径说明，避免后来者误读")
    ap.add_argument("--n-requested", type=int, default=60)
    ap.add_argument("--topk", type=int, default=2)
    ap.add_argument("--sel", default=SEL_CSV, help="扫描产物 CSV（默认 sel_v10term.csv）")
    ap.add_argument("--skip-arm-gen", action="store_true",
                    help="复用已有的 arm_gen_v10term.csv（不重算尺寸/strict）")
    args = ap.parse_args()

    if not os.path.exists(args.sel):
        print("找不到 %s —— 扫描还没跑完？" % args.sel)
        return 1
    rows = list(csv.DictReader(io.open(args.sel, encoding="utf-8-sig")))
    if not rows:
        print("扫描 CSV 为空")
        return 1
    cands = [(r["ckpt"], r.get("run_dir") or os.path.join(
        KEEP, os.path.splitext(os.path.basename(r["ckpt"]))[0])) for r in rows]

    with io.open(os.path.join(NIGHT, "decide_v10term_report.md"), "w", encoding="utf-8") as fh:
        log("# v10term 生成级判定（自动生成，判据 §5.41 预登记）\n", fh)
        log("> ⚠ **口径提示**：官方基线臂 = `%s`\n> 口径 = %s。\n"
            "> 若要与候选做**严格同协议**比较，请用 `run_s1_experiment.py` 产出的同协议官方臂"
            "（n=60/beam50/steps40，种子 2024/2025）替换 `--baseline`。\n"
            "> 尺寸/strict 的比较对协议敏感（steps 决定尺寸上限），**跨协议只能当参考**。\n"
            % (args.baseline, args.baseline_protocol), fh)
        arms = ["official=%s" % args.baseline]
        for name, d in cands:
            tag = os.path.splitext(os.path.basename(name))[0]
            arms.append("%s=%s" % (tag, d))
        out_csv = os.path.join(NIGHT, "arm_gen_v10term.csv")
        cmd = [PY, os.path.join(ROOT, "src", "scripts", "arm_gen_metrics.py"),
               "--arms"] + arms + ["--baseline", "official", "--out", out_csv]
        log("```\n%s\n```\n" % " ".join(cmd), fh)
        if args.skip_arm_gen and os.path.exists(out_csv):
            log("（--skip-arm-gen：复用已有 arm CSV，未重算）\n", fh)
        else:
            r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, errors="ignore")
            log("```\n%s\n```\n" % (r.stdout or "")[-4000:], fh)
        if not os.path.exists(out_csv):
            log("arm_gen_metrics 未产出 CSV，终止（见上面的 stdout）", fh)
            return 1

        m = {x["arm"]: x for x in csv.DictReader(io.open(out_csv, encoding="utf-8-sig"))}
        base = m.get("official")
        if base is None:
            log("CSV 里没有 official 行", fh)
            return 1
        b_strict = float(base.get("strict_rate") or 0)

        passed, failed = [], []
        for name, _d in cands:
            tag = os.path.splitext(os.path.basename(name))[0]
            x = m.get(tag)
            if x is None:
                failed.append((tag, "未评测"))
                continue
            mw = float(x.get("mw_med") or 0)
            ha = float(x.get("ha_med") or 0)
            sr = float(x.get("strict_rate") or 0)
            n = int(float(x.get("n_finished") or x.get("n") or 0))   # CSV 列名是 n_finished
            size_ok = mw >= MW_MIN and ha >= HA_MIN
            strict_ok = bool(b_strict) and (sr / b_strict) >= STRICT_RATIO_MIN
            count_ok = n >= args.n_requested
            note = "尺寸%s strict%s(%.1f%% vs 基线 %.1f%%, 比值 %.2f) 生成数%d%s" % (
                "OK" if size_ok else "**不足**", "OK" if strict_ok else "**不足**",
                100 * sr, 100 * b_strict, (sr / b_strict if b_strict else 0), n,
                "" if count_ok else " **低于请求值**")
            (passed if (size_ok and strict_ok) else failed).append((tag, note))
            log("- %-10s MW %6.1f HA %5.1f | %s" % (tag, mw, ha, note), fh)

        log("\n## 结论\n", fh)
        if not passed:
            log("**无候选同时满足 ①尺寸 ②strict** → 本档（V1，常数 pos_weight=2.0）"
                "不能解决尺寸/口径问题 → 按升级阶梯进入 **V2**"
                "（`configs/train_gpcr_mass_v10term_bal.yml`，逐样本平衡）。", fh)
            return 0
        log("通过预登记判据（单种子 seed 2024）的候选：%s" % ", ".join(t for t, _ in passed), fh)
        log("\n⚠ **单种子只能过滤、不能选点**（§5.33：单种子既不能给检查点排序，也不能代表其水平）。"
            "故下一步先做**多种子稳定性**（§5.33 处方①：K≥4，判据 `340≤MW≤430` 且跨种子极差 ≤10%），"
            "**再做对接级 A/B**。\n", fh)
        for tag, _ in passed[:args.topk]:
            log("### 候选 `%s`\n" % tag, fh)
            log("**第 1 步（GPU 采样，K=4 种子）**：多种子尺寸稳定性 + 生成级四项\n"
                "```powershell\n"
                "& '%s' src\\scripts\\eval_checkpoint_targets.py --ckpt <候选权重.pt> "
                "--targets A2A --seeds 2024 2025 2026 2027 --num-samples 60 --beam 50 --max-steps 40 "
                "--parallel 6\n```\n" % PY, fh)
            log("判据：**MW 中位落在 340–430** 且 **跨种子极差 ≤10%**；命中者才算「可晋级候选」。\n", fh)
            log("**第 2 步（CPU 对接，复用已采分子，不再占 GPU）**：\n"
                "```powershell\n"
                "& '%s' outputs\\night_pipeline_v10\\dock_v10term_candidates.py --force\n```\n" % PY, fh)
            log("同种子成对 vs 官方；**双口径**（原始 Vina + LE）+ 尺寸 15% 可比性检查。\n", fh)
        log("\n⚠ 另按 §5.41 补一条：对最终候选**顺带扫阈值**（改损失权重会改终止校准）；"
            "并注意 §5.24 的生成数硬规则假阳性。", fh)
    print("报告: %s" % os.path.join(NIGHT, "decide_v10term_report.md"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
