# -*- coding: utf-8 -*-
"""A′ 自适应阈值采样：把"**掉进小模式**"的采样批次自动拉回尺寸达标区间。

为什么需要它（实测依据）
------------------------
1. 同一检查点 v10@5500 的生成尺寸随**种子**在"正常/停早"两个模式间跳变
   （n=60/seed2024 → MW 397；n=17/seed2024 → 350.9；n=17/seed2025 → 308.9）。
   全项目 v7–v11 的相邻检查点尺寸 |ΔMW| 中位 59.7、lag-1 自相关 ≈ −0.74；
   按 val loss 选点会**系统性**选到停早检查点（v9/v10/v11/v11b 连续 4 次命中）。
2. 朴素的"多检查点池化"**不能**修复：组分自己就掉进小模式时，池子同样小
   （实测：3×17 池化 → MW 253/263、ΔVina +0.54/+2.11，两臂尺寸比 0.68/0.70）。
3. 但**阈值是主旋钮**：同一失败组合下 thr 0 → MW 308.9，**thr −0.25 → 348.9**
   （重原子 25.5、QED 0.814、SA 2.18），thr −0.5 → 455.5（过头）。
   → 于是有本脚本：**先采一小批探针，按 MW 中位把阈值往反方向调，直到落进达标带**。

策略（A′，预登记）
------------------
- 起点阈值取自**权重自己的推荐值**（`src/paths.py: recommended_frontier_threshold`；
  训练检查点通常无侧车条目 → 0.0；v3 → −0.5），CLI 可显式覆盖；来源会写进 CSV。
- 探针规模 n0=20；达标带 **MW 中位 ∈ [340, 430]**（下限 = §5.2e 预登记尺寸门槛 340；
  上限 = §5.33 处方① 的 430）。
- 不在带内则 **thr ∓ step**（step=0.25；据 §5.33：dMW/dthr ≈ −100…−200，阈值↓则分子↑）；
  最多 `--max-tries` 次调整（默认 2）。
- **无效批次门禁**：子进程 rc≠0 或没产出 SMILES 或 session 早于本次启动 → 标 invalid，
  **不参与调参**（同阈值重试），单列 `invalid_attempts`（1.0 版把失败当"MW 偏大"调错方向，已修）。
- 每完成一个种子**立即**写一行 CSV（增量落盘；1.0 版末尾一次性写，崩溃即全丢）。

v1.1（2026-09-29）修复清单（来自独立校验子代理的审查）
------------------------------------------------------
① stdout 重定向到文件时 Windows 用 cp936，旧版打印 `✅` 触发 UnicodeEncodeError 导致两臂崩溃
   → 现在强制 `sys.stdout.reconfigure(encoding="utf-8")`，且日志只用 ASCII 标记 `[IN]/[ADJ]`；
② CSV 改增量写 + flush；③ 无效批次门禁（见上）；④ 补 ckpt/target/采样参数/脚本版本/时间/阈值来源列；
⑤ try 目录加运行戳 + 校验 session 新鲜度（防静默复用旧批次）；⑥ 末尾另报**跨种子 MW 极差**
   （§5.33 判据 ≤10%）；⑦ top-up 记录 rc 与自身落区状态。

判据（用于"是否把 A′ 作为出货采样策略"）
----------------------------------------
- **落区率 ≥90%**（K≥4 个种子）+ **跨种子 MW 极差 ≤10%**（§5.33 处方①）；
- 通过后再做**对接级双种子 A/B**（`dock_smiles_ab.py`，同 n、同 beam/步数，报 Vina+LE，尺寸差>15% 报警）。

⚠ 显存纪律：本机 8 GB；**训练运行时不要跑本脚本**（12:17 已因并发 OOM 崩过一次训练）。

用法
----
    python src/scripts/adaptive_threshold_sample.py --ckpt logs/<run>/checkpoints/5500.pt \
        --target A2A --seeds 2024 2025 2026 2027 --n0 20 --final-n 50 \
        --work outputs/adaptA --out outputs/adaptA_summary.csv
"""
import argparse
import csv
import glob
import json
import os
import statistics as st
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from paths import ROOT, OUTPUTS, ensure_dir, recommended_frontier_threshold  # noqa: E402

SELECT = os.path.join(ROOT, "src", "scripts", "select_ckpt_by_generation.py")
BAND = (340.0, 430.0)
XSEED_SPREAD_MAX = 10.0        # §5.33 判据：跨种子 MW 极差 ≤10%
SCRIPT_VER = "adaptA-1.1(20260929)"

# stdout 重定向到文件时 Windows 默认 cp936，任何非 GBK 字符都会 UnicodeEncodeError
# （v1.0 就因打印 ✅ 而崩，两臂在 16:01 双双挂掉）→ 强制 UTF-8，且日志只用 ASCII 标记。
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

FIELDS = ["seed", "start_thr", "thr_source", "attempts", "invalid_attempts", "final_thr",
          "in_band", "n", "mw_med", "ha_med", "qed_med", "sa_med", "strict_rate",
          "topup_n", "topup_mw", "topup_ha", "topup_in_band", "topup_rc", "topup_session",
          "final_session", "ckpt", "target", "n0", "beam", "max_steps", "step", "max_tries",
          "script_ver", "run_utc", "traj"]


def log(msg):
    print(msg, flush=True)


def measure(smis):
    """返回 (n, mw_med, ha_med, qed_med, sa_med, strict_rate)。口径与全项目一致。"""
    sys.path.insert(0, os.path.join(ROOT, "src", "scripts"))
    from sweep_frontier_threshold import metrics
    m = metrics(smis)
    if not m:
        return dict(n=0, mw_med=None, ha_med=None, qed_med=None, sa_med=None, strict_rate=None)
    return dict(n=m["n"], mw_med=round(m["mw"], 1), ha_med=round(m["ha"], 1),
                qed_med=round(m["qed"], 4), sa_med=round(m["sa"], 3),
                strict_rate=round(m["pass_pct"] / 100.0, 4))


def sample_once(ckpt, target, n, beam, steps, seed, thr, keep_dirtag):
    """跑一次采样，返回 (smis, session_dir, rc, started_at)。

    keep_dirtag 里带运行戳，避免重跑时**静默复用**上一次的旧 session
    （v1.0 只在 tag 里写 try 序号与阈值，重跑会 glob 到旧目录）。
    """
    ensure_dir(keep_dirtag)
    started = time.time()
    cfg_csv = os.path.join(keep_dirtag, "gen.csv")
    cmd = [sys.executable, SELECT, "--ckpts", ckpt, "--target", target,
           "--num-samples", str(n), "--beam", str(beam), "--max-steps", str(steps),
           "--seed", str(seed), "--rank-by", "both", "--keep-dir", keep_dirtag,
           "--frontier-threshold", str(thr), "--out", cfg_csv]
    with open(os.path.join(keep_dirtag, "run.log"), "a", encoding="utf-8") as lf:
        lf.write("\n$ " + " ".join(cmd) + "\n")
        lf.flush()
        rc = subprocess.run(cmd, cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT).returncode
    sess = sorted(d for d in glob.glob(os.path.join(keep_dirtag, "**", "samples", "*"),
                                       recursive=True) if os.path.isdir(d))
    smis = []
    sessdir = sess[-1] if sess else None
    # session 必须比本次启动新，否则视为旧批次 → 判无效，避免"拿旧结果当本次"
    if sessdir is not None and os.path.getmtime(sessdir) < started - 5:
        sessdir = None
    if sessdir is not None:
        p = os.path.join(sessdir, "SMILES.txt")
        if os.path.exists(p):
            smis = [x.strip() for x in open(p, encoding="utf-8", errors="ignore") if x.strip()]
    return smis, sessdir, rc


def emit(path, row):
    """增量落盘：每完成一个种子就追加一行（崩溃也能保住已完成结果）。"""
    new = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerow(row)
        fh.flush()


def main():
    ap = argparse.ArgumentParser(description="A′ 自适应阈值采样（把尺寸拉回达标带）")
    ap.add_argument("--ckpt", required=True, help="检查点路径（版本标签沿用其训练版本，如 v10@5500）")
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--seeds", nargs="+", type=int, default=[2024, 2025, 2026, 2027])
    ap.add_argument("--n0", type=int, default=20, help="探针批大小")
    ap.add_argument("--final-n", type=int, default=0,
                    help=">0 时，落区后用最终阈值再采一批到 final-n（对接 A/B 用）")
    ap.add_argument("--beam", type=int, default=50)
    ap.add_argument("--max-steps", type=int, default=40)
    ap.add_argument("--step", type=float, default=0.25, help="每次调整的阈值步长")
    ap.add_argument("--max-tries", type=int, default=2, help="最多额外重采次数")
    ap.add_argument("--start-thr", type=float, default=None,
                    help="起点阈值；不传=用权重推荐值（paths.recommended_frontier_threshold）")
    ap.add_argument("--work", default=os.path.join(OUTPUTS, "adaptA"))
    ap.add_argument("--out", default=os.path.join(OUTPUTS, "adaptA_summary.csv"))
    args = ap.parse_args()

    # peek_ckpt=False：只读侧车表（按文件名），不 torch.load —— 训练检查点通常查不到，
    # 因此来源会如实写成"默认(原版行为)"（数值=0.0，与 v10 的评测口径一致）。
    rec_thr, rec_src = recommended_frontier_threshold(args.ckpt, peek_ckpt=False)
    start_thr = rec_thr if args.start_thr is None else float(args.start_thr)
    run_utc = time.strftime("%Y-%m-%dT%H:%M:%S")
    run_stamp = time.strftime("%Y%m%d_%H%M%S")
    log("=" * 84)
    log("A' adaptive threshold sampling (v1.1) | ckpt=%s | target=%s | seeds=%s" % (
        os.path.basename(args.ckpt), args.target, args.seeds))
    log("  start_thr %+.2f (source: %s) | n0=%d | band MW[%.0f, %.0f] | step %.2f | max adjusts %d"
        % (start_thr, rec_src, args.n0, BAND[0], BAND[1], args.step, args.max_tries))
    log("  out=%s | run_utc=%s" % (os.path.relpath(args.out, ROOT), run_utc))
    log("=" * 84)

    rows = []
    for seed in args.seeds:
        thr = start_thr
        traj = []
        final = None
        invalid = 0
        for attempt in range(args.max_tries + 1):
            tag = os.path.join(args.work, "s%d" % seed,
                               "run%s_try%d_thr%+.2f" % (run_stamp, attempt, thr))
            smis, sess, rc = sample_once(args.ckpt, args.target, args.n0, args.beam,
                                         args.max_steps, seed, thr, tag)
            m = measure(smis)
            valid = (rc == 0 and m["n"] > 0 and sess is not None)
            in_band = bool(valid and m["mw_med"] is not None and BAND[0] <= m["mw_med"] <= BAND[1])
            if not valid:
                invalid += 1
            traj.append(dict(attempt=attempt, thr=thr, rc=rc, valid=bool(valid),
                             session=sess, in_band=bool(in_band), **m))
            log("  seed %d | try %d | thr %+.2f | n=%s | MW %s | HA %s | QED %s | SA %s | strict %s | %s"
                % (seed, attempt, thr, m["n"], m["mw_med"], m["ha_med"], m["qed_med"],
                   m["sa_med"], m["strict_rate"],
                   "IN-BAND" if in_band else ("INVALID" if not valid else "ADJUST")))
            final = traj[-1]
            if in_band or attempt == args.max_tries:
                break
            if not valid:
                log("    [WARN] invalid batch (rc=%s n=%s session=%s) -> threshold NOT changed, retry same"
                    % (rc, m["n"], sess))
                continue
            thr = thr - args.step if m["mw_med"] < BAND[0] else thr + args.step

        topup = None
        if args.final_n and final and final["in_band"] and final["n"] < args.final_n:
            tag = os.path.join(args.work, "s%d" % seed, "run%s_topup_thr%+.2f" % (run_stamp, final["thr"]))
            smis, sess, rc = sample_once(args.ckpt, args.target, args.final_n, args.beam,
                                         args.max_steps, seed, final["thr"], tag)
            mm = measure(smis)
            topup = dict(**mm, session=sess, thr=final["thr"], rc=rc,
                         in_band=bool(rc == 0 and mm["n"] > 0 and mm["mw_med"] is not None
                                      and BAND[0] <= mm["mw_med"] <= BAND[1]))
            log("  seed %d | top-up thr %+.2f rc=%s n=%d | MW %s | HA %s | strict %s | %s"
                % (seed, final["thr"], rc, topup["n"], topup["mw_med"], topup["ha_med"],
                   topup["strict_rate"], "IN-BAND" if topup["in_band"] else "OUT"))
        row = dict(seed=seed, start_thr=start_thr, thr_source=rec_src,
                   attempts=len(traj), invalid_attempts=invalid,
                   final_thr=final["thr"], in_band=final["in_band"],
                   n=final["n"], mw_med=final["mw_med"], ha_med=final["ha_med"],
                   qed_med=final["qed_med"], sa_med=final["sa_med"],
                   strict_rate=final["strict_rate"],
                   topup_n=(topup or {}).get("n"), topup_mw=(topup or {}).get("mw_med"),
                   topup_ha=(topup or {}).get("ha_med"),
                   topup_in_band=(topup or {}).get("in_band"), topup_rc=(topup or {}).get("rc"),
                   topup_session=(topup or {}).get("session"),
                   final_session=final["session"],
                   ckpt=os.path.relpath(args.ckpt, ROOT), target=args.target, n0=args.n0,
                   beam=args.beam, max_steps=args.max_steps, step=args.step,
                   max_tries=args.max_tries, script_ver=SCRIPT_VER, run_utc=run_utc,
                   traj=json.dumps(traj, ensure_ascii=False))
        rows.append(row)
        emit(args.out, row)          # 增量落盘

    valid_rows = [r for r in rows if r["mw_med"] is not None]
    rate = (sum(1 for r in valid_rows if r["in_band"]) / len(valid_rows)) if valid_rows else 0.0
    mws = [r["mw_med"] for r in valid_rows]
    spread = (max(mws) - min(mws)) / st.median(mws) * 100.0 if mws else None
    base_in = sum(1 for r in valid_rows
                  if r["traj"] and json.loads(r["traj"])[0].get("in_band"))
    log("\n" + "=" * 84)
    log("in-band rate = %d/%d = %.1f%% (criterion >=90%%) -> %s"
        % (sum(1 for r in valid_rows if r["in_band"]), len(valid_rows), 100 * rate,
           "PASS" if rate >= 0.90 else "FAIL"))
    log("baseline (try0, no adaptation) in-band = %d/%d = %.1f%%"
        % (base_in, len(valid_rows), 100.0 * base_in / len(valid_rows) if valid_rows else 0.0))
    if spread is not None:
        log("cross-seed MW median %.1f | spread %.1f (%.1f%%) (criterion <=%.0f%%) -> %s"
            % (st.median(mws), max(mws) - min(mws), spread, XSEED_SPREAD_MAX,
               "PASS" if spread <= XSEED_SPREAD_MAX else "FAIL"))
    if any(r["invalid_attempts"] for r in rows):
        log("note: %d seed(s) had invalid batches (rc!=0 / no product); retried same threshold, "
            "not used for tuning" % sum(1 for r in rows if r["invalid_attempts"]))
    for r in rows:
        log("  seed %d: final thr %+.2f (adjusts %d, invalid %d) | MW %s | HA %s | strict %s"
            % (r["seed"], r["final_thr"], r["attempts"] - 1, r["invalid_attempts"],
               r["mw_med"], r["ha_med"], r["strict_rate"]))
    log("summary: %s (incremental)" % os.path.relpath(args.out, ROOT))
    return 0 if rate >= 0.90 else 1


if __name__ == "__main__":
    sys.exit(main())
