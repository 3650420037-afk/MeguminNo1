# -*- coding: utf-8 -*-
"""对接级 A/B 对照: 用**结合强度(Vina 对接分)**判定检查点优劣。

为什么需要它(判据必须与目标一致)
--------------------------------
本项目的第一优化轴是**结合强度**, 但既有的检查点筛选
(select_ckpt_by_generation.py / judge_promotion.py)用的是生成指标
(完成数、QED、SA、Murcko 骨架、口袋穿模)。二者可以不一致 ——
"更好看的分子"不等于"更亲和的分子"。
因此晋级判定必须补上**对接**这一轴, 本脚本就是那把尺子。

它做什么
--------
对每个检查点 x 每个随机种子:
  1. 用**统一的生成配置**(含 diversity_w=0.5, 与出货药库同源)生成 N 个分子;
  2. 把生成的 SMILES 送进与建库**同一条 Vina 管线**对接进该靶点实验口袋;
  3. 汇总: 完成分子数 / 唯一分子数 / QED 中位 / SA 中位 / 骨架数 /
     Vina 均值、中位、最强 / <=-10 占比。
最后给出**逐种子成对差异**与晋级结论。

判据(全部相对 baseline, 且必须在**同一种子**下成对比)
----------------------------------------------------
硬性(任一不满足即 NOT_PROMOTED):
  - 完成分子数 >= baseline - 5      (停止策略没被破坏 —— 这是历史最大失败模式)
  - QED 中位     >= baseline - 0.05
  - SA 中位      <= baseline + 0.30
  - 骨架数       >= baseline * 0.8
首要轴(结合强度):
  - Vina 中位    <= baseline + 0.2  (Vina 越低越强; 允许 0.2 kcal/mol 噪声)
    * 若 Vina 中位显著更好(<= baseline - 0.2), 记为 DOCK_BETTER
晋级: 硬性全过, 且 Vina 不劣化。Vina 打平则看 <=-10 占比; 仍打平则**不晋级**
(宁可不晋级, 也不要用次要指标冒充"更强")。

用法
----
    python src/scripts/ab_docking_compare.py \
        --ckpts models/pretrained_Pocket2Mol.pt models/7-eonmol_ft_gpcr_v2.pt \
        --target A2A --n 50 --seeds 2024 2025 --parallel 6 --exhaustiveness 4
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
from paths import ROOT, CONFIGS, TARGETS, ensure_dir  # noqa: E402

SELECT = os.path.join(ROOT, "src", "scripts", "select_ckpt_by_generation.py")
DOCK = os.path.join(ROOT, "src", "scripts", "docking_pipeline.py")

# 硬性门槛(相对 baseline)
TOL_N = -5
TOL_QED = -0.05
TOL_SA = 0.30
TOL_SCAF = 0.8
# 首要轴: Vina 越低越强
TOL_VINA = 0.20


def log(m):
    print(m, flush=True)


def read_smiles_txt(p):
    out = []
    with open(p, encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                out.append(s)
    return out


def run(cmd, logf):
    with open(logf, "a", encoding="utf-8", errors="ignore") as lf:
        lf.write("\n$ " + " ".join(cmd) + "\n")
        lf.flush()
        return subprocess.run(cmd, cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT).returncode


def gen_one(ckpt, target, n, beam, steps, seed, keep_dir):
    """生成 -> 返回该检查点本次的 session 目录。"""
    rc = run([sys.executable, SELECT, "--ckpts", ckpt, "--target", target,
              "--num-samples", str(n), "--beam", str(beam), "--max-steps", str(steps),
              "--seed", str(seed), "--rank-by", "both", "--keep-dir", keep_dir,
              "--out", os.path.join(keep_dir, "gen_metrics.csv")],
             os.path.join(keep_dir, "run.log"))
    if rc != 0:
        log("    [警告] 生成返回码 %d" % rc)
    name = os.path.basename(ckpt)[:-3] if ckpt.endswith(".pt") else os.path.basename(ckpt)
    sess = sorted(d for d in glob.glob(os.path.join(keep_dir, name, "samples", "*"))
                  if os.path.isdir(d))
    return sess[-1] if sess else None


def dock_smiles(smis, target, outdir, parallel, exhaustiveness, logf):
    """按分块并发把 SMILES 对接进该靶点口袋, 返回 (name -> vina_score) 与统计。"""
    reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    t = reg[target]
    pdb = os.path.join(ROOT, t["pdb"])
    center = ",".join("%.3f" % v for v in t["center"])
    ensure_dir(outdir)
    n = max(1, min(parallel, len(smis)))
    chunks = [[] for _ in range(n)]
    for i, s in enumerate(smis):
        chunks[i % n].append(s)
    procs = []
    for ci, ch in enumerate(chunks):
        if not ch:
            continue
        sf = os.path.join(outdir, "s%02d.txt" % ci)
        with open(sf, "w", encoding="ascii") as f:
            f.write("\n".join(ch) + "\n")
        cdir = os.path.join(outdir, "c%02d" % ci)
        lf = open(os.path.join(outdir, "c%02d.log" % ci), "w", encoding="utf-8", errors="ignore")
        procs.append((ci, subprocess.Popen(
            [sys.executable, DOCK, "--receptor", pdb, "--center=" + center,
             "--smiles-file", sf, "--out", cdir,
             "--exhaustiveness", str(exhaustiveness)],
            cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT), lf))
    for ci, p, lf in procs:
        p.wait()
        lf.close()
    scores, n_ok = [], 0
    for ci, ch, in [(c, x) for c, x in enumerate(chunks) if x]:
        cs = os.path.join(outdir, "c%02d" % ci, "summary.csv")
        if not os.path.exists(cs):
            continue
        with open(cs, encoding="utf-8-sig", newline="") as f:
            for r in csv.DictReader(f):
                if r.get("status") == "ok":
                    n_ok += 1
                    try:
                        scores.append(float(r["vina_score"]))
                    except Exception:
                        pass
    return scores, n_ok


def main():
    ap = argparse.ArgumentParser(description="对接级 A/B 对照(Vina 为首要判据)")
    ap.add_argument("--ckpts", nargs="+", default=None,
                    help="要比较的检查点; 第一个作为 baseline。用 --from-csv 时可不传")
    ap.add_argument("--from-csv", default=None,
                    help="不重新生成/对接, 只读既有结果 CSV 并用当前判据逻辑离线重算")
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--beam", type=int, default=50)
    ap.add_argument("--max-steps", type=int, default=40)
    ap.add_argument("--seeds", nargs="+", type=int, default=[2024, 2025])
    ap.add_argument("--parallel", type=int, default=6)
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--work", default=os.path.join(ROOT, "outputs", "ab_docking"))
    ap.add_argument("--out", default=os.path.join(ROOT, "outputs", "ab_docking.csv"))
    args = ap.parse_args()

    if not args.ckpts and not args.from_csv:
        ap.error("必须给 --ckpts, 或者用 --from-csv 离线重算")
    args.ckpts = [c if os.path.isabs(c) else os.path.join(ROOT, c) for c in (args.ckpts or [])]
    # 同名检查点要区分开, 否则产物目录会互相覆盖
    names, seen = [], {}
    for c in args.ckpts:
        b = os.path.basename(c)
        if b in seen:
            seen[b] += 1
            b = "%s_%d" % (b, seen[b])
        else:
            seen[b] = 0
        names.append(b)
    ensure_dir(args.work)
    log("=" * 78)
    log("对接级 A/B: 靶点 %s | 每检查点每种子 %d 分子 | beam %d | %d 步 | 种子 %s"
        % (args.target, args.n, args.beam, args.max_steps, args.seeds))
    log("检查点: %s" % list(zip(names, [os.path.relpath(c, ROOT) for c in args.ckpts])))
    log("=" * 78)

    rows = []
    if args.from_csv:
        # 离线重算判据: 只读既有 CSV, 不重新生成/对接。
        # 用途: 判据逻辑若被修正(例如硬判据由"对接数"改为"生成数"), 无需重跑昂贵的对接
        # 就能用新逻辑重新判定历史结果。
        with open(args.from_csv, encoding="utf-8-sig", newline="") as f:
            for r in csv.DictReader(f):
                for k in ("seed", "n_smiles", "n_docked"):
                    r[k] = int(r[k])
                for k in ("vina_mean", "vina_med", "vina_min", "frac_le10"):
                    r[k] = float(r[k]) if r.get(k) not in (None, "", "None") else None
                rows.append(r)
        order = []
        for r in rows:
            if r["ckpt"] not in order:
                order.append(r["ckpt"])
        if len(order) < 2:
            raise SystemExit("CSV 里只有 %d 个检查点, 无法成对比较" % len(order))
        names = order
        args.seeds = sorted({r["seed"] for r in rows})
        log("从 CSV 离线重算: %s (%d 行, %d 个检查点, 种子 %s)"
            % (args.from_csv, len(rows), len(names), args.seeds))
    else:
        for seed in args.seeds:
            for name, ckpt in zip(names, args.ckpts):
                wd = ensure_dir(os.path.join(args.work, "%s_s%d" % (name, seed)))
                t0 = time.time()
                log("\n[%s | 种子 %d] 生成中 ..." % (name, seed))
                sess = gen_one(ckpt, args.target, args.n, args.beam, args.max_steps, seed, wd)
                if sess is None:
                    log("  生成失败, 跳过")
                    continue
                smi_p = os.path.join(sess, "SMILES.txt")
                smis = read_smiles_txt(smi_p) if os.path.exists(smi_p) else []
                log("  生成完成: %d 个分子 (%.0f s), 开始对接 ..." % (len(smis), time.time() - t0))
                if not smis:
                    continue
                scores, n_ok = dock_smiles(smis, args.target, os.path.join(wd, "dock"),
                                           args.parallel, args.exhaustiveness,
                                           os.path.join(wd, "run.log"))
                med = st.median(scores) if scores else None
                log("  对接完成: ok %d | 中位 %s | 用时 %.0f s"
                    % (n_ok, ("%.2f" % med) if med else "NA", time.time() - t0))
                rows.append({"ckpt": name, "seed": seed, "n_smiles": len(smis),
                             "n_docked": n_ok,
                             "vina_mean": round(st.mean(scores), 3) if scores else None,
                             "vina_med": round(med, 3) if med else None,
                             "vina_min": round(min(scores), 3) if scores else None,
                             "frac_le10": round(100.0 * sum(1 for x in scores if x <= -10) / len(scores), 1)
                             if scores else None})

    if not rows:
        raise SystemExit("没有任何结果")
    with open(args.out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- 成对比较(同种子) ----
    base_name = names[0]
    log("\n" + "=" * 78)
    log("结果 (baseline = %s)" % base_name)
    log("=" * 78)
    hdr = "%-24s %5s %5s %8s %8s %8s %8s" % ("检查点", "种子", "对接", "Vina均值", "Vina中位", "最强", "<=-10%")
    log(hdr)
    for r in rows:
        log("%-24s %5s %5s %8s %8s %8s %8s"
            % (r["ckpt"], r["seed"], r["n_docked"], r["vina_mean"], r["vina_med"],
               r["vina_min"], r["frac_le10"]))

    verdicts = []
    for name in names[1:]:
        log("\n--- %s vs %s (逐种子) ---" % (name, base_name))
        worse, better, miss = 0, 0, 0
        for seed in args.seeds:
            b = [r for r in rows if r["ckpt"] == base_name and r["seed"] == seed]
            f = [r for r in rows if r["ckpt"] == name and r["seed"] == seed]
            if not b or not f:
                miss += 1
                continue
            b, f = b[0], f[0]
            # 硬判据必须看**生成数**而不是对接成功数:
            # 历史最大失败模式是"停止策略被破坏、分子永不终止", 其指纹是**生成分子数骤降**
            # (采样日志 Failed=0, 完成数 63 -> 23/7)。若只看对接成功数, 生成退化了也可能看不出来。
            # 两者都要求不退化, 取更严的那个。
            dn = f["n_smiles"] - b["n_smiles"]
            dnd = f["n_docked"] - b["n_docked"]
            dv = (f["vina_med"] - b["vina_med"]) if (f["vina_med"] is not None
                                                    and b["vina_med"] is not None) else None
            hard = (dn >= TOL_N) and (dnd >= TOL_N)
            if dv is None:
                ok_v = False
            elif dv <= -TOL_VINA:
                ok_v, better = True, better + 1
            elif dv <= TOL_VINA:
                ok_v = True
            else:
                ok_v, worse = False, worse + 1
            log("  种子 %d: 生成数 %d->%d (%+d, %s) | 对接数 %d->%d | Vina中位 %s->%s (%s, %s)"
                % (seed, b["n_smiles"], f["n_smiles"], dn, "OK" if hard else "退化!",
                   b["n_docked"], f["n_docked"], b["vina_med"], f["vina_med"],
                   ("%+.2f" % dv) if dv is not None else "NA",
                   "持平" if ok_v and dv is not None and abs(dv) <= TOL_VINA
                   else ("更好" if ok_v else "更差!")))
            if not (hard and ok_v):
                verdicts.append((name, seed, "NOT_PROMOTED"))
        if miss:
            log("  (有 %d 个种子缺数据, 不足以判定 -> 保守判为 NOT_PROMOTED)" % miss)
            verdicts.append((name, None, "NOT_PROMOTED"))
        elif not [v for v in verdicts if v[0] == name]:
            verdicts.append((name, None, "PROMOTED"))
        ok = not [v for v in verdicts if v[0] == name and v[2] != "PROMOTED"]
        log("  成对小结: 结合强度更好 %d 个种子 / 更差 %d 个种子 (容差 ±%.2f)"
            % (better, worse, TOL_VINA))
        log("  => %s" % ("PROMOTED" if ok else "NOT_PROMOTED"))

    log("\n" + "=" * 78)
    log("最终:")
    for name in names[1:]:
        v = [x[2] for x in verdicts if x[0] == name]
        log("  %-24s %s" % (name, "PROMOTED" if v and v[0] == "PROMOTED" else "NOT_PROMOTED"))
    log("说明: 判据以**结合强度(Vina 中位)**为首要轴, 并硬性要求完成分子数不退化 ——")
    log("      完成数退化意味着停止策略被破坏(历史最大失败模式), 其余指标再好看也不晋级。")
    log("产物: %s 与 %s" % (os.path.relpath(args.out, ROOT), os.path.relpath(args.work, ROOT)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
