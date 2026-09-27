# -*- coding: utf-8 -*-
"""夜间模式: 把 4 个靶点的全量 ChEMBL 活性分子对接进各自实验口袋, 建成训练集。

为什么需要它
------------
已实测证明「每口袋样本密度」是微调成败的决定因素(1.04/口袋 -> 分子永不终止;
40-127/口袋 -> 完成数 61-64)。此前每靶点只有 99-191 个活性分子, 现已抓到
17,177 个(data/ligands/*_活性配体_full.csv), 本脚本把它们对接成 (口袋, 位姿) 训练对。

可断点续跑
----------
按 (靶点, 分块号) 落盘 summary.csv 与位姿文件; 重跑时**已完成的分块直接跳过**,
因此夜间中断/断网后可原样续跑。

用法
----
    # 1) 全量对接(可断点续跑)
    python src/scripts/mass_dock_actives.py --targets A2A B2AR D3 5HT2B --parallel 12
    # 2) 布局转换 + 校验(mass_dock -> build_docked_dataset 期望的布局)
    python src/scripts/adapt_mass_dock.py --src outputs/mass_dock --dst outputs/mass_dock_build
    # 3) 建训练集
    python src/scripts/build_docked_dataset.py --work outputs/mass_dock_build \
        --reuse-dock --out data/gpcr_mass_v1

注意: 本脚本**只负责对接**, 不建数据集。分块是 chunks[i % n], 因此换 --parallel 或
--cap 会整体改变分子到块的映射; 续跑时会逐块比对当时排队的分子名单, 不一致就归档重跑,
避免"静默少对接一批分子"。
"""
import argparse
import csv
import glob
import json
import os
import shutil
import statistics as st
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, CONFIGS, DATA, LIGANDS, ensure_dir  # noqa: E402

SUPPORTED = {6, 7, 8, 9, 15, 16, 17}     # C N O F P S Cl
MW_LO, MW_HI, LOGP_MAX, RINGS_MIN = 200.0, 600.0, 6.0, 2
WORK_ROOT = os.path.join(ROOT, "outputs", "mass_dock")


def log(m):
    print(m, flush=True)


def load_queue(target, cap):
    """读取该靶点的全量活性分子, 按词表/类药性过滤, 按效价降序取前 cap 个。"""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    p = os.path.join(LIGANDS, "%s_活性配体_full.csv" % target)
    if not os.path.exists(p):
        return []
    rows = []
    with open(p, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            smi = (r.get("canonical_smiles") or "").strip()
            if not smi:
                continue
            try:
                pv = float(r.get("pchembl_value") or 0)
            except ValueError:
                pv = 0.0
            rows.append((pv, smi))
    rows.sort(key=lambda x: -x[0])
    kept, seen, rej = [], set(), {"elem": 0, "mw": 0, "logp": 0, "rings": 0, "bad": 0}
    for pv, smi in rows:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            rej["bad"] += 1
            continue
        if {a.GetAtomicNum() for a in m.GetAtoms()} - SUPPORTED:
            rej["elem"] += 1
            continue
        if not (MW_LO <= Descriptors.MolWt(m) <= MW_HI):
            rej["mw"] += 1
            continue
        if Crippen.MolLogP(m) > LOGP_MAX:
            rej["logp"] += 1
            continue
        if rdMolDescriptors.CalcNumRings(m) < RINGS_MIN:
            rej["rings"] += 1
            continue
        c = Chem.MolToSmiles(m)
        if c in seen:
            continue
        seen.add(c)
        kept.append(c)
        if cap and len(kept) >= cap:
            break
    log("  [%s] 候选 %d -> 过滤后 %d (淘汰: %s)" % (target, len(rows), len(kept), rej))
    return kept


def run_target(target, args):
    reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    t = reg[target]
    pdb = os.path.join(ROOT, t["pdb"])
    center = ",".join("%.3f" % v for v in t["center"])
    work = ensure_dir(os.path.join(WORK_ROOT, target))

    smiles = load_queue(target, args.cap)
    if not smiles:
        log("  [%s] 无可用分子, 跳过" % target)
        return None
    n = max(1, args.parallel)
    chunks = [[] for _ in range(n)]
    for i, s in enumerate(smiles):
        chunks[i % n].append(s)

    # ---- 续跑一致性校验(必须) ----------------------------------------------
    # 分块是 chunks[i % n], 所以 **n(或 cap/队列)一变, 分子到块的映射整体改变**。
    # 若仅凭 summary.csv 是否存在来判断"这块已完成", 换过 --parallel/--cap 之后就会
    # 把旧映射的块误判为新映射的完成块而跳过。
    # 实测后果: 先用 --parallel 6 冒烟跑 300 个(每块 50), 再用 --parallel 12 全量跑,
    # 日志出现「共 12 块, 已完成 6 块, 待跑 6 块」—— 于是第 300~3000 号活性分子
    # **永远不会被对接**(静默丢数据, 且后续建库不会报错)。
    # 因此判断"已完成"的必要条件是: 该块当时排队的分子名单与当前**逐字一致**。
    todo, stale = [], []
    for ci, ch in enumerate(chunks):
        if not ch:
            continue
        cdir = os.path.join(work, "c%02d" % ci)
        sf = os.path.join(work, "s%02d.txt" % ci)
        has_sum = os.path.exists(os.path.join(cdir, "summary.csv"))
        same = False
        if os.path.exists(sf):
            try:
                with open(sf, encoding="ascii") as f:
                    same = f.read().split() == ch
            except Exception:
                same = False
        if has_sum and not same:
            stale.append(ci)
        todo.append((ci, ch, cdir, has_sum and same))

    if stale:
        stamp = time.strftime("%Y%m%d_%H%M%S")
        arch = ensure_dir(os.path.join(work, "_stale_" + stamp))
        for ci in stale:
            for p in ([os.path.join(work, "c%02d" % ci), os.path.join(work, "s%02d.txt" % ci)]):
                if os.path.exists(p):
                    shutil.move(p, os.path.join(arch, os.path.basename(p)))
            log("  [%s] 块 %02d 的分子名单与当前分块不一致 -> 已归档到 _stale_%s, 将重新对接"
                % (target, ci, stamp))
        # 顶层合并产物同样是基于旧分块生成的, 一并归档, 避免残留旧条目
        for p in (os.path.join(work, "docking"), os.path.join(work, "summary_merged.csv")):
            if os.path.exists(p):
                shutil.move(p, os.path.join(arch, os.path.basename(p)))

    # 记录本次运行的参数, 便于追溯"这份产物是用什么 cap/并发跑出来的"
    with open(os.path.join(work, "run_info.json"), "w", encoding="utf-8") as f:
        json.dump({"target": target, "cap": int(args.cap), "parallel": n,
                   "exhaustiveness": int(args.exhaustiveness),
                   "n_queue": len(smiles),
                   "chunk_sizes": [len(c) for c in chunks]},
                  f, ensure_ascii=False, indent=1)

    ndone = sum(1 for *_, d in todo if d)
    log("  [%s] 共 %d 块, 已完成 %d 块, 待跑 %d 块" % (target, len(todo), ndone, len(todo) - ndone))

    running = []
    for ci, ch, cdir, done in todo:
        if done:
            continue
        cf = os.path.join(work, "s%02d.txt" % ci)
        with open(cf, "w", encoding="ascii") as f:
            f.write("\n".join(ch) + "\n")
        cmd = [sys.executable, os.path.join(ROOT, "src", "scripts", "docking_pipeline.py"),
               "--receptor", pdb, "--center=" + center, "--smiles-file", cf,
               "--out", cdir, "--exhaustiveness", str(args.exhaustiveness)]
        lf = open(os.path.join(work, "c%02d.log" % ci), "w", encoding="utf-8", errors="ignore")
        running.append((ci, subprocess.Popen(cmd, cwd=ROOT, stdout=lf,
                                             stderr=subprocess.STDOUT), lf))
    log("  [%s] 已启动 %d 个对接进程 (每块 %d 个分子)"
        % (target, len(running), len(chunks[0]) if chunks else 0))
    t0 = time.time()
    while running:
        alive = []
        for ci, p, lf in running:
            if p.poll() is None:
                alive.append((ci, p, lf))
            else:
                lf.close()
                log("  [%s] 块 %02d 完成 (%.0f 秒)"
                    % (target, ci, time.time() - t0))
        running = alive
        if running:
            time.sleep(30)
    log("  [%s] 全部对接完成, 用时 %.1f 分钟" % (target, (time.time() - t0) / 60))

    # 合并位姿文件到统一目录, 供建库阶段使用
    alldock = ensure_dir(os.path.join(work, "docking"))
    for ci, ch, cdir, done in todo:
        for pq in glob.glob(os.path.join(cdir, "docking", "*_out.pdbqt")):
            base = "c%02d_%s" % (ci, os.path.basename(pq))
            dst = os.path.join(alldock, base)
            if not os.path.exists(dst):
                shutil.copyfile(pq, dst)
    # 合并 summary
    merged = []
    for ci, ch, cdir, done in todo:
        cs = os.path.join(cdir, "summary.csv")
        if not os.path.exists(cs):
            continue
        for r in csv.DictReader(open(cs, encoding="utf-8-sig")):
            r = dict(r)
            r["name"] = "c%02d_%s" % (ci, r.get("name", ""))
            merged.append(r)
    with open(os.path.join(work, "summary_merged.csv"), "w",
              encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "smiles", "vina_score", "status"])
        w.writeheader()
        w.writerows(merged)
    ok = [r for r in merged if r.get("status") == "ok"]
    vs = sorted(float(r["vina_score"]) for r in ok if r.get("vina_score"))
    log("  [%s] 汇总: %d 条, ok %d%s"
        % (target, len(merged), len(ok),
           (", Vina 中位 %.2f 最强 %.2f 强于-10 %.1f%%"
            % (st.median(vs), vs[0], 100.0 * sum(1 for x in vs if x <= -10) / len(vs)))
           if vs else ""))
    return work


def main():
    ap = argparse.ArgumentParser(description="夜间模式: 全量活性分子对接")
    ap.add_argument("--targets", nargs="+", default=["A2A", "B2AR", "D3", "5HT2B"])
    ap.add_argument("--parallel", type=int, default=12, help="每靶点并发分块数")
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--cap", type=int, default=6000, help="每靶点最多对接多少分子(按效价降序)")
    args = ap.parse_args()

    ensure_dir(WORK_ROOT)
    log("=" * 74)
    log("夜间模式: 全量活性分子对接 (并发 %d/靶点, 每靶点上限 %d)"
        % (args.parallel, args.cap))
    log("=" * 74)
    for tg in args.targets:
        log("[%s] 开始" % tg)
        run_target(tg, args)
        log("[%s] 结束" % tg)
    log("=" * 74)
    log("全部靶点对接完成。产物: outputs/mass_dock/<靶点>/summary_merged.csv + docking/")


if __name__ == "__main__":
    main()
