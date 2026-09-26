# -*- coding: utf-8 -*-
"""对一个已有的化合物库目录跑 AutoDock Vina 对接, 并写出两种结果文件。

为什么单独成脚本: docking_pipeline.py 是**单进程串行**的(实测约 1.7 个/分钟),
直接对几百个分子跑要几小时。本脚本把配体切成 N 块并发调用它, 再合并结果,
实测 521 个分子约 7 分钟。

产物(都写在库目录下, 供 predict.py 的 ENRICH_SPECS 自动并入清单):
    docking_summary.csv   smiles, vina_kcal, receptor, center, exhaustiveness
    top_candidates.csv    rank, 综合分, vina_kcal, qed, sa, mw, murcko, smiles
                          (predict.py 读它的 vina_kcal 列 -> 清单的对接分)

用法:
    python src/scripts/dock_library.py --library results --target A2A
    python src/scripts/dock_library.py --library <lib> --target B2AR --parallel 6
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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, CONFIGS, OUTPUTS, ensure_dir  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="对化合物库跑并行对接")
    ap.add_argument("--library", required=True, help="含 compounds.csv 的库目录")
    ap.add_argument("--target", required=True, help="靶点名(取自 configs/targets.json)")
    ap.add_argument("--parallel", type=int, default=6, help="并发分块数")
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--work", default=None, help="中间目录, 默认 outputs/dock_lib_<target>")
    args = ap.parse_args()

    lib = os.path.abspath(args.library)
    comp = os.path.join(lib, "compounds.csv")
    if not os.path.exists(comp):
        raise SystemExit("缺少 %s" % comp)
    reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    if args.target not in reg:
        raise SystemExit("未知靶点 %s; 可选 %s" % (args.target, sorted(reg)))
    t = reg[args.target]
    pdb = os.path.join(ROOT, t["pdb"])
    center = ",".join("%.3f" % v for v in t["center"])

    rows = list(csv.DictReader(open(comp, encoding="utf-8-sig")))
    smis = [r["smiles"] for r in rows if r.get("smiles")]
    print("靶点 %s | 受体 %s | 中心 %s | 待对接 %d 个分子"
          % (args.target, os.path.basename(pdb), center, len(smis)), flush=True)

    work = args.work or os.path.join(OUTPUTS, "dock_lib_%s" % args.target)
    shutil.rmtree(work, ignore_errors=True)
    ensure_dir(work)

    n = max(1, min(args.parallel, len(smis)))
    chunks = [[] for _ in range(n)]
    for i, s in enumerate(smis):
        chunks[i % n].append(s)

    procs = []
    for ci, ch in enumerate(chunks):
        if not ch:
            continue
        cf = os.path.join(work, "smiles_p%d.txt" % ci)
        with open(cf, "w", encoding="ascii") as f:
            f.write("\n".join(ch) + "\n")
        cdir = os.path.join(work, "p%d" % ci)
        cmd = [sys.executable, os.path.join(ROOT, "src", "scripts", "docking_pipeline.py"),
               "--receptor", pdb, "--center=" + center, "--smiles-file", cf,
               "--out", cdir, "--exhaustiveness", str(args.exhaustiveness)]
        lf = open(os.path.join(work, "p%d.log" % ci), "w", encoding="utf-8", errors="ignore")
        procs.append((subprocess.Popen(cmd, cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT),
                      cdir, lf))
        print("  分块 %d: %d 个分子" % (ci, len(ch)), flush=True)

    for p, cdir, lf in procs:
        p.wait()
        lf.close()

    best = {}
    for p, cdir, lf in procs:
        cs = os.path.join(cdir, "summary.csv")
        if not os.path.exists(cs):
            continue
        for r in csv.DictReader(open(cs, encoding="utf-8-sig")):
            if r.get("status") != "ok":
                continue
            v = float(r["vina_score"])
            s = r["smiles"]
            if s not in best or v < best[s]:
                best[s] = v
    print("对接成功(去重): %d / %d" % (len(best), len(smis)), flush=True)

    with open(os.path.join(lib, "docking_summary.csv"), "w",
              encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["smiles", "vina_kcal", "receptor", "center", "exhaustiveness"])
        for s, v in sorted(best.items(), key=lambda x: x[1]):
            w.writerow([s, round(v, 3), os.path.basename(pdb), center, args.exhaustiveness])

    by = {r["smiles"]: r for r in rows}
    with open(os.path.join(lib, "top_candidates.csv"), "w",
              encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "综合分", "vina_kcal", "qed", "sa", "mw", "murcko", "smiles"])
        for i, (s, v) in enumerate(sorted(best.items(), key=lambda x: x[1]), 1):
            r = by.get(s, {})
            w.writerow([i, "", round(v, 3), r.get("qed", ""), r.get("sa", ""),
                        r.get("mw", ""), r.get("murcko_scaffold", ""), s])

    vs = sorted(best.values())
    if vs:
        print("Vina: n=%d 中位 %.2f 最强 %.2f 最弱 %.2f | 强于 -10 的 %d (%.1f%%)"
              % (len(vs), st.median(vs), vs[0], vs[-1],
                 sum(1 for x in vs if x <= -10),
                 100.0 * sum(1 for x in vs if x <= -10) / len(vs)), flush=True)
    print("已写: docking_summary.csv, top_candidates.csv", flush=True)


if __name__ == "__main__":
    main()
