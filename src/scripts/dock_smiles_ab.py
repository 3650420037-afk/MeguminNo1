# -*- coding: utf-8 -*-
"""按**对接分**比较任意两组分子(不重新生成)。

为什么需要它
------------
`ab_docking_compare.py` 内置了生成步骤, 用于比较**检查点**; 但很多实验产生的是
**已经采好的分子集合**(例如 `sweep_frontier_threshold.py` 各阈值档的产物、
不同引导参数下的产物)。要判断"这组分子是否真的更亲和", 必须能直接把它们送进
同一条 Vina 管线做受控比较, 而不是重新采样一遍。

纪律(与本仓库其它 A/B 工具一致)
-------------------------------
- **同时报原始 Vina 与配体效率 LE**, 因为 Vina 打分与分子大小强相关;
  两臂尺寸差超过 15% 时显式警告, 防止把"分子更小"误读成"亲和力更差"。
- 硬性要求**分子数不退化**(不再生成, 故比对输入分子数)。
- 复用 ab_docking_compare.dock_smiles, 保证对接协议与其它实验完全一致。

用法
----
    python src/scripts/dock_smiles_ab.py \
        --arms "阈值0=outputs/thr_sweep_v3/thr0p0/samples/*/SMILES.txt" \
               "阈值-0.25=outputs/thr_fine_v3/thr-0p25/samples/*/SMILES.txt" \
        --target A2A --parallel 6
"""
import argparse
import glob
import os
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, ensure_dir  # noqa: E402
from scripts.ab_docking_compare import dock_smiles  # noqa: E402

TOL_VINA = 0.20


def log(m):
    print(m, flush=True)


def main():
    ap = argparse.ArgumentParser(description="按对接分比较任意两组分子")
    ap.add_argument("--arms", nargs="+", required=True,
                    help="每项形如 名称=glob路径(可含通配); 第一项作为 baseline")
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--parallel", type=int, default=6)
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--work", default=os.path.join(ROOT, "outputs", "dock_smiles_ab"))
    args = ap.parse_args()

    parsed = []
    for spec in args.arms:
        if "=" not in spec:
            raise SystemExit("--arms 每项必须是 名称=路径, 收到: %s" % spec)
        name, pat = spec.split("=", 1)
        smis = []
        for f in glob.glob(os.path.normpath(pat)):
            smis += [l.strip() for l in open(f, encoding="utf-8") if l.strip()]
        if not smis:
            raise SystemExit("该臂没有分子: %s (%s)" % (name, pat))
        parsed.append((name, smis))
        log("[%s] %d 个分子" % (name, len(smis)))

    rows = []
    for name, smis in parsed:
        wd = ensure_dir(os.path.join(args.work, name.replace("/", "_").replace(" ", "_")))
        log("\n[%s] 对接中 (%d 个) ..." % (name, len(smis)))
        pairs, n_ok = dock_smiles(smis, args.target, os.path.join(wd, "dock"),
                                  args.parallel, args.exhaustiveness,
                                  os.path.join(wd, "dock.log"))
        if not pairs:
            log("  没有可用结果")
            continue
        vs = [v for v, _h in pairs]
        has = [h for _v, h in pairs]
        le = [-v / h for v, h in pairs]
        rows.append(dict(name=name, n_in=len(smis), n_docked=n_ok,
                         vina_med=st.median(vs), ha_med=st.median(has), le_med=st.median(le),
                         vina_min=min(vs),
                         frac_le10=100.0 * sum(1 for x in vs if x <= -10) / len(vs)))
        log("  ok %d | Vina中位 %.2f | 重原子中位 %.1f | LE中位 %.3f | <=-10 %.1f%%"
            % (n_ok, rows[-1]["vina_med"], rows[-1]["ha_med"], rows[-1]["le_med"],
               rows[-1]["frac_le10"]))

    if not rows:
        raise SystemExit("没有任何结果")

    log("\n" + "=" * 78)
    log("结果 (baseline = %s)" % rows[0]["name"])
    log("=" * 78)
    log("%-16s %5s %5s %9s %7s %8s %8s" %
        ("臂", "输入", "对接", "Vina中位", "重原子", "LE", "<=-10%"))
    for r in rows:
        log("%-16s %5d %5d %9.2f %7.1f %8.3f %7.1f%%" %
            (r["name"], r["n_in"], r["n_docked"], r["vina_med"], r["ha_med"],
             r["le_med"], r["frac_le10"]))

    base = rows[0]
    log("\n--- 与 baseline 对比 ---")
    log("⚠ Vina 打分与分子大小强相关, 必须同时看 LE; 尺寸差 >15% 时原始分不可直接比。")
    for r in rows[1:]:
        dn = r["n_in"] - base["n_in"]
        dv = r["vina_med"] - base["vina_med"]
        dle = r["le_med"] - base["le_med"]
        ratio = r["ha_med"] / base["ha_med"] if base["ha_med"] else None
        log("  %-16s 分子数 %+d (%s) | 原始Vina %+.2f (%s) | LE %+.3f | 尺寸比 %s%s"
            % (r["name"], dn, "OK" if dn >= -5 else "退化!", dv,
               "更好" if dv <= -TOL_VINA else ("打平" if dv <= TOL_VINA else "更差!"),
               dle, ("%.2f" % ratio) if ratio else "NA",
               "  ⚠尺寸不可比" if (ratio and (ratio < 0.85 or ratio > 1.15)) else ""))
    out = os.path.join(args.work, "dock_smiles_ab.csv")
    import csv as _csv
    with open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    log("\n产物: %s" % os.path.relpath(out, ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
