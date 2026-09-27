# -*- coding: utf-8 -*-
"""从既有采样产物计算「尺寸 / 类药性 / strict 通过率」口径, 供 v9 预登记验收判据使用。

为什么需要它
------------
记忆文件 §5.2e 为 v9 预登记了两条硬指标, 而 `select_ckpt_by_generation.py`
与 `ab_docking_compare.py` 都不直接给出这两条:
  1. 生成分子的 MW/重原子中位落在训练分布附近(**MW>=340 / 重原子>=25**);
  2. 药库 **strict 过滤通过率不得低于官方基线的 90%**。
本脚本只做"读既有采样产物 + 算口径", **不重新生成、不重新对接**,
因此可以在任何 A/B 跑完后离线复核, 也便于把同样的尺子套到官方/v3/v9 各臂上。

口径来源(单一真源, 不另立标准)
------------------------------
- strict 档位: `build_gpcr_v3.py` 的 `DRUGLIKE_TIERS['strict']`
  = MW∈[250,500] / LogP≤5.0 / QED≥0.40 / SA≤6 / 环数≥1
- 计算实现直接 `import` 自 `sweep_frontier_threshold.py` 的 `metrics()`,
  避免第三份重复实现导致口径漂移。

用法
----
    python src/scripts/arm_gen_metrics.py --arms official=outputs/ab_docking/<name> \\
                                                v9=outputs/sel_v9/7500 \\
                                                v3=outputs/sel_v3
    # 也接受直接指向 SMILES.txt
    python src/scripts/arm_gen_metrics.py --arms v9=outputs/sel_v9/7500/samples/<session>/SMILES.txt

每个 arm 的参数是 `名字=路径`; 路径下会递归找**最新的**一个 `SMILES.txt`
(采样会话目录按时间戳命名, 取最新的即可)。
"""
import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, OUTPUTS, ensure_dir  # noqa: E402

SIZE_MW_MIN = 340.0   # §5.2e 预登记: MW 中位下限
SIZE_HA_MIN = 25.0    # §5.2e 预登记: 重原子中位下限
STRICT_RATIO_MIN = 0.90  # §5.2e 预登记: strict 通过率 >= 官方基线的 90%


def find_smiles(path):
    """返回该路径下最新的 SMILES.txt。"""
    if os.path.isfile(path):
        return path
    cands = glob.glob(os.path.join(path, "**", "SMILES.txt"), recursive=True)
    if not cands:
        return None
    # 会话目录名带时间戳, 按目录名排序即按时间排序
    return sorted(cands)[-1]


def main():
    ap = argparse.ArgumentParser(description="既有采样产物的尺寸/类德行口径")
    ap.add_argument("--arms", nargs="+", required=True, help="名字=路径 (至少两个)")
    ap.add_argument("--baseline", default=None,
                    help="作为 strict 通过率基准的 arm 名; 默认第一个")
    ap.add_argument("--out", default=os.path.join(OUTPUTS, "arm_gen_metrics.csv"))
    args = ap.parse_args()

    sys.path.insert(0, os.path.join(ROOT, "src", "scripts"))
    from sweep_frontier_threshold import metrics  # 单一真源的口径实现
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    rows = []
    for spec in args.arms:
        if "=" not in spec:
            raise SystemExit("arm 参数必须是 名字=路径: %s" % spec)
        name, path = spec.split("=", 1)
        sp = find_smiles(path)
        if not sp:
            print("[跳过] %s: 在 %s 下找不到 SMILES.txt" % (name, path))
            continue
        smis = [s.strip() for s in open(sp, encoding="utf-8", errors="ignore") if s.strip()]
        m = metrics(smis)
        if not m:
            print("[跳过] %s: SMILES 全部无法解析" % name)
            continue
        mols = [x for x in (Chem.MolFromSmiles(s) for s in smis) if x is not None]
        n = len(mols)
        n_mw = sum(1 for x in mols if x.GetNumHeavyAtoms() >= 1)
        rows.append(dict(
            arm=name, smiles_file=os.path.relpath(sp, ROOT), n_finished=len(smis), n_parsed=n,
            mw_med=round(m["mw"], 1), ha_med=round(m["ha"], 1),
            qed_med=round(m["qed"], 4), sa_med=round(m["sa"], 3), n_scaffold=m["n_scaf"],
            # sweep_frontier_threshold.metrics 只返回 pass_pct(百分比), 这里换回计数
            strict_pass=round(m["pass_pct"] * n / 100.0),
            strict_rate=round(m["pass_pct"] / 100.0, 4),
            frac_mw_ge340=round(sum(1 for x in mols
                                    if _mw(x) >= SIZE_MW_MIN) / n, 4),
            frac_ha_ge25=round(sum(1 for x in mols
                                   if x.GetNumHeavyAtoms() >= SIZE_HA_MIN) / n, 4),
            _n_mw=n_mw))
        print("[%s] n=%d | MW中位 %.1f | 重原子中位 %.1f | QED %.3f | SA %.2f | 骨架 %d | "
              "strict %d (%.1f%%) | MW>=340 %.0f%% | HA>=25 %.0f%%"
              % (name, len(smis), m["mw"], m["ha"], m["qed"], m["sa"], m["n_scaf"],
                 round(m["pass_pct"] * n / 100.0), 100.0 * round(m["pass_pct"] * n / 100.0) / n,
                 100.0 * sum(1 for x in mols if _mw(x) >= SIZE_MW_MIN) / n,
                 100.0 * sum(1 for x in mols if x.GetNumHeavyAtoms() >= SIZE_HA_MIN) / n))

    if not rows:
        raise SystemExit("没有任何可用 arm")

    base = args.baseline or rows[0]["arm"]
    b = next((r for r in rows if r["arm"] == base), rows[0])
    print("\n" + "=" * 96)
    print("预登记判据复核 (基准 arm = %s)" % b["arm"])
    print("=" * 96)
    verdicts = {}
    for r in rows:
        ok_size = (r["mw_med"] >= SIZE_MW_MIN) and (r["ha_med"] >= SIZE_HA_MIN)
        ratio = r["strict_rate"] / b["strict_rate"] if b["strict_rate"] else None
        ok_strict = (ratio is None) or (ratio >= STRICT_RATIO_MIN)
        verdicts[r["arm"]] = dict(size_ok=bool(ok_size), strict_ratio=ratio,
                                  strict_ok=bool(ok_strict))
        print("%-10s MW中位 %6.1f %s | 重原子中位 %5.1f %s | strict通过率 %5.1f%% "
              "(占基准 %s) %s"
              % (r["arm"], r["mw_med"], "OK" if r["mw_med"] >= SIZE_MW_MIN else "**低**",
                 r["ha_med"], "OK" if r["ha_med"] >= SIZE_HA_MIN else "**低**",
                 100.0 * r["strict_rate"],
                 ("%.0f%%" % (100 * ratio)) if ratio is not None else "NA",
                 "OK" if ok_strict else "**不足 90%**"))
    ensure_dir(os.path.dirname(args.out))
    with open(args.out, "w", encoding="utf-8-sig", newline="") as f:
        import csv
        w = csv.DictWriter(f, fieldnames=[k for k in rows[0].keys()])
        w.writeheader()
        for r in rows:
            w.writerow({k: v for k, v in r.items() if k in w.fieldnames})
    print("\n判据结论(JSON): %s" % json.dumps(verdicts, ensure_ascii=False))
    print("明细: %s" % os.path.relpath(args.out, ROOT))
    return 0


def _mw(mol):
    from rdkit.Chem import Descriptors
    return Descriptors.MolWt(mol)


if __name__ == "__main__":
    sys.exit(main())
