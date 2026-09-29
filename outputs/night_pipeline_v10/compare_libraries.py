# -*- coding: utf-8 -*-
"""新旧药库对比：读**旧库**（results/，2026-09 由 v2 权重生成、已备份）与**新库**（outputs/library_v10/，v10 生成），
输出一页对照表（分子数 / QED / SA / MW / 重原子 / strict 通过率 / 对接中位 / 唯一骨架），供 results/README、ModelCard 与附件3 引用。

用法：python outputs/night_pipeline_v10/compare_libraries.py [--targets A2A B2AR D3 5HT2B]
"""
import argparse
import csv
import glob
import io
import os
import statistics as st
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import Descriptors, QED, rdMolDescriptors  # noqa: E402
from rdkit.Chem.Scaffolds import MurckoScaffold  # noqa: E402

RDLogger.DisableLog("rdApp.*")
STRICT = dict(mw=(250.0, 500.0), logp_max=5.0, qed_min=0.40, sa_max=6.0)


def find_col(header, *keys):
    for i, h in enumerate(header):
        for k in keys:
            if k.lower() in h.lower():
                return i
    return None


def load_smiles(path):
    out = []
    if not os.path.exists(path):
        return out
    with io.open(path, encoding="utf-8-sig", errors="ignore") as fh:
        rd = csv.reader(fh)
        header = next(rd, None)
        if header is None:
            return out
        ci = find_col(header, "smiles", "SMILES")
        if ci is None:
            return out
        for row in rd:
            if len(row) > ci and row[ci].strip():
                m = Chem.MolFromSmiles(row[ci].strip())
                if m:
                    out.append(m)
    return out


def strict_rate(mols, sa):
    n = 0
    for m, s in zip(mols, sa):
        mw = Descriptors.MolWt(m)
        if (STRICT["mw"][0] <= mw <= STRICT["mw"][1] and Descriptors.MolLogP(m) <= STRICT["logp_max"]
                and QED.qed(m) >= STRICT["qed_min"] and s <= STRICT["sa_max"]
                and rdMolDescriptors.CalcNumRings(m) >= 1):
            n += 1
    return n / len(mols) if mols else float("nan")


def summarise(mols, sa):
    if not mols:
        return None
    mw = sorted(Descriptors.MolWt(m) for m in mols)
    ha = sorted(m.GetNumHeavyAtoms() for m in mols)
    qs = sorted(QED.qed(m) for m in mols)
    # 注意：Murcko 骨架要用 rdkit.Chem.Scaffolds.MurckoScaffold（rdMolDescriptors 里没有该函数，自测踩到）
    scaf = {MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols}
    return dict(n=len(mols), mw=st.median(mw), ha=st.median(ha), qed=st.median(qs),
                sa=st.median(sa), scaf=len(scaf), strict=strict_rate(mols, sa))


def sa_of(mols):
    from rdkit.Chem import RDConfig
    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer
    return [sascorer.calculateScore(m) for m in mols]


def doc_median(csv_path):
    """从对接结果里取 Vina 中位（兼容列名差异）。"""
    if not os.path.exists(csv_path):
        return None
    vals = []
    for r in csv.DictReader(io.open(csv_path, encoding="utf-8-sig")):
        for k, v in r.items():
            if "vina" in k.lower() and v not in (None, "", "nan"):
                try:
                    vals.append(float(v))
                except Exception:
                    pass
                break
    return (st.median(vals), len(vals)) if vals else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", default=["A2A", "B2AR", "D3", "5HT2B"])
    args = ap.parse_args()
    print("%-6s %-10s %6s %8s %7s %7s %7s %8s %9s" %
          ("靶点", "库", "n", "MW中位", "HA", "QED", "SA", "strict", "对接中位"))
    for t in args.targets:
        # A2A 旧库在 results/ 根目录（results/compounds.csv），其余靶点在 results/<T>/
        old_a2a = os.path.join(ROOT, "results", "compounds.csv")
        old_t = os.path.join(ROOT, "results", t, "compounds.csv")
        old = old_t if os.path.exists(old_t) else old_a2a
        pairs = [("旧库(v2)", old),
                 ("新库(v10)", os.path.join(ROOT, "outputs", "library_v10", t, "library", "compounds.csv"))]
        for tag, path in pairs:
            mols = load_smiles(path)
            if not mols:
                print("%-6s %-10s %6s" % (t, tag, "（缺）"))
                continue
            m = summarise(mols, sa_of(mols))
            d = doc_median(os.path.join(os.path.dirname(path), "docking_summary.csv"))
            print("%-6s %-10s %6d %8.1f %7.1f %7.3f %7.2f %7.3f %9s"
                  % (t, tag, m["n"], m["mw"], m["ha"], m["qed"], m["sa"], m["strict"],
                     ("%.2f (n=%d)" % d) if d else "（未对接）"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
