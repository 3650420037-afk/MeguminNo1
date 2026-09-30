# -*- coding: utf-8 -*-
"""汇总新库（outputs/library_v10/）的关键数字，供 results/README 与文档直接引用。

输出：每靶点的入库数、QED/SA/MW 中位、唯一 Murcko 骨架、对接中位/最强/强于−10 占比、
top-50 对接中位（与 README 既有口径对齐）、ADMET 分级（若有）。
"""
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
from rdkit.Chem import Descriptors, QED  # noqa: E402
from rdkit.Chem.Scaffolds import MurckoScaffold  # noqa: E402

RDLogger.DisableLog("rdApp.*")
STAGE = os.path.join(ROOT, "outputs", "library_v10")


def load(path):
    out = []
    if not os.path.exists(path):
        return out
    with io.open(path, encoding="utf-8-sig", errors="ignore") as fh:
        rd = csv.DictReader(fh)
        key = next((k for k in rd.fieldnames if "smiles" in k.lower()), None)
        if not key:
            return out
        for r in rd:
            m = Chem.MolFromSmiles((r.get(key) or "").strip())
            if m:
                out.append(m)
    return out


def dock(target, lib):
    p = os.path.join(lib, "docking_summary.csv")
    if not os.path.exists(p):
        return None
    v = []
    for r in csv.DictReader(io.open(p, encoding="utf-8-sig")):
        k = next((c for c in r if "vina" in c.lower()), None)
        try:
            v.append(float(r[k]))
        except Exception:
            pass
    if not v:
        return None
    v.sort()
    top50 = st.median(v[:50]) if len(v) >= 50 else None
    return dict(n=len(v), med=st.median(v), best=min(v),
                p10=sum(1 for x in v if x <= -10) / len(v), top50=top50)


def admet(target):
    for cand in (os.path.join(STAGE, target, "library", "admet_screen.csv"),
                 os.path.join(STAGE, target, "results", "admet_screen.csv")):
        if os.path.exists(cand):
            c = {}
            for r in csv.DictReader(io.open(cand, encoding="utf-8-sig")):
                k = next((x for x in r if "ADMET" in x or "分级" in x), None)
                if k:
                    c[r[k]] = c.get(r[k], 0) + 1
            return c
    return None


def main():
    print("%-6s %5s %7s %6s %6s %7s %8s %8s %8s %8s %s"
          % ("靶点", "入库", "MW中位", "HA", "QED", "SA", "骨架", "对接中位", "top50", "≤-10占比", "ADMET"))
    for t in ("A2A", "B2AR", "D3", "5HT2B"):
        lib = os.path.join(STAGE, t, "library")
        mols = load(os.path.join(lib, "compounds.csv"))
        if not mols:
            print("%-6s （缺）" % t)
            continue
        mw = st.median(Descriptors.MolWt(m) for m in mols)
        ha = st.median(m.GetNumHeavyAtoms() for m in mols)
        q = st.median(QED.qed(m) for m in mols)
        scaf = len({MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols})
        d = dock(t, lib)
        a = admet(t)
        print("%-6s %5d %7.1f %6.1f %6.3f %7s %8d %8s %8s %8s %s"
              % (t, len(mols), mw, ha, q, "-", scaf,
                 ("%.2f" % d["med"]) if d else "-",
                 ("%.2f" % d["top50"]) if d and d["top50"] else "-",
                 ("%.1f%%" % (100 * d["p10"])) if d else "-",
                 ("A%d/B%d/C%d" % (a.get("A", 0), a.get("B", 0), a.get("C", 0))) if a else "-"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
