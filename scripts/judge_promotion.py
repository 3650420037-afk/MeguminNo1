# -*- coding: utf-8 -*-
"""微调模型晋级判定: 对比 baseline vs finetuned 的 finished/QED/SA/骨架多样性."""
import sys, os, glob
sys.path.insert(0, r"D:\MMModel\Pocket2Mol")
import numpy as np
from rdkit import RDLogger, Chem
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import QED, Descriptors, Crippen
import utils.guidance  # rdkit.six shim
from evaluation.sascorer import calculateScore

def load_run(path):
    fs = glob.glob(os.path.join(path, "**", "SMILES.txt"), recursive=True)
    if not fs:
        return None
    smis = sorted(set(s for s in open(fs[0], encoding="utf-8").read().splitlines() if s.strip()))
    ms = [m for m in (Chem.MolFromSmiles(s) for s in smis) if m]
    if not ms:
        return None
    return {
        "n": len(smis),
        "qed": np.mean([QED.qed(m) for m in ms]),
        "sa": np.mean([calculateScore(m) for m in ms]),
        "scaf": len({MurckoScaffold.MurckoScaffoldSmiles(smiles=s) for s in smis}),
    }

if __name__ == "__main__":
    b = load_run(sys.argv[1])
    f = load_run(sys.argv[2])
    out = sys.argv[3]
    lines = []
    lines.append("baseline: %s" % b)
    lines.append("finetuned: %s" % f)
    promoted = False
    if b and f:
        # 计划表判据: 有效率/合法率/QED/SA/多样性不低于基线(允许小幅容差)
        promoted = (f["n"] >= b["n"] - 5) and (f["qed"] >= b["qed"] - 0.05) \
                   and (f["sa"] <= b["sa"] + 0.3) and (f["scaf"] >= b["scaf"] * 0.8)
    lines.append("DECISION: " + ("PROMOTED" if promoted else "NOT_PROMOTED"))
    open(out, "w", encoding="utf-8").write("\n".join(lines))
    print("\n".join(lines))
