# -*- coding: utf-8 -*-
"""循环清洗: 审计源数据集孤立原子 -> 剔除 -> 重建 -> 合并 -> 审计合并产物, 直到零脏."""
import os, pickle, random, shutil, subprocess, sys
sys.path.insert(0, r"D:\MMModel\Pocket2Mol")
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
import torch

BASE = r"D:\MMModel\Pocket2Mol\data"
PY = r"D:\Miniconda3\envs\Pocket2Mol\python.exe"
SRC = ["gpcr_a2a", "gpcr_b2ar", "gpcr_d3", "gpcr_5ht2b"]
TOTAL_REMOVED = {}

def audit_dir(d, index):
    bad = []
    for pocket, lig, _, _ in index:
        m = next(iter(Chem.SDMolSupplier(os.path.join(d, lig), removeHs=False, sanitize=False)), None)
        if m is None:
            bad.append(lig); continue
        bonded = set()
        for b in m.GetBonds():
            bonded.add(b.GetBeginAtomIdx()); bonded.add(b.GetEndAtomIdx())
        if any(i not in bonded for i in range(m.GetNumAtoms())):
            bad.append(lig)
    return bad

def rebuild_split(d, entries):
    random.seed(2021); es = entries[:]; random.shuffle(es)
    sp = max(1, int(len(es) * 0.8))
    torch.save({"train": [(e[0], e[1]) for e in es[:sp]], "test": [(e[0], e[1]) for e in es[sp:]]},
               os.path.join(d, "split_by_name.pt"))

# 1) 循环清洗源数据集
for ds in SRC:
    d = os.path.join(BASE, ds)
    for round_i in range(10):
        index = pickle.load(open(os.path.join(d, "index.pkl"), "rb"))
        bad = audit_dir(d, index)
        if not bad:
            print("%s round%d: 干净 (%d 对)" % (ds, round_i, len(index)))
            break
        clean = [e for e in index if e[1] not in bad]
        print("%s round%d: %d -> %d (剔 %s)" % (ds, round_i, len(index), len(clean), bad))
        pickle.dump(clean, open(os.path.join(d, "index.pkl"), "wb"))
        rebuild_split(d, clean)
        for suf in ["_processed.lmdb", "_name2id.pt"]:
            p = os.path.join(BASE, ds + suf)
            if os.path.exists(p): os.remove(p)
        TOTAL_REMOVED[ds] = TOTAL_REMOVED.get(ds, 0) + len(bad)

# 2) 删合并产物旧 lmdb 后重新合并
for suf in ["gpcr_multitarget_v2_processed.lmdb", "gpcr_multitarget_v2_name2id.pt"]:
    p = os.path.join(BASE, suf)
    if os.path.exists(p): os.remove(p)
r = subprocess.run([PY, os.path.join(r"D:\MMModel\Pocket2Mol\scripts", "merge_gpcr_datasets.py"),
                    "--inputs"] + [os.path.join(BASE, s) for s in SRC] +
                   ["--output", os.path.join(BASE, "gpcr_multitarget_v2")],
                   stdout=open(os.path.join(BASE, "merge_log.txt"), "w"), stderr=subprocess.STDOUT)
print("merge exit:", r.returncode)

# 3) 审计合并产物, 若有脏则从 index 剔并重建 split (SDF 文件冗余无碍)
d = os.path.join(BASE, "gpcr_multitarget_v2")
for round_i in range(5):
    index = pickle.load(open(os.path.join(d, "index.pkl"), "rb"))
    bad = audit_dir(d, index)
    if not bad:
        print("multitarget_v2 round%d: 干净 (%d 对)" % (round_i, len(index)))
        break
    clean = [e for e in index if e[1] not in bad]
    print("multitarget_v2 round%d: %d -> %d (剔 %s)" % (round_i, len(index), len(clean), bad))
    pickle.dump(clean, open(os.path.join(d, "index.pkl"), "wb"))
    rebuild_split(d, clean)
    for suf in ["_processed.lmdb", "_name2id.pt"]:
        p = os.path.join(d[:-0] if False else BASE, "gpcr_multitarget_v2" + suf)
        if os.path.exists(p): os.remove(p)
print("TOTAL_REMOVED:", TOTAL_REMOVED)
