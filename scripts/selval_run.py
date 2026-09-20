# -*- coding: utf-8 -*-
"""任务7 灵敏度验证驱动: 四组对接 (A2A选择性/A1选择性 × A2A受体/A1受体)"""
import os, subprocess, sys, time

REPO = r"D:\MMModel\Pocket2Mol"
PY = r"D:\Miniconda3\envs\Pocket2Mol\python.exe"
KNOWN = r"D:\MMModel\已知药物库"
os.chdir(REPO)

A2A_PDB = r"D:\MMModel\靶点结构\4EIY_A2A受体.pdb"
A1_PDB = r"D:\MMModel\靶点结构\5UEN.pdb"
JOBS = [
    ("sel_a2a", A2A_PDB, "-0.4,8.5,17.1", os.path.join(KNOWN, "selectivity_a2a_selective.txt")),
    ("sel_a2a", A1_PDB, "55.969,58.897,143.624", os.path.join(KNOWN, "selectivity_a2a_selective.txt")),
    ("sel_a1", A2A_PDB, "-0.4,8.5,17.1", os.path.join(KNOWN, "selectivity_a1_selective.txt")),
    ("sel_a1", A1_PDB, "55.969,58.897,143.624", os.path.join(KNOWN, "selectivity_a1_selective.txt")),
]

prog = os.path.join(REPO, "outputs", "selval_progress.txt")
RCODE = {A2A_PDB: "4EIY", A1_PDB: "5UEN"}  # 输出目录只用 ASCII 受体代码 (RDKit 中文路径坑)
for g, rec, center, smif in JOBS:
    rcode = RCODE[rec]
    outdir = os.path.join(REPO, "outputs", "selval_%s_%s" % (g, rcode))
    logf = open(os.path.join(os.environ["TEMP"], "selval_%s_%s.log" % (g, rcode)), "w", encoding="utf-8")
    if os.path.exists(os.path.join(outdir, "summary.csv")):
        print("skip (已有):", outdir)
        continue
    r = subprocess.run([PY, "scripts/docking_pipeline.py", "--receptor", rec,
                        "--center=" + center, "--smiles-file", smif,
                        "--out", outdir], cwd=REPO, stdout=logf, stderr=subprocess.STDOUT)
    logf.close()
    line = "%s %s %s EXIT=%d" % (time.strftime("%H:%M"), g, rcode, r.returncode)
    with open(prog, "a", encoding="utf-8") as f:
        f.write(line + "\n")
    print(line)
print("SELVAL_DONE")
