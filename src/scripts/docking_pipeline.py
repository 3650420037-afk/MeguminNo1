#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pocket2Mol docking pipeline: AutoDock Vina + Open Babel.

Workflow:
  1. Receptor PDB  -> PDBQT   (obabel -ipdb -opdbqt -xr; waters/HETATM stripped first)
  2. Ligand SMILES -> RDKit add-H + 3D embed + MMFF/UFF optimize -> SDF
                             -> obabel -> PDBQT
     (or: --sdf-dir uses existing SDF files directly)
  3. Vina config   (center, size 22.5 A default, exhaustiveness 8)
  4. vina --receptor --ligand --config --out ; parse lowest binding energy
  5. Summary CSV   (smiles, vina_score, status)

IMPORTANT (this machine): subprocess stdout/stderr are redirected to FILE
handles, never pipes (named pipes are forbidden by the local sandbox -> EPERM).

Usage (run with the Pocket2Mol conda python; paths relative to the repo root):
  python src/scripts/docking_pipeline.py ^
      --receptor data/targets/rec.pdb ^
      --center=-0.4,8.5,17.1 ^
      --smiles-file smiles.txt ^
      --out outputs/docking_out ^
      [--size 22.5] [--exhaustiveness 8] [--sdf-dir DIR] [--smiles "CCO" ...]
      [--seed 42] [--vina-timeout 1800]

  --center accepts three numbers ("-0.4 8.5 17.1") or one comma string; when
  the first number is negative use the "--center=x,y,z" form so the shell does
  not treat it as an option.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import OBABEL, OBABEL_DATA, VINA, src_on_path
src_on_path()

# ---------------------------------------------------------------------------
# Tool locations (来自 paths.py, 可用环境变量覆盖)
# ---------------------------------------------------------------------------
VINA_EXE = VINA
OBABEL_EXE = OBABEL
BABEL_DATADIR = OBABEL_DATA

POCKET_SIZE = 22.5          # Angstrom, x=y=z
EXHAUSTIVENESS = 8
VINA_TIMEOUT_S = 1800       # per-ligand wall clock limit

# affinity rows: "   1       -8.853          0          0"
# (rmsd columns may be plain integers for the best mode, e.g. "0")
SCORE_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$")
SANITIZE_RE = re.compile(r"[^A-Za-z0-9_\-]+")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def tool_env():
    env = os.environ.copy()
    env["BABEL_DATADIR"] = BABEL_DATADIR
    return env


def run_cmd(cmd, log_path, err_path, timeout=None):
    """Run a command; stdout/stderr -> files. NO capture_output/pipes."""
    with open(log_path, "w", encoding="utf-8", errors="replace") as fo, \
         open(err_path, "w", encoding="utf-8", errors="replace") as fe:
        try:
            proc = subprocess.run(cmd, stdout=fo, stderr=fe, timeout=timeout)
            return proc.returncode, False
        except subprocess.TimeoutExpired:
            return -999, True


def read_text(path):
    # utf-8-sig: 兼容带 BOM 的文件。旧写法 utf-8 + line.strip() 无法去除 U+FEFF
    # (它不是空白字符), 会让首行 smiles 带着 BOM 进入 summary.csv, 与库 CSV 主键
    # 不一致而静默丢样 (实测每靶点丢 1 个候选)。
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        return f.read()


def parse_center(tokens):
    parts = []
    for tok in tokens:
        parts.extend(t for t in re.split(r"[,\s]+", tok.strip()) if t)
    if len(parts) != 3:
        raise ValueError("--center needs exactly 3 numbers, got: %r" % tokens)
    return tuple(float(p) for p in parts)


# ---------------------------------------------------------------------------
# receptor preparation
# ---------------------------------------------------------------------------
def clean_receptor_pdb(src_pdb, dst_pdb):
    """Keep ATOM/TER records only (drop HETATM waters/ligands/ions, altloc != A)."""
    kept = dropped = 0
    with open(src_pdb, "r", encoding="utf-8", errors="replace") as fi, \
         open(dst_pdb, "w", encoding="utf-8", newline="") as fo:
        for line in fi:
            rec = line[:6].strip()
            if rec == "ATOM":
                altloc = line[16:17]
                if altloc and altloc not in (" ", "A"):
                    continue
                fo.write(line if line.endswith("\n") else line + "\n")
                kept += 1
            elif rec == "TER":
                fo.write("TER\n")
            elif rec == "HETATM":
                dropped += 1
    return kept, dropped


def pdbqt_stats(path):
    """Return (n_atoms, n_H, sum_charge) of a PDBQT file."""
    n = n_h = 0
    qsum = 0.0
    for line in read_text(path).splitlines():
        if line.startswith(("ATOM", "HETATM")):
            n += 1
            elem = line[77:79].strip().upper()
            name = line[12:16].strip().upper()
            if elem.startswith("H") or (not elem and name.startswith("H")):
                n_h += 1
            try:
                qsum += float(line[70:76])
            except ValueError:
                pass
    return n, n_h, qsum


def prepare_receptor(pdb_path, workdir):
    """PDB -> PDBQT. Strip waters/HETATM first; warn (not fail) on poor H/charge."""
    rec_pdbqt = os.path.join(workdir, "receptor.pdbqt")
    logdir = os.path.join(workdir, "logs")
    notes = []
    txt = read_text(pdb_path)
    n_hetatm = sum(1 for l in txt.splitlines() if l.startswith("HETATM"))
    n_hoh = sum(1 for l in txt.splitlines()
                if l.startswith("HETATM") and l[17:20].strip() == "HOH")
    if n_hoh or n_hetatm:
        clean_pdb = os.path.join(workdir, "receptor_clean.pdb")
        kept, dropped = clean_receptor_pdb(pdb_path, clean_pdb)
        notes.append("stripped %d HETATM (%d water) from receptor; %d ATOM kept"
                     % (dropped, n_hoh, kept))
        src = clean_pdb
    else:
        src = pdb_path

    cmd = [OBABEL_EXE, "-ipdb", src, "-opdbqt", "-xr", "-O", rec_pdbqt]
    rc, _ = run_cmd(cmd,
                    os.path.join(logdir, "receptor_obabel.log"),
                    os.path.join(logdir, "receptor_obabel.err"))
    if rc != 0 or not os.path.exists(rec_pdbqt) or \
            os.path.getsize(rec_pdbqt) == 0:
        notes.append("first obabel attempt failed (rc=%s); retry with -p 7.4"
                     % rc)
        cmd = [OBABEL_EXE, "-ipdb", src, "-opdbqt", "-xr", "-p", "7.4",
               "-O", rec_pdbqt]
        rc, _ = run_cmd(cmd,
                        os.path.join(logdir, "receptor_obabel_retry.log"),
                        os.path.join(logdir, "receptor_obabel_retry.err"))
        if rc != 0 or not os.path.exists(rec_pdbqt) or \
                os.path.getsize(rec_pdbqt) == 0:
            return None, notes + ["receptor conversion FAILED"]

    n_atoms, n_h, qsum = pdbqt_stats(rec_pdbqt)
    if n_h == 0:
        notes.append("WARNING: receptor PDBQT has no hydrogen atoms "
                     "(obabel H/charge quality uncertain); continuing")
    if abs(qsum) > 5.0:
        notes.append("WARNING: receptor net charge %.2f far from 0; "
                     "check charges" % qsum)
    notes.append("receptor PDBQT: %d atoms (%d H), net charge %.2f"
                 % (n_atoms, n_h, qsum))
    return rec_pdbqt, notes


# ---------------------------------------------------------------------------
# ligand preparation
# ---------------------------------------------------------------------------
def prep_ligand_smiles(smiles, name, ligdir, logdir):
    """RDKit: add H, 3D embed, MMFF/FF optimize -> SDF -> obabel -> PDBQT."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "invalid_smiles"
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0x5EED) == -1:
        if AllChem.EmbedMolecule(mol, randomSeed=0x5EED,
                                 useRandomCoords=True) == -1:
            return None, "embed_failed"
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        ff = "MMFF"
    elif AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        ff = "UFF"
    else:
        ff = "none"

    mol.SetProp("_Name", name)
    sdf = os.path.join(ligdir, name + ".sdf")
    writer = Chem.SDWriter(sdf)
    writer.write(mol)
    writer.close()

    return convert_sdf_to_pdbqt(sdf, name, ligdir, logdir, ff)


def convert_sdf_to_pdbqt(sdf, name, ligdir, logdir, ff_tag="asis"):
    pdbqt = os.path.join(ligdir, name + ".pdbqt")
    rc, _ = run_cmd([OBABEL_EXE, "-isdf", sdf, "-opdbqt", "-O", pdbqt],
                    os.path.join(logdir, name + "_obabel.log"),
                    os.path.join(logdir, name + "_obabel.err"))
    if rc != 0 or not os.path.exists(pdbqt) or os.path.getsize(pdbqt) == 0:
        return None, "obabel_failed"
    return pdbqt, "ok(%s)" % ff_tag


# ---------------------------------------------------------------------------
# vina
# ---------------------------------------------------------------------------
def write_vina_config(path, rec_pdbqt, center, size, exhaustiveness, seed):
    lines = [
        "receptor = %s" % rec_pdbqt,
        "center_x = %.3f" % center[0],
        "center_y = %.3f" % center[1],
        "center_z = %.3f" % center[2],
        "size_x = %.1f" % size,
        "size_y = %.1f" % size,
        "size_z = %.1f" % size,
        "exhaustiveness = %d" % exhaustiveness,
        "num_modes = 9",
    ]
    if seed is not None:
        lines.append("seed = %d" % seed)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def parse_vina_score(log_path):
    """Lowest affinity (kcal/mol) from vina stdout table; None if absent."""
    best = None
    in_table = False
    for line in read_text(log_path).splitlines():
        if line.startswith("-----+"):
            in_table = True
            continue
        if in_table:
            m = SCORE_ROW_RE.match(line)
            if m:
                aff = float(m.group(2))
                best = aff if best is None else min(best, aff)
    return best


def run_vina(rec_pdbqt, lig_pdbqt, cfg_path, out_pdbqt, name, logdir,
             timeout):
    cmd = [VINA_EXE, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--config", cfg_path, "--out", out_pdbqt]
    vlog = os.path.join(logdir, name + "_vina.log")
    verr = os.path.join(logdir, name + "_vina.err")
    rc, timed_out = run_cmd(cmd, vlog, verr, timeout=timeout)
    if timed_out:
        return None, "timeout"
    if rc != 0:
        return None, "vina_failed(rc=%d)" % rc
    score = parse_vina_score(vlog)
    if score is None:
        return None, "vina_no_score"
    return score, "ok"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="AutoDock Vina docking pipeline (Vina + Open Babel)")
    ap.add_argument("--receptor", required=True, help="receptor PDB path")
    ap.add_argument("--center", nargs="+", required=True,
                    help="pocket center: x y z | x,y,z (use --center=x,y,z "
                         "when x is negative)")
    ap.add_argument("--smiles", action="append", default=[],
                    help="one ligand SMILES (repeatable)")
    ap.add_argument("--smiles-file",
                    help="text file, one SMILES per line")
    ap.add_argument("--sdf-dir", help="directory of prepared SDF ligands")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--size", type=float, default=POCKET_SIZE,
                    help="box size x=y=z in Angstrom (default 22.5)")
    ap.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--vina-timeout", type=int, default=VINA_TIMEOUT_S)
    args = ap.parse_args(argv)

    center = parse_center(args.center)

    # ---- collect ligand inputs ------------------------------------------
    ligands = []  # (name, smiles_or_None, sdf_or_None)
    if args.smiles_file:
        for i, line in enumerate(read_text(args.smiles_file).splitlines(), 1):
            s = line.strip().lstrip("\ufeff")      # 双保险: 去 BOM
            if s:
                ligands.append(("lig_%03d" % len(ligands), s, None))
    for s in args.smiles:
        s = s.strip()
        if s:
            ligands.append(("lig_%03d" % len(ligands), s, None))
    if args.sdf_dir:
        for fn in sorted(os.listdir(args.sdf_dir)):
            if fn.lower().endswith(".sdf"):
                stem = SANITIZE_RE.sub("_",
                                       os.path.splitext(fn)[0]) or "sdf"
                ligands.append(("sdf_" + stem, None,
                                os.path.join(args.sdf_dir, fn)))
    if not ligands:
        ap.error("no ligands: give --smiles / --smiles-file / --sdf-dir")

    # ---- workdirs --------------------------------------------------------
    outdir = os.path.abspath(args.out)
    for sub in ("", "receptor", "ligands", "docking", "logs"):
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)
    logdir = os.path.join(outdir, "logs")
    ligdir = os.path.join(outdir, "ligands")

    log("receptor: %s" % args.receptor)
    log("center:   (%.3f, %.3f, %.3f)  size: %.1f  exhaustiveness: %d"
        % (center + (args.size, args.exhaustiveness)))
    log("ligands:  %d" % len(ligands))

    # ---- receptor --------------------------------------------------------
    rec_pdbqt, notes = prepare_receptor(args.receptor, outdir)
    for n in notes:
        log("receptor: " + n)
    if rec_pdbqt is None:
        log("FATAL: receptor preparation failed, aborting")
        return 2

    # ---- vina config -----------------------------------------------------
    cfg_path = os.path.join(outdir, "vina_config.txt")
    write_vina_config(cfg_path, rec_pdbqt, center, args.size,
                      args.exhaustiveness, args.seed)
    log("vina config: %s" % cfg_path)

    # ---- per-ligand ------------------------------------------------------
    from rdkit import Chem  # noqa: F401  (fail fast before heavy work)
    from rdkit.Chem import AllChem  # noqa: F401

    rows = []
    t0 = time.time()
    for idx, (name, smiles, sdf) in enumerate(ligands, 1):
        status = ""
        score = None
        if sdf is not None:
            pdbqt, status = convert_sdf_to_pdbqt(sdf, name, ligdir, logdir)
            # 用 open()+MolFromMolBlock 读 SDF: RDKit C++ 文件 API 在 Windows 上打不开
            # 含中文的路径 (旧写法 MolFromMolFile 抛 OSError, 被 except 吞掉后 smiles="",
            # 使 summary.csv 无法与库 CSV 关联)。
            # 同时 MolToSmiles 前必须 RemoveHs: 否则输出 [H]OC([H])([H])... 显式氢 SMILES,
            # 与库中的规范 SMILES 永久不匹配。
            try:
                with open(sdf, encoding="utf-8", errors="ignore") as fh:
                    block = fh.read()
                m = Chem.MolFromMolBlock(block, removeHs=False, sanitize=True)
                if m is None:
                    smiles = ""
                    log("  WARN: %s SDF 解析失败(MolFromMolBlock=None), smiles 留空" % name)
                else:
                    smiles = Chem.MolToSmiles(Chem.RemoveHs(m))
            except Exception as e:
                smiles = ""
                log("  WARN: %s 读取 SDF 异常 %s: %s" % (name, type(e).__name__, e))
        else:
            pdbqt, status = prep_ligand_smiles(smiles, name, ligdir, logdir)
        if pdbqt is None:
            log("[%d/%d] %s  prep FAILED: %s" % (idx, len(ligands), name,
                                                 status))
            rows.append((name, smiles, "", status))
            continue

        out_pdbqt = os.path.join(outdir, "docking", name + "_out.pdbqt")
        score, status = run_vina(rec_pdbqt, pdbqt, cfg_path, out_pdbqt,
                                 name, logdir, args.vina_timeout)
        log("[%d/%d] %s  score=%s  (%s)"
            % (idx, len(ligands), name,
               ("%.2f" % score) if score is not None else "n/a", status))
        rows.append((name, smiles,
                     "" if score is None else "%.2f" % score, status))

    # ---- summary CSV ------------------------------------------------------
    csv_path = os.path.join(outdir, "summary.csv")
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "smiles", "vina_score", "status"])
        w.writerows(rows)

    ok = [r for r in rows if r[3] == "ok" and r[2] != ""]
    log("done in %.1fs -> %s" % (time.time() - t0, csv_path))
    log("summary: %d ok / %d total" % (len(ok), len(rows)))
    for name, smiles, sc, st in rows:
        log("  %-10s %-10s %-8s %s" % (name, sc or "-", st, smiles[:60]))
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
