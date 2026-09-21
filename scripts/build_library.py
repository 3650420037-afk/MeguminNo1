# -*- coding: utf-8 -*-
"""build_library.py -- Pocket2Mol 产出过滤入库脚本.

读取 outputs/<run_name>/<session_dir>/SMILES.txt, 全库 RDKit 规范化去重,
按 PAINS/BRENK、MW、LogP、QED、SA、元素、环数过滤后入库:
  <library>/compounds.csv   列: smiles, 来源runs, qed, sa, mw, logp, tpsa,
                            hbd, hba, rings, murcko_scaffold, pdbqt_ready
  <library>/library.db      SQLite 表 compounds (同结构) + 索引 (+rejected/meta 表)
  <library>/sdf/            每个分子一个 3D SDF (ETKDG + MMFF94s/UFF 优化),
                            文件名 mol_<骨架hash8>_<分子hash12>.sdf
  <library>/rejected.csv    被剔除分子及原因
  <library>/summary.txt     总数/过滤漏斗/骨架数/QED 直方图等

用法 (python 一律用 conda 全路径):
  D:\\Miniconda3\\envs\\Pocket2Mol\\python.exe D:\\MMModel\\Pocket2Mol\\scripts\\build_library.py ^
      --runs D:\\MMModel\\Pocket2Mol\\outputs\\run_a D:\\MMModel\\Pocket2Mol\\outputs\\run_b ^
      [--library D:\\MMModel\\化合物库] [--no-filter] [--append] [--workers 4]

  --runs      一个或多个 run 目录 (空格分隔, 顺序无关)。也可传 outputs 下的
              批次容器目录 (如 outputs\\overnight_a2a), 其下每个含 SMILES.txt
              的子目录会被当作独立 run 计入来源。
  --library   输出库目录, 默认 D:\\MMModel\\化合物库
  --no-filter 跳过过滤, 所有规范化去重后的分子直接入库 (仍计算描述符)
  --append    增量入库: 合并已有 library.db 的旧记录 (来源取并集), 否则全量重建
  --workers   SDF 3D 构象生成 multiprocessing 进程数, 默认 4

说明:
  * pdbqt_ready=True 表示该分子的 3D 构象 SDF 已成功生成 (sdf/ 下有对应文件),
    可直接交给 OpenBabel/AutoDockTools 类工具转 pdbqt。
  * SA score 通过仓库自带 evaluation/sascorer 计算: 先 sys.path.insert(0, 仓库根),
    再 import utils.guidance 借其 rdkit.six shim, 最后
    from evaluation.sascorer import calculateScore。
  * 全程不使用 subprocess 管道捕获; 所有文件 utf-8, 文件名 ASCII。
"""

import argparse
import csv
import glob as _glob
import hashlib
import os
import sqlite3
import sys
import time
from collections import Counter

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_LIBRARY = r"D:\MMModel\化合物库"

ALLOWED_ELEMENTS = ("C", "N", "O", "F", "P", "S", "Cl")
MW_MIN, MW_MAX = 250.0, 500.0
LOGP_MIN, LOGP_MAX = 1.0, 5.0
QED_MIN = 0.4
SA_MAX = 6.0
MIN_RINGS = 1
ETKDG_SEED = 0xF00D

CSV_HEADER = ["smiles", "来源runs", "qed", "sa", "mw", "logp", "tpsa",
              "hbd", "hba", "rings", "murcko_scaffold", "pdbqt_ready"]
REJECTED_HEADER = ["smiles", "来源runs", "reason"]

_SA_CALC = None


# --------------------------------------------------------------------------
# SA score (借 utils.guidance 的 rdkit.six shim)
# --------------------------------------------------------------------------
def _load_sa_score():
    global _SA_CALC
    if _SA_CALC is None:
        if REPO_ROOT not in sys.path:
            sys.path.insert(0, REPO_ROOT)
        import utils.guidance  # noqa: F401  仅借其 rdkit.six shim
        from evaluation.sascorer import calculateScore
        _SA_CALC = calculateScore
    return _SA_CALC


# --------------------------------------------------------------------------
# 输入收集
# --------------------------------------------------------------------------
def resolve_run_dir(raw):
    p = os.path.abspath(raw)
    if os.path.isdir(p):
        return p
    alt = os.path.join(REPO_ROOT, "outputs", raw)
    if os.path.isdir(alt):
        return alt
    raise SystemExit("[ERROR] run 目录不存在: %s" % raw)


def collect_runs(run_dirs):
    """返回 [(run_label, [SMILES.txt 路径...]), ...]。

    支持三种输入形态:
      run 目录本身含 <session>/SMILES.txt  -> label = 目录名
      run 目录直接含 SMILES.txt            -> label = 目录名
      批次容器目录 (子目录才是 run)        -> 每个含 SMILES.txt 的子目录一个 label
    """
    out = []
    for raw in run_dirs:
        p = resolve_run_dir(raw)
        sessions = sorted(_glob.glob(os.path.join(p, "*", "SMILES.txt")))
        if sessions:
            out.append((os.path.basename(p), sessions))
            continue
        direct = sorted(_glob.glob(os.path.join(p, "SMILES.txt")))
        if direct:
            out.append((os.path.basename(p), direct))
            continue
        sub = []
        for child in sorted(os.listdir(p)):
            cp = os.path.join(p, child)
            if not os.path.isdir(cp):
                continue
            cf = (sorted(_glob.glob(os.path.join(cp, "*", "SMILES.txt")))
                  or sorted(_glob.glob(os.path.join(cp, "SMILES.txt"))))
            if cf:
                sub.append((child, cf))
        if not sub:
            raise SystemExit("[ERROR] 该目录下没有 SMILES.txt: %s" % p)
        out.extend(sub)
    return out


# --------------------------------------------------------------------------
# 规范化 / 描述符
# --------------------------------------------------------------------------
def canonical_smiles(text):
    from rdkit import Chem
    m = Chem.MolFromSmiles(text)
    if m is None:
        return None
    return Chem.MolToSmiles(m, isomericSmiles=True)


def compute_props(mol):
    """计算入库描述符; 单项失败置 None, 不抛异常。"""
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold

    p = {"qed": None, "sa": None, "mw": None, "logp": None, "tpsa": None,
         "hbd": None, "hba": None, "rings": None, "murcko_scaffold": ""}
    try:
        p["mw"] = Descriptors.MolWt(mol)
        p["logp"] = Crippen.MolLogP(mol)
        p["tpsa"] = rdMolDescriptors.CalcTPSA(mol)
        p["hbd"] = Lipinski.NumHDonors(mol)
        p["hba"] = Lipinski.NumHAcceptors(mol)
        p["rings"] = rdMolDescriptors.CalcNumRings(mol)
    except Exception:
        pass
    try:
        p["qed"] = QED.qed(mol)
    except Exception:
        pass
    try:
        p["sa"] = _load_sa_score()(mol)
    except Exception:
        pass
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        p["murcko_scaffold"] = Chem.MolToSmiles(scaf) if scaf.GetNumAtoms() > 0 else ""
    except Exception:
        pass
    return p


def build_catalog(catalog_id):
    from rdkit.Chem import FilterCatalog
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(catalog_id)
    return FilterCatalog.FilterCatalog(params)


def check_elements(mol):
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in ALLOWED_ELEMENTS:
            return sym
    return None


# --------------------------------------------------------------------------
# 3D 构象 (multiprocessing worker, 每个 worker 进程 rdBase 独立)
# --------------------------------------------------------------------------
def _sdf_init():
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")


def _build_3d(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, "parse_fail"
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = ETKDG_SEED
    ps.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(m, ps) == -1:
        ps.useRandomCoords = True
        ps.maxIterations = 5000
        if AllChem.EmbedMolecule(m, ps) == -1:
            return None, "embed_fail"
    try:
        if AllChem.MMFFHasAllMoleculeParams(m):
            AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=2000)
        else:
            AllChem.UFFOptimizeMolecule(m, maxIters=2000)
    except Exception:
        pass  # 力场优化失败不影响 ETKDG 3D 坐标有效性
    return m, None


def _sdf_worker(task):
    smiles, out_path = task
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    try:
        m, err = _build_3d(smiles)
        if m is None:
            return (smiles, False, err)
        m.SetProp("_Name", "mol_" + hashlib.md5(smiles.encode("utf-8")).hexdigest()[:12])
        block = Chem.MolToMolBlock(m, kekulize=True)
        text = block + ">  <smiles>\n%s\n\n$$$$\n" % smiles
        with open(out_path, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        return (smiles, True, "")
    except Exception as exc:  # 任何异常都不允许拖垮整个池
        return (smiles, False, "exception: %s" % exc)


def generate_sdfs(rows, workers):
    """为 rows 中尚未有 3D SDF 的分子生成构象文件, 返回 {smiles: (ok, msg)}。"""
    import multiprocessing as mp

    tasks = []
    for smi, row in rows.items():
        out_path = row["_sdf_path"]
        if out_path and os.path.isfile(out_path):  # 已存在 (append/重复构建) 直接复用
            row["pdbqt_ready"] = True
            continue
        tasks.append((smi, out_path))
    if not tasks:
        return {}
    results = {}
    workers = max(1, min(int(workers), len(tasks)))
    with mp.Pool(processes=workers, initializer=_sdf_init) as pool:
        for i, (smi, ok, msg) in enumerate(
                pool.imap_unordered(_sdf_worker, tasks, chunksize=8), 1):
            results[smi] = (ok, msg)
            if i % 250 == 0 or i == len(tasks):
                print("  [sdf] %d / %d done" % (i, len(tasks)), flush=True)
    return results


# --------------------------------------------------------------------------
# 库写出
# --------------------------------------------------------------------------
def sdf_names(canon, scaffold):
    """SDF 文件名: 骨架 hash 前缀(按骨架聚类) + 分子 hash(唯一且幂等)。"""
    mol_h = hashlib.md5(canon.encode("utf-8")).hexdigest()[:12]
    scaf_h = hashlib.md5(scaffold.encode("utf-8")).hexdigest()[:8] if scaffold else "chain"
    return "mol_%s_%s.sdf" % (scaf_h, mol_h)


def _fmt(v, nd):
    return "" if v is None else round(float(v), nd)


def write_outputs(lib_dir, rows, rejected, funnel_lines, stats_lines, meta):
    csv_path = os.path.join(lib_dir, "compounds.csv")
    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_HEADER)
        for smi in sorted(rows):
            r = rows[smi]
            w.writerow([smi, ";".join(sorted(r["sources"])),
                        _fmt(r["qed"], 4), _fmt(r["sa"], 3), _fmt(r["mw"], 2),
                        _fmt(r["logp"], 3), _fmt(r["tpsa"], 2),
                        r["hbd"] if r["hbd"] is not None else "",
                        r["hba"] if r["hba"] is not None else "",
                        r["rings"] if r["rings"] is not None else "",
                        r["murcko_scaffold"],
                        "True" if r["pdbqt_ready"] else "False"])

    rej_path = os.path.join(lib_dir, "rejected.csv")
    with open(rej_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(REJECTED_HEADER)
        for text, labels, reason in rejected:
            w.writerow([text, ";".join(sorted(labels)), reason])

    db_path = os.path.join(lib_dir, "library.db")
    conn = sqlite3.connect(db_path)
    try:
        conn.execute(
            'CREATE TABLE IF NOT EXISTS compounds ('
            ' smiles TEXT NOT NULL PRIMARY KEY,'
            ' "来源runs" TEXT NOT NULL,'
            ' qed REAL, sa REAL, mw REAL, logp REAL, tpsa REAL,'
            ' hbd INTEGER, hba INTEGER, rings INTEGER,'
            ' murcko_scaffold TEXT,'
            ' pdbqt_ready INTEGER NOT NULL DEFAULT 0)')
        conn.execute("DROP TABLE IF EXISTS rejected")
        conn.execute(
            'CREATE TABLE rejected ('
            ' smiles TEXT NOT NULL, "来源runs" TEXT NOT NULL, reason TEXT NOT NULL)')
        conn.execute("DROP TABLE IF EXISTS meta")
        conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
        conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_compounds_smiles"
                     " ON compounds(smiles)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_compounds_scaffold"
                     " ON compounds(murcko_scaffold)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_compounds_qed"
                     " ON compounds(qed)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_compounds_sa"
                     " ON compounds(sa)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_compounds_pdbqt"
                     " ON compounds(pdbqt_ready)")
        conn.executemany(
            'INSERT OR REPLACE INTO compounds (smiles, "来源runs", qed, sa, mw,'
            ' logp, tpsa, hbd, hba, rings, murcko_scaffold, pdbqt_ready)'
            ' VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
            [(smi, ";".join(sorted(rows[smi]["sources"])),
              rows[smi]["qed"], rows[smi]["sa"], rows[smi]["mw"],
              rows[smi]["logp"], rows[smi]["tpsa"], rows[smi]["hbd"],
              rows[smi]["hba"], rows[smi]["rings"], rows[smi]["murcko_scaffold"],
              1 if rows[smi]["pdbqt_ready"] else 0) for smi in sorted(rows)])
        conn.executemany(
            'INSERT INTO rejected (smiles, "来源runs", reason) VALUES (?,?,?)',
            [(t, ";".join(sorted(l)), r) for t, l, r in rejected])
        conn.executemany(
            'INSERT OR REPLACE INTO meta (key, value) VALUES (?,?)', meta)
        conn.commit()
    finally:
        conn.close()

    sum_path = os.path.join(lib_dir, "summary.txt")
    with open(sum_path, "w", encoding="utf-8", newline="\n") as fh:
        for k, v in meta:
            fh.write("%-24s: %s\n" % (k, v))
        fh.write("\n" + "\n".join(funnel_lines) + "\n\n")
        fh.write("\n".join(stats_lines) + "\n")
    return csv_path, db_path, rej_path, sum_path


def load_old_rows(db_path):
    """append 模式: 读取已有 compounds 表 -> {smiles: row dict}。"""
    if not os.path.isfile(db_path):
        return {}
    try:
        conn = sqlite3.connect(db_path)
        cur = conn.execute(
            'SELECT smiles, "来源runs", qed, sa, mw, logp, tpsa, hbd, hba,'
            ' rings, murcko_scaffold, pdbqt_ready FROM compounds')
        rows = {}
        for (smi, src, qed, sa, mw, logp, tpsa, hbd, hba, rings, scaf, ready) in cur:
            rows[smi] = {"sources": set(x for x in (src or "").split(";") if x),
                         "qed": qed, "sa": sa, "mw": mw, "logp": logp,
                         "tpsa": tpsa, "hbd": hbd, "hba": hba, "rings": rings,
                         "murcko_scaffold": scaf or "", "pdbqt_ready": bool(ready),
                         "_sdf_path": None}
        conn.close()
        return rows
    except Exception as exc:
        print("[WARN] 读取旧库失败, 按全新构建处理: %s" % exc)
        return {}


# --------------------------------------------------------------------------
# 汇总文本
# --------------------------------------------------------------------------
def qed_histogram(qeds):
    counts = [0] * 10
    n = 0
    for q in qeds:
        if q is None:
            continue
        counts[min(int(q * 10), 9)] += 1
        n += 1
    if n == 0:
        return ["---- QED distribution (无可统计分子) ----"]
    mx = max(counts) if counts else 1
    if mx <= 0:                      # 防御: 全 None 时 max(counts)=0 会导致下方除零
        mx = 1
    lines = ["---- QED distribution (n=%d) ----" % n]
    for i, c in enumerate(counts):
        bar = "#" * int(round(40.0 * c / mx))
        lines.append("  %.1f-%.1f | %-40s %d" % (i / 10.0, (i + 1) / 10.0, bar, c))
    return lines


def mini_stats(name, values):
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return "  %-6s: n/a" % name
    med = vals[len(vals) // 2]
    return ("  %-6s: min=%.3f  median=%.3f  max=%.3f  n=%d"
            % (name, vals[0], med, vals[-1], len(vals)))


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Pocket2Mol 产出过滤入库: SMILES.txt -> 化合物库")
    ap.add_argument("--runs", nargs="+", required=True, metavar="RUN_DIR",
                    help="outputs 下一个或多个 run 目录 (空格分隔)")
    ap.add_argument("--library", default=DEFAULT_LIBRARY, metavar="DIR",
                    help="输出库目录 (默认 %(default)s)")
    ap.add_argument("--no-filter", action="store_true",
                    help="跳过 PAINS/BRENK/MW/LogP/QED/SA/元素/环数 过滤")
    ap.add_argument("--append", action="store_true",
                    help="合并已有 library.db 旧记录 (来源取并集)")
    ap.add_argument("--workers", type=int, default=4, metavar="N",
                    help="SDF 生成进程数 (默认 4)")
    args = ap.parse_args(argv)

    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    import rdkit

    t0 = time.perf_counter()
    lib_dir = os.path.abspath(args.library)
    sdf_dir = os.path.join(lib_dir, "sdf")
    os.makedirs(sdf_dir, exist_ok=True)

    # ---- 1. 读取全部 SMILES.txt, 规范化, 全库去重 ----
    run_files = collect_runs(args.runs)
    print("== build_library ==")
    print("rdkit %s | python %s" % (rdkit.__version__, sys.version.split()[0]))
    print("runs:")
    for label, files in run_files:
        print("  %-40s %d SMILES.txt" % (label, len(files)))

    raw_lines, invalid = 0, 0
    per_run_raw = Counter()
    records = {}   # canon_smiles -> set(run labels)
    rejected = []  # (text, labels, reason)
    t_read = time.perf_counter()
    for label, files in run_files:
        for path in files:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    text = line.strip()
                    if not text:
                        continue
                    raw_lines += 1
                    per_run_raw[label] += 1
                    canon = canonical_smiles(text)
                    if canon is None:
                        invalid += 1
                        rejected.append((text, {label}, "invalid_smiles"))
                        continue
                    records.setdefault(canon, set()).add(label)
    t_read = time.perf_counter() - t_read
    deduped = len(records)

    # ---- 2. 过滤 (漏斗: 各环节用桶计数, 最后统一汇总) ----
    kill = Counter()  # 桶: reparse/pains/brenk/mw/logp/qed/sa/element/no_ring
    rows = {}         # canon -> row dict
    if not args.no_filter:
        from rdkit.Chem import FilterCatalog
        pains = build_catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        brenk = build_catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        allowed = set(ALLOWED_ELEMENTS)
        t_f = time.perf_counter()
        for smi in sorted(records):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:  # canonical round-trip 理论不会失败, 防御
                kill["reparse"] += 1
                rejected.append((smi, records[smi], "reparse_fail"))
                continue
            props = compute_props(mol)

            entry = pains.GetFirstMatch(mol)
            if entry is not None:
                kill["pains"] += 1
                rejected.append((smi, records[smi], "pains: %s" % entry.GetDescription()))
                continue
            entry = brenk.GetFirstMatch(mol)
            if entry is not None:
                kill["brenk"] += 1
                rejected.append((smi, records[smi], "brenk: %s" % entry.GetDescription()))
                continue
            if props["mw"] is None:
                kill["mw"] += 1
                rejected.append((smi, records[smi], "prop_error: mw"))
                continue
            if not (MW_MIN <= props["mw"] <= MW_MAX):
                kill["mw"] += 1
                rejected.append((smi, records[smi], "mw_out_of_range (%.1f)" % props["mw"]))
                continue
            if props["logp"] is None:
                kill["logp"] += 1
                rejected.append((smi, records[smi], "prop_error: logp"))
                continue
            if not (LOGP_MIN <= props["logp"] <= LOGP_MAX):
                kill["logp"] += 1
                rejected.append((smi, records[smi], "logp_out_of_range (%.2f)" % props["logp"]))
                continue
            if props["qed"] is None or props["qed"] < QED_MIN:
                kill["qed"] += 1
                rejected.append((smi, records[smi],
                                 "qed_below_%g (%s)" % (QED_MIN, props["qed"])))
                continue
            if props["sa"] is None or props["sa"] > SA_MAX:
                kill["sa"] += 1
                rejected.append((smi, records[smi],
                                 "sa_above_%g (%s)" % (SA_MAX, props["sa"])))
                continue
            bad = check_elements(mol)
            if bad is not None:
                kill["element"] += 1
                rejected.append((smi, records[smi], "element_not_allowed: %s" % bad))
                continue
            if props["rings"] is None or props["rings"] < MIN_RINGS:
                kill["no_ring"] += 1
                rejected.append((smi, records[smi], "rings_below_%d" % MIN_RINGS))
                continue
            rows[smi] = dict(props, sources=set(records[smi]),
                             pdbqt_ready=False, _sdf_path=None)
        t_filter = time.perf_counter() - t_f
    else:
        t_filter = 0.0
        for smi in sorted(records):
            mol = Chem.MolFromSmiles(smi)
            base = {"qed": None, "sa": None, "mw": None, "logp": None, "tpsa": None,
                    "hbd": None, "hba": None, "rings": None, "murcko_scaffold": ""}
            if mol is not None:
                base.update(compute_props(mol))
            rows[smi] = dict(base, sources=set(records[smi]),
                             pdbqt_ready=False, _sdf_path=None)

    # ---- 3. append 合并旧库 ----
    if args.append:
        old = load_old_rows(os.path.join(lib_dir, "library.db"))
        for smi, orow in old.items():
            if smi in rows:
                rows[smi]["sources"] |= orow["sources"]
            else:
                rows[smi] = orow

    # ---- 4. SDF 3D 构象 (Pool(4), 每 worker rdBase 独立) ----
    for smi, row in rows.items():
        row["_sdf_path"] = os.path.join(sdf_dir, sdf_names(smi, row["murcko_scaffold"]))
    t_sdf = time.perf_counter()
    if rows:
        print("generating 3D SDFs for %d molecules (%d workers)..."
              % (len(rows), args.workers), flush=True)
        sdf_results = generate_sdfs(rows, args.workers)
        for smi, (ok, msg) in sdf_results.items():
            rows[smi]["pdbqt_ready"] = bool(ok)
            if not ok:
                rejected.append((smi, rows[smi]["sources"], "sdf_fail: %s" % msg))
    t_sdf = time.perf_counter() - t_sdf

    # ---- 5. 漏斗/统计 + 写出 ----
    funnel_lines = ["---- filter funnel ----"]
    funnel_lines.append("  %-38s %6d" % ("raw SMILES lines read", raw_lines))
    funnel_lines.append("  %-38s %6d" % ("invalid/unparseable (dropped)", invalid))
    funnel_lines.append("  %-38s %6d" % ("unique canonical SMILES (deduped)", deduped))
    if args.no_filter:
        funnel_lines.append("  %-38s %6d" % ("filtering SKIPPED (--no-filter)", deduped))
    else:
        stages = [
            ("after PAINS", kill["pains"] + kill["reparse"]),
            ("after BRENK", kill["brenk"]),
            ("after MW in [%g, %g]" % (MW_MIN, MW_MAX), kill["mw"]),
            ("after Crippen LogP in [%g, %g]" % (LOGP_MIN, LOGP_MAX), kill["logp"]),
            ("after QED >= %g" % QED_MIN, kill["qed"]),
            ("after SA score <= %g" % SA_MAX, kill["sa"]),
            ("after elements %s" % "".join(ALLOWED_ELEMENTS), kill["element"]),
            ("after rings >= %d" % MIN_RINGS, kill["no_ring"]),
        ]
        cum = deduped
        for name, killed in stages:
            cum -= killed
            funnel_lines.append("  %-38s %6d  (-%d, %.1f%% of unique)"
                                % (name, cum, killed, 100.0 * cum / deduped if deduped else 0.0))

    n_ready = sum(1 for r in rows.values() if r["pdbqt_ready"])
    scaffolds = Counter(r["murcko_scaffold"] for r in rows.values() if r["murcko_scaffold"])

    funnel_lines.append("  %-38s %6d" % ("final library compounds", len(rows)))
    funnel_lines.append("  %-38s %6d" % ("with 3D SDF (pdbqt_ready)", n_ready))

    stats_lines = ["---- library stats ----"]
    stats_lines.append("unique murcko scaffolds: %d" % len(scaffolds))
    stats_lines.append("top scaffolds:")
    for scaf, c in scaffolds.most_common(10):
        stats_lines.append("  %-58s %4d" % (scaf[:58], c))
    stats_lines.append(mini_stats("qed", [r["qed"] for r in rows.values()]))
    stats_lines.append(mini_stats("sa", [r["sa"] for r in rows.values()]))
    stats_lines.append(mini_stats("mw", [r["mw"] for r in rows.values()]))
    stats_lines.append(mini_stats("logp", [r["logp"] for r in rows.values()]))
    stats_lines.extend(qed_histogram([r["qed"] for r in rows.values()]))
    stats_lines.append("---- per run (raw lines / final kept+merged / rejected) ----")
    for label, _files in run_files:
        contrib = sum(1 for r in rows.values() if label in r["sources"])
        rej = sum(1 for _t, l, _r in rejected if label in l)
        stats_lines.append("  %-40s raw=%-6d kept=%-5d rejected=%d"
                           % (label, per_run_raw[label], contrib, rej))
    total_s = time.perf_counter() - t0
    stats_lines.append("---- timing ----")
    stats_lines.append("  read+normalize %.1fs | filter+props %.1fs | sdf %.1fs | total %.1fs"
                       % (t_read, t_filter, t_sdf, total_s))

    meta = [("script", os.path.abspath(__file__)),
            ("rdkit", rdkit.__version__),
            ("runs", ";".join(label for label, _ in run_files)),
            ("mode", "append" if args.append else "rebuild"),
            ("filter", "off" if args.no_filter else
             "PAINS+BRENK, MW[%g,%g], LogP[%g,%g], QED>=%g, SA<=%g, elements=%s, rings>=%d"
             % (MW_MIN, MW_MAX, LOGP_MIN, LOGP_MAX, QED_MIN, SA_MAX,
                "".join(ALLOWED_ELEMENTS), MIN_RINGS)),
            ("raw_lines", str(raw_lines)),
            ("invalid_smiles", str(invalid)),
            ("unique_canonical", str(deduped)),
            ("kept", str(len(rows))),
            ("pdbqt_ready", str(n_ready)),
            ("elapsed_sec", "%.1f" % total_s)]

    csv_path, db_path, rej_path, sum_path = write_outputs(
        lib_dir, rows, rejected, funnel_lines, stats_lines, meta)

    print("")
    print("\n".join(funnel_lines))
    print("rejected (reason recorded): %d -> %s" % (len(rejected), rej_path))
    print("library: %d compounds (%d pdbqt_ready) | scaffolds: %d"
          % (len(rows), n_ready, len(scaffolds)))
    print("outputs:")
    print("  %s" % csv_path)
    print("  %s" % db_path)
    print("  %s" % sum_path)
    print("total %.1fs" % total_s)


if __name__ == "__main__":
    main()
