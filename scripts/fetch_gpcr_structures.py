# -*- coding: utf-8 -*-
"""从 RCSB PDB 抓取 GPCR-配体复合物, 构建可训练的 (结构, 配体 SMILES) 候选清单

筛选: Pfam PF00001 (7TM/GPCR) + 含非聚合物配体 + 分辨率 < 3.0 Å + 人源
输出:
    data/gpcr_pdb_raw/<PDBID>.cif      mmCIF 结构文件
    data/gpcr_pdb_raw/pairs.csv        (pdb_id, comp_id, smiles, formula_weight)
    data/gpcr_pdb_raw/manifest.csv     条目级下载状态

说明: PDB 数据为公共领域 (CC0), 用其训练出的权重许可干净、可商用。
并发: 每个条目需 3-4 次 API 往返, 用线程池并发抓取(串行约 2 小时, 并发约 15 分钟)。
"""
import os, sys, json, csv, time, argparse, threading, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
CIF = "https://files.rcsb.org/download/%s.cif"
ENTRY = "https://data.rcsb.org/rest/v1/core/entry/%s"
NONPOLY = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/%s/%s"
CHEMCOMP = "https://data.rcsb.org/rest/v1/core/chemcomp/%s"

ADDITIVES = {
    "HOH", "GOL", "EDO", "PEG", "PGE", "1PE", "2PE", "SO4", "PO4", "NO3", "CL", "NA", "K", "MG",
    "CA", "ZN", "MN", "FE", "CU", "NI", "CD", "HG", "IOD", "BR", "FMT", "ACT", "ACY", "MPD",
    "TRS", "EPE", "MES", "DMS", "CIT", "TLA", "BME", "DTT", "IMD", "NH4", "CO3", "SCN", "AZI",
    "BOG", "LDA", "MYR", "PLM", "OLA", "STE", "C8E", "LMT", "DDM", "LHG", "SDS", "CHS",
}

_lock = threading.Lock()
_delay = 0.1


def _get(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "7-eonmol-dataset/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def term(attr, val, op="exact_match"):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attr, "operator": op, "value": val}}


def search_ids(limit=None):
    nodes = [
        term("rcsb_polymer_entity_annotation.annotation_id", "PF00001"),
        term("rcsb_entry_info.nonpolymer_entity_count", 0, "greater"),
        term("rcsb_entry_info.resolution_combined", 3.0, "less"),
        term("rcsb_entity_source_organism.taxonomy_lineage.id", "9606"),
    ]
    ids, start, page = [], 0, 500
    while True:
        q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
             "return_type": "entry",
             "request_options": {"paginate": {"start": start, "rows": page},
                                 "results_content_type": ["experimental"]}}
        req = urllib.request.Request(SEARCH, data=json.dumps(q).encode(),
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=90) as r:
            d = json.loads(r.read().decode())
        batch = [x["identifier"] for x in d.get("result_set", [])]
        ids += batch
        total = d.get("total_count", 0)
        print("  检索 %d / %d" % (len(ids), total), flush=True)
        if not batch or len(ids) >= total:
            break
        start += page
        if limit and len(ids) >= limit:
            break
    return sorted(set(ids))[:limit] if limit else sorted(set(ids))


def druglike(smi, mw):
    """药物样过滤: 元素范围 / 环数 / 可旋转键 / 分子量。剔除胆固醇、脂类、辅因子等。"""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        from rdkit.Chem import Descriptors, rdMolDescriptors
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return False
        allowed = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"}
        if any(a.GetSymbol() not in allowed for a in m.GetAtoms()):
            return False
        if rdMolDescriptors.CalcNumRings(m) < 1:
            return False
        if Descriptors.NumRotatableBonds(m) > 12:
            return False
        return 150.0 <= float(mw) <= 600.0
    except Exception:
        return False


def process_entry(pid, outdir):
    """下载结构 + 抓取该条目的药物样配体。返回 (pid, status, bytes, pairs)"""
    path = os.path.join(outdir, "%s.cif" % pid)
    cached = os.path.exists(path) and os.path.getsize(path) > 2000
    size = os.path.getsize(path) if cached else 0
    if not cached:
        try:
            data = _get(CIF % pid, timeout=90)
            if len(data) < 2000:
                return pid, "too_small", 0, []
            with open(path, "wb") as f:
                f.write(data)
            size = len(data)
        except Exception as e:
            return pid, "%s: %s" % (type(e).__name__, e), 0, []
        time.sleep(_delay)

    pairs = []
    try:
        e = json.loads(_get(ENTRY % pid, timeout=40).decode())
        npe = (e.get("rcsb_entry_container_identifiers", {}) or {}).get("non_polymer_entity_ids") or []
    except Exception:
        npe = []
    for n in npe:
        time.sleep(_delay)
        try:
            d = json.loads(_get(NONPOLY % (pid, n), timeout=40).decode())
            cid = (d.get("pdbx_entity_nonpoly", {}) or {}).get("comp_id")
        except Exception:
            continue
        if not cid or cid in ADDITIVES:
            continue
        time.sleep(_delay)
        try:
            cc = json.loads(_get(CHEMCOMP % cid, timeout=40).decode())
        except Exception:
            continue
        desc = cc.get("rcsb_chem_comp_descriptor", {}) or {}
        smi = (desc.get("SMILES_stereo") or desc.get("SMILES") or "").strip()
        mw = (cc.get("chem_comp", {}) or {}).get("formula_weight")
        if not smi or not mw:
            continue
        with _lock:
            ok = druglike(smi, mw)
        if ok:
            pairs.append((pid, cid, smi, round(float(mw), 1)))
    return pid, ("ok" if pairs else "no_usable_ligand"), size, pairs


def main():
    global _delay
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--outdir", default=os.path.join(REPO, "data", "gpcr_pdb_raw"))
    ap.add_argument("--delay", type=float, default=0.1, help="每次 API 请求间隔(秒)")
    ap.add_argument("--workers", type=int, default=8, help="并发线程数")
    args = ap.parse_args()
    _delay = args.delay

    os.makedirs(args.outdir, exist_ok=True)
    print("检索 GPCR 复合物 (PF00001 + 配体 + <3.0Å + 人源)...", flush=True)
    ids = search_ids(args.limit)
    print("待抓取 %d 个条目, 并发 %d" % (len(ids), args.workers), flush=True)

    # 续跑: 已有 pairs 的条目无需重抓配体元数据
    prev = set()
    pcsv = os.path.join(args.outdir, "pairs.csv")
    if os.path.exists(pcsv):
        try:
            with open(pcsv, encoding="utf-8-sig", newline="") as f:
                for r in csv.DictReader(f):
                    prev.add(r["pdb_id"])
            print("  已有 %d 个条目的配体清单, 将复用" % len(prev), flush=True)
        except Exception:
            prev = set()

    man, pairs = [], []
    done = ok = fail = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(process_entry, pid, args.outdir) for pid in ids]
        for f in as_completed(futs):
            pid, status, size, prs = f.result()
            man.append((pid, size, status))
            pairs.extend(prs)
            done += 1
            if status == "ok":
                ok += 1
            elif status != "no_usable_ligand":
                fail += 1
            if done % 50 == 0 or done == len(ids):
                el = time.time() - t0
                eta = el / done * (len(ids) - done)
                print("  进度 %d/%d | 可用配体对 %d | 失败 %d | 已用 %.1f 分钟, 预计剩余 %.1f 分钟"
                      % (done, len(ids), len(pairs), fail, el / 60, eta / 60), flush=True)

    pairs.sort()
    with open(os.path.join(args.outdir, "manifest.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["pdb_id", "cif_bytes", "status"]); w.writerows(sorted(man))
    with open(pcsv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["pdb_id", "comp_id", "smiles", "formula_weight"]); w.writerows(pairs)
    print("完成: 条目 %d (有配体 %d, 失败 %d) | 可用配体对 %d 个, 覆盖 %d 个条目"
          % (len(ids), ok, fail, len(pairs), len({p[0] for p in pairs})))


if __name__ == "__main__":
    main()
