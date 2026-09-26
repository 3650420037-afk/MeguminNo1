# -*- coding: utf-8 -*-
"""从抓取的 GPCR-配体复合物构建训练数据集 (与现有 gpcr_multitarget_v2 格式一致)

对每个 (pdb_id, comp_id) 对:
    1. 解析 mmCIF 一次, 同时得到: 各配体的实验质心 + 全蛋白 PDB
       (项目 PDBProtein 只解析 PDB 格式, 故需转换; 合并为单次解析以省时间)
    2. 以该配体实验质心为口袋中心 (真实结合位点, 而非人工指定)
    3. 按半径提取口袋残基 -> pocket PDB
    4. 配体取**沉积结构的实验坐标**(用该配体的 CCD SMILES 校正键级) -> ligand SDF
       实验坐标不可用时默认丢弃该样本(--on-mismatch skip), 避免"重建构象只平移质心"
       带来的配体穿模; 如需旧行为可显式 --on-mismatch fallback

并发: 按结构并行(ProcessPool), 否则 600+ 个 mmCIF 串行需 1 小时以上。

输出: <outdir>/{pocket_*.pdb, ligands/*.sdf, index.pkl, split_by_name.pt, build_report.txt}
用法: python scripts/build_gpcr_v3.py [--workers 6]
"""
import os, sys, csv, pickle, random, argparse, warnings, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import DATA, src_on_path
src_on_path()

import numpy as np

warnings.simplefilter("ignore")

_GROUPS = {}
_ARGS = None


def parse_cif_once(cif_path, pdb_path):
    """解析 mmCIF 一次: 返回 ({comp_id: 质心}, {comp_id: [原子]}) 并写出全蛋白 PDB

    同时取出配体的**实验原子坐标** —— 这是关键: 早期实现只把 SMILES 重建的构象
    平移到实验质心, 扭转角与实验构象不同, 配体会插进蛋白壁(实测 16% 样本穿模,
    最近重原子距离低至 0.25-0.62 A, 而实验位姿为 2.5-3.3 A)。改用实验坐标可以
    从根上消除穿模。
    """
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    from Bio.PDB import MMCIFParser, PDBIO
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.simplefilter("ignore", PDBConstructionWarning)

    cents, lig_atoms = {}, {}
    s = MMCIFParser(QUIET=True).get_structure("x", cif_path)
    for res in s.get_residues():
        if not str(res.id[0]).strip().startswith("H_"):
            continue
        name = res.get_resname().strip()
        if name in cents:
            continue
        atoms = []
        for a in res.get_atoms():
            el = (getattr(a, "element", "") or "").strip()
            if not el or el.upper() == "H":
                continue
            atoms.append((a.get_name(), el, tuple(float(v) for v in a.coord)))
        if len(atoms) >= 5:
            cents[name] = np.array([x[2] for x in atoms], dtype=np.float64).mean(axis=0)
            lig_atoms[name] = atoms
    if not os.path.exists(pdb_path) or os.path.getsize(pdb_path) < 5000:
        io = PDBIO()
        io.set_structure(s)
        io.save(pdb_path)
    return cents, lig_atoms


def heavy_atom_count(smiles):
    """SMILES 的重原子数(用于诊断报告), 解析失败返回 -1。"""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles) if smiles else None
    return m.GetNumAtoms() if m is not None else -1


def build_from_experimental(atoms, smiles):
    """用实验原子坐标构造配体分子 (并用 SMILES 模板校正键级)。

    返回 (mol, status); status 取值:
        ok                    成功 —— mol 的坐标就是实验坐标(不穿模)
        bad_template          SMILES 无法解析
        pdb_parse_fail        PDB 块无法解析
        composition_mismatch  沉积原子组成与 SMILES 不一致
        template_fail         模板键级映射失败
        sanitize_fail         价态非法

    失败时 mol 为 None, 由调用方按 --on-mismatch 处理。

    为什么 composition_mismatch 必须单独区分: 它在真实数据里并不罕见 —— 例如
    5UEN 的 DU1, 沉积坐标漏建了磺酰氟上的 F(沉积 35 个重原子, CCD 定义 36 个)。
    这类条目下 "使用实验坐标" 与 "保持 SMILES 指定的完整分子" 不可兼得:
      - 回退到生成构象 -> 配体重新插进蛋白壁(正是本次要修的穿模)
      - 丢弃该样本     -> 少一条数据, 但留下的一律不穿模
    故默认丢弃, 并把这个选择权交给 --on-mismatch。
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None:
        return None, "bad_template"
    lines = []
    for i, (nm, el, (x, y, z)) in enumerate(atoms, 1):
        # 严格按 PDB 固定列位输出(列错位会让 RDKit 解析失败或读错元素)
        lines.append(
            "HETATM%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
            % (i, str(nm)[:4], " ", "LIG", "A", 1, " ", x, y, z, 1.00, 0.00,
               str(el).rjust(2)))
    block = "\n".join(lines) + "\nEND\n"
    mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True)
    if mol is None:
        return None, "pdb_parse_fail"
    # 元素组成必须与模板一致, 否则说明沉积结构里该配体不完整(缺原子/多原子),
    # 无法同时满足 "用实验坐标" 与 "SMILES 指定的分子"
    if sorted(a.GetAtomicNum() for a in mol.GetAtoms()) != \
       sorted(a.GetAtomicNum() for a in tmpl.GetAtoms()):
        return None, "composition_mismatch"
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(tmpl, mol)
    except Exception:
        return None, "template_fail"
    if mol is None or mol.GetNumConformers() == 0:
        return None, "template_fail"
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None, "sanitize_fail"
    return mol, "ok"


def empty_stats():
    return {"experimental": 0, "fallback": 0, "mismatch": 0,
            "unusable": {}, "mismatch_detail": []}


# 常见非类药配体(脂类/去垢剂/辅料/辅酶/核苷酸): 它们在共晶结构里很常见,
# 但不是药物分子, 用来训练会教模型生成胆固醇一类的油腻分子。
JUNK_CCD = {
    "CLR", "OLA", "OLB", "OLC", "PLM", "MYR", "STE", "LFA", "PC1", "Y01",
    "SOG", "BOG", "LMT", "C8E", "PEG", "PGE", "PG4", "1PE", "MYS",
    "FMN", "FAD", "NAD", "NAP", "GTP", "GDP", "GNP", "ATP", "ADP", "AMP", "ANP",
    "SO4", "PO4", "GOL", "EDO", "ACT", "DMS", "IMD", "TRS", "MPD", "CIT", "EPE",
    "NAG", "BMA", "MAN", "FUC", "GAL", "GLC", "NDG", "SIA", "XYP", "RIB",
}

# 模型元素词表(与 src/scripts/build_gpcr_dataset.py 一致): C N O F P S Cl。
# 含词表外元素的配体必须在建库阶段剔除, 否则训练时 transforms 会
# assert (cls_real.sum(-1) > 0).all() 失败 —— 实测曾在 iter 43 崩掉整次训练。
SUPPORTED_ELEMENTS = {6, 7, 8, 9, 15, 16, 17}

# 类药性过滤档位: (MW 下限, MW 上限, LogP 上限, QED 下限, 最少环数)
DRUGLIKE_TIERS = {
    "off": None,
    "loose": (200.0, 600.0, 6.5, 0.25, 2),
    "medium": (250.0, 600.0, 6.0, 0.30, 2),
    "strict": (250.0, 500.0, 5.0, 0.40, 1),
    # compact: 收紧分子量窗口以对齐目标输出尺寸。实测 medium 档训练集的 MW 中位
    # 387.9 但 P90 达 489.5(44% >=400), 微调后模型把尺寸先验推到了 MW 中位 496.6
    # (官方权重生成的是 401.5) —— 大分子更耗步数、QED 更低, 同时拉低产量与 QED。
    "compact": (250.0, 420.0, 5.0, 0.45, 2),
}


def curate_entries(entries, outdir, tier, max_per_ligand, skipped):
    """按类药性与冗余度筛选样本; 被淘汰的样本连同其文件一并移除。

    动机(实测): 738 个 RCSB GPCR 复合物里, 共晶配体有 31.4% 是胆固醇(CLR),
    另有去垢剂/辅酶/核苷酸等; CLR(253)+TEP(77) 两个分子就占了 41% 的样本。
    直接用这种数据微调会把类药分子生成器推向油腻/单一化学空间。

    返回 (kept_entries, stats)。stats 记录各淘汰原因计数。
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import QED, Descriptors, Crippen, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    stats = {"druglike_reject": 0, "junk_ccd": 0, "duplicate_capped": 0,
             "unreadable": 0, "bad_element": 0,
             "unique_ligands_before": 0, "unique_ligands_after": 0}
    spec = DRUGLIKE_TIERS.get(tier)

    # 1) 逐样本计算配体描述符与规范化 SMILES
    rows = []
    for e in entries:
        lpath = os.path.join(outdir, e[1])
        comp = os.path.basename(e[0])[len("pocket_"):-len(".pdb")].rpartition("_")[2]
        m = Chem.MolFromMolFile(lpath, sanitize=False, removeHs=True)
        if m is None:
            stats["unreadable"] += 1
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            stats["unreadable"] += 1
            continue
        # 元素词表检查(先于一切描述符; 词表外元素会让训练直接崩)
        bad = sorted({a.GetAtomicNum() for a in m.GetAtoms()} - SUPPORTED_ELEMENTS)
        if bad:
            stats["bad_element"] += 1
            skipped.append((os.path.basename(e[1]),
                            "词表外元素 %s (仅支持 C/N/O/F/P/S/Cl)" % bad))
            continue
        canon = Chem.MolToSmiles(m)
        if spec and comp in JUNK_CCD:
            stats["junk_ccd"] += 1
            skipped.append((os.path.basename(e[1]), "非类药配体(%s: 脂类/去垢剂/辅酶)" % comp))
            continue
        if spec:
            lo, hi, lpmax, qmin, rmin = spec
            mw = Descriptors.MolWt(m)
            lp = Crippen.MolLogP(m)
            q = QED.qed(m)
            nr = rdMolDescriptors.CalcNumRings(m)
            if not (lo <= mw <= hi):
                stats["druglike_reject"] += 1
                skipped.append((os.path.basename(e[1]),
                                "类药性: MW %.0f 不在 [%.0f,%.0f]" % (mw, lo, hi)))
                continue
            if lp > lpmax:
                stats["druglike_reject"] += 1
                skipped.append((os.path.basename(e[1]), "类药性: LogP %.2f > %.1f" % (lp, lpmax)))
                continue
            if q < qmin:
                stats["druglike_reject"] += 1
                skipped.append((os.path.basename(e[1]), "类药性: QED %.3f < %.2f" % (q, qmin)))
                continue
            if nr < rmin:
                stats["druglike_reject"] += 1
                skipped.append((os.path.basename(e[1]), "类药性: 环数 %d < %d" % (nr, rmin)))
                continue
        rows.append({"entry": e, "canon": canon, "comp": comp,
                     "pid": os.path.basename(e[0])[len("pocket_"):-len(".pdb")].rpartition("_")[0]})

    stats["unique_ligands_before"] = len({r["canon"] for r in rows})

    # 2) 同配体去冗余: 按 (规范化 SMILES) 分组, 每组最多保留 max_per_ligand 个
    #    按 pid 排序保证确定性(与随机种子无关, 结果可复现)
    by_lig = {}
    for r in sorted(rows, key=lambda x: (x["canon"], x["pid"], x["entry"][0])):
        by_lig.setdefault(r["canon"], []).append(r)
    kept = []
    for canon, group in by_lig.items():
        keep_n = max_per_ligand if max_per_ligand and max_per_ligand > 0 else len(group)
        kept.extend(group[:keep_n])
        for r in group[keep_n:]:
            stats["duplicate_capped"] += 1
            skipped.append((os.path.basename(r["entry"][1]),
                            "同配体冗余: %s 已达上限 %d" % (canon[:40], max_per_ligand)))
    stats["unique_ligands_after"] = len({r["canon"] for r in kept})

    # 3) 淘汰样本的文件删除, 避免留下无人引用的垃圾文件
    kept_ids = {r["entry"][0] for r in kept}
    for e in entries:
        if e[0] not in kept_ids:
            for rel in (e[0], e[1]):
                try:
                    os.remove(os.path.join(outdir, rel))
                except OSError:
                    pass
    return [r["entry"] for r in kept], stats


def process_pdb(args):
    """处理一个 PDB 条目下的所有配体对。返回 (entries, skipped, stats)"""
    (pid, items, outdir, work, radius, seed, base, on_mismatch) = args
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    entries, skipped = [], []
    stats = empty_stats()
    cif = os.path.join(base, "%s.cif" % pid)
    if not os.path.exists(cif):
        return [], [(pid, "cif missing")], stats
    pdb = os.path.join(work, "%s.pdb" % pid)
    try:
        cents, lig_atoms = parse_cif_once(cif, pdb)
    except Exception as e:
        return [], [(pid, "cif parse: %s" % e)], stats
    try:
        from utils.protein_ligand import PDBProtein
        protein = PDBProtein(pdb)
    except Exception as e:
        return [], [(pid, "PDBProtein: %s" % e)], stats
    from scripts.build_gpcr_dataset import build_conformer

    for k, (comp_id, smiles) in enumerate(items):
        center = cents.get(comp_id)
        if center is None:
            skipped.append(("%s/%s" % (pid, comp_id), "no experimental centroid"))
            continue
        try:
            residues = protein.query_residues_radius(center, radius, criterion="center_of_mass")
        except Exception as e:
            skipped.append(("%s/%s" % (pid, comp_id), "pocket: %s" % e))
            continue
        if not residues or len(residues) < 8:
            skipped.append(("%s/%s" % (pid, comp_id), "too few pocket residues"))
            continue
        pname = "pocket_%s_%s.pdb" % (pid, comp_id)
        try:
            with open(os.path.join(outdir, pname), "w", encoding="ascii") as h:
                h.write(protein.residues_to_pdb_block(residues, name="GPCR_POCKET"))
        except Exception as e:
            skipped.append(("%s/%s" % (pid, comp_id), "pocket write: %s" % e))
            continue
        # 优先使用实验坐标(不穿模)。失败的处置见 build_from_experimental 的说明:
        # 默认丢弃该样本, 而不是回退到已知会穿模的"生成构象 + 质心平移"。
        mol, status = None, "no_experimental_atoms"
        atoms = lig_atoms.get(comp_id)
        if atoms:
            try:
                mol, status = build_from_experimental(atoms, smiles)
            except Exception as e:
                mol, status = None, "exception: %s" % e
        if mol is not None:
            stats["experimental"] += 1
        else:
            stats["unusable"][status] = stats["unusable"].get(status, 0) + 1
            if status == "composition_mismatch":
                stats["mismatch"] += 1
                stats["mismatch_detail"].append(
                    (pid, comp_id, len(atoms or []), heavy_atom_count(smiles)))
            if on_mismatch == "fallback":
                stats["fallback"] += 1
                try:
                    mol = build_conformer(smiles, center, seed + base_index(pid, k))
                except Exception as e:
                    skipped.append(("%s/%s" % (pid, comp_id), "conformer: %s" % e))
                    try:
                        os.remove(os.path.join(outdir, pname))
                    except Exception:
                        pass
                    continue
            else:
                skipped.append(("%s/%s" % (pid, comp_id),
                                "ligand pose unusable: %s (已弃用, 不回退以免穿模)"
                                % status))
                try:
                    os.remove(os.path.join(outdir, pname))
                except Exception:
                    pass
                continue
        lname = os.path.join("ligands", "lig_%s_%s.sdf" % (pid, comp_id))
        try:
            w = Chem.SDWriter(os.path.join(outdir, lname))
            w.write(mol)
            w.close()
        except Exception as e:
            skipped.append(("%s/%s" % (pid, comp_id), "sdf write: %s" % e))
            continue
        entries.append((pname, lname, None, "0.0"))
    return entries, skipped, stats


def base_index(pid, k):
    """条目的稳定随机偏移。

    不能用内置 hash(): 它对 str 带进程随机盐(PYTHONHASHSEED), 会让回退构象的
    随机种子每次运行都不同, 同一份输入建出的数据集随之变化, 破坏可复现性。
    """
    import zlib
    return (zlib.crc32(str(pid).encode("utf-8")) % 100000) + k


def init_worker(groups, args_dict):
    global _GROUPS, _ARGS
    _GROUPS, _ARGS = groups, args_dict


def worker(pid):
    a = _ARGS
    return process_pdb((pid, _GROUPS[pid], a["out"], a["work"], a["radius"],
                        a["seed"], a["src"], a["on_mismatch"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(DATA, "gpcr_pdb_raw"))
    ap.add_argument("--out", default=os.path.join(DATA, "gpcr_v3"))
    ap.add_argument("--pocket-radius", type=float, default=12.0)
    ap.add_argument("--val-ratio", type=float, default=0.15)
    ap.add_argument("--split-by", choices=["structure", "pair"], default="structure",
                    help="划分单位: structure=按结构(PDB ID)划分(默认, 防同一口袋泄漏到验证集); "
                         "pair=按配体对随机划分(旧行为, 存在同口袋跨集泄漏)")
    ap.add_argument("--drug-like", choices=["off", "loose", "medium", "strict", "compact"],
                    default="medium",
                    help="共晶配体的类药性过滤档位(默认 medium): 落在 RCSB GPCR 复合物里的"
                         "配体有很大比例是胆固醇/去垢剂/辅酶, 不筛掉会教模型生成油腻分子。"
                         "compact 进一步收紧分子量窗口以对齐目标输出尺寸")
    ap.add_argument("--max-per-ligand", type=int, default=8,
                    help="同一配体(按规范化 SMILES)最多保留多少个不同口袋的样本, "
                         "0 表示不限制。实测 CLR 253 个、TEP 77 个, 两个分子占 41%% 样本")
    ap.add_argument("--seed", type=int, default=2021)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--on-mismatch", choices=["skip", "fallback"], default="skip",
                    help="配体无法用实验坐标构建时: skip=丢弃该样本(默认, 不引入穿模); "
                         "fallback=回退到生成构象+质心平移(与修复前行为一致, 会穿模)")
    ap.add_argument("--min-samples", type=int, default=10,
                    help="有效样本数低于此值即报错退出(便于小规模冒烟测试)")
    args = ap.parse_args()

    pairs_csv = os.path.join(args.src, "pairs.csv")
    if not os.path.exists(pairs_csv):
        raise SystemExit("缺少 %s" % pairs_csv)
    groups = {}
    with open(pairs_csv, encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            groups.setdefault(row["pdb_id"], []).append((row["comp_id"], row["smiles"]))
    n_pairs = sum(len(v) for v in groups.values())
    print("输入: %d 个结构, %d 个候选配体对; 并发 %d" % (len(groups), n_pairs, args.workers), flush=True)

    os.makedirs(os.path.join(args.out, "ligands"), exist_ok=True)
    work = os.path.join(args.src, "_pdb_cache")
    os.makedirs(work, exist_ok=True)

    adict = {"out": args.out, "work": work, "radius": args.pocket_radius,
             "seed": args.seed, "src": args.src, "on_mismatch": args.on_mismatch}
    entries, skipped, done = [], [], 0
    total = empty_stats()
    mismatch_rows = []
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx,
                             initializer=init_worker, initargs=(groups, adict)) as ex:
        futs = [ex.submit(worker, pid) for pid in sorted(groups)]
        for f in as_completed(futs):
            try:
                e, s, st = f.result()
            except Exception as e2:
                skipped.append(("worker", "%s" % e2))
                continue
            entries.extend(e)
            skipped.extend(s)
            total["experimental"] += st["experimental"]
            total["fallback"] += st["fallback"]
            total["mismatch"] += st["mismatch"]
            for k2, v2 in st["unusable"].items():
                total["unusable"][k2] = total["unusable"].get(k2, 0) + v2
            mismatch_rows.extend(st["mismatch_detail"])
            done += 1
            if done % 25 == 0 or done == len(groups):
                print("  结构 %d/%d | 已构建样本 %d | 跳过 %d"
                      % (done, len(groups), len(entries), len(skipped)), flush=True)

    if len(entries) < args.min_samples:
        raise SystemExit("有效样本过少 (%d < %d)" % (len(entries), args.min_samples))

    # ---- 数据治理: 类药性过滤 + 同配体去冗余 ----
    # 这一步在建库阶段就做, 而不是留给训练脚本, 以保证交付的数据集本身就是干净的
    entries, cur = curate_entries(entries, args.out, args.drug_like,
                                  args.max_per_ligand, skipped)
    print("数据治理(%s, 同配体上限 %s): 保留 %d | 非类药配体 -%d | 描述符不合格 -%d | "
          "词表外元素 -%d | 冗余裁剪 -%d | 不可解析 -%d"
          % (args.drug_like, args.max_per_ligand, len(entries), cur["junk_ccd"],
             cur["druglike_reject"], cur["bad_element"], cur["duplicate_capped"],
             cur["unreadable"]), flush=True)
    print("  唯一配体: 治理前 %d -> 治理后 %d" % (cur["unique_ligands_before"],
                                                  cur["unique_ligands_after"]), flush=True)
    if len(entries) < args.min_samples:
        raise SystemExit("类药性治理后有效样本过少 (%d < %d)" % (len(entries), args.min_samples))

    # ---- 划分 ----
    # 默认按**结构(PDB ID)**划分, 而不是按配体对随机划分。
    # 理由: 一个结构常含多个配体对(实测 738 结构 / 1009 对), 按对随机会把**同一个
    # 口袋**同时放进训练集和验证集, 验证损失因"见过这个口袋"而虚低 —— 这是数据泄漏。
    # 按结构划分保证验证集的口袋在训练集中从未出现, 验证指标才反映真实泛化。
    pid_of = {}
    for e in entries:
        # 口袋文件名形如 pocket_<pid>_<comp>.pdb
        name = os.path.basename(e[0])
        parts = name[len("pocket_"):-len(".pdb")].split("_") if name.startswith("pocket_") else []
        pid_of[e[0]] = parts[0] if parts else name

    random.seed(args.seed)
    if args.split_by == "structure":
        pids = sorted(set(pid_of.values()))
        random.shuffle(pids)
        n_val_pid = max(1, int(round(len(pids) * args.val_ratio)))
        val_pids = set(pids[:n_val_pid])
        train_entries = [e for e in entries if pid_of[e[0]] not in val_pids]
        val_entries = [e for e in entries if pid_of[e[0]] in val_pids]
    else:
        shuffled = entries[:]
        random.shuffle(shuffled)
        n_val = max(1, int(len(shuffled) * args.val_ratio))
        val_entries, train_entries = shuffled[:n_val], shuffled[n_val:]

    split = {"train": [tuple(e[:2]) for e in train_entries],
             "test": [tuple(e[:2]) for e in val_entries]}
    # 泄漏自检: 训练集与验证集的结构集合必须不相交
    tr_pids = {pid_of[e[0]] for e in train_entries}
    va_pids = {pid_of[e[0]] for e in val_entries}
    leaked = sorted(tr_pids & va_pids)
    if args.split_by == "structure" and leaked:
        raise SystemExit("划分自检失败: 训练/验证结构重叠 %d 个: %s" % (len(leaked), leaked[:5]))

    with open(os.path.join(args.out, "index.pkl"), "wb") as f:
        pickle.dump(entries, f)
    import torch
    torch.save(split, os.path.join(args.out, "split_by_name.pt"))
    with open(os.path.join(args.out, "build_report.txt"), "w", encoding="utf-8") as f:
        f.write("sources=rcsb_pfam_PF00001_human_ligand_bound\n")
        f.write("total=%d\ntrain=%d\nval=%d\n" % (len(entries), len(split["train"]), len(split["test"])))
        f.write("split_by=%s (seed=%d, val_ratio=%.3f)\n"
                % (args.split_by, args.seed, args.val_ratio))
        f.write("train_structures=%d\nval_structures=%d\nstructure_overlap=%d\n"
                % (len(tr_pids), len(va_pids), len(leaked)))
        f.write("expected_pairs_per_structure=%.2f\n"
                % (len(entries) / max(1, len({pid_of[e[0]] for e in entries}))))
        f.write("on_mismatch=%s\n" % args.on_mismatch)
        f.write("drug_like=%s\n" % args.drug_like)
        f.write("max_per_ligand=%s\n" % args.max_per_ligand)
        f.write("curation: junk_ccd=%d druglike_reject=%d bad_element=%d duplicate_capped=%d unreadable=%d\n"
                % (cur["junk_ccd"], cur["druglike_reject"], cur["bad_element"],
                   cur["duplicate_capped"], cur["unreadable"]))
        f.write("unique_ligands: before=%d after=%d\n"
                % (cur["unique_ligands_before"], cur["unique_ligands_after"]))
        f.write("ligand_source=experimental_coords(default)_or_generated_fallback\n")
        f.write("experimental_pose=%d\n" % total["experimental"])
        f.write("fallback_generated=%d\n" % total["fallback"])
        f.write("composition_mismatch=%d\n" % total["mismatch"])
        f.write("unusable_reasons=%s\n" % (sorted(total["unusable"].items()) or "none"))
        f.write("structures=%d\npocket_radius=%.1f\nskipped=%d\n"
                % (len({e[0] for e in entries}), args.pocket_radius, len(skipped)))
        for s, why in skipped[:60]:
            f.write("  skip %s: %s\n" % (s, why))
    # 沉积结构不完整的配体明细(供人工复核: 差了几个重原子)
    csv_path = os.path.join(args.out, "unusable_ligands.csv")
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pdb_id", "comp_id", "deposited_heavy_atoms", "smiles_heavy_atoms", "delta"])
        for pid, cid, nd, ns in sorted(mismatch_rows):
            w.writerow([pid, cid, nd, ns, ns - nd if ns >= 0 else ""])
    print("完成: %d 对样本 (train %d / val %d), 覆盖 %d 个结构; 跳过 %d"
          % (len(entries), len(split["train"]), len(split["test"]),
             len({e[0] for e in entries}), len(skipped)))
    print("配体位姿: 实验坐标 %d | 回退生成 %d | 组成不匹配 %d"
          % (total["experimental"], total["fallback"], total["mismatch"]))
    if total["unusable"]:
        print("不可用原因: %s" % sorted(total["unusable"].items()))
    print("输出: %s" % args.out)


if __name__ == "__main__":
    main()
