# -*- coding: utf-8 -*-
"""用 AutoDock Vina 把各靶点的 ChEMBL 活性分子对接进**实验口袋**, 构建
"每口袋高密度 + 物理可行位姿"的微调数据集。

为什么这样做
------------
微调破坏生长/终止回路的根因已定位为「每口袋约 1 个训练样本」(见
docs/微调实验_v3_报告.md §四)。旧数据位姿虽坏(31.5% 粗穿模), 但「4 个口袋 x ~96 样本」
让模型把每个口袋的"多大多满就停"记了下来, 产量因此保住了(完成 67 vs 基线 63)。
本脚本用对接一次性拿到两者:
  - 每口袋上百样本  -> 恢复"每口袋密度"
  - Vina 的位姿不含原子重叠 -> 物理可行(用与 diagnose_dataset 完全相同的判据复核)

产物(与 utils.datasets.pl.PocketLigandDataset 兼容)
--------------------------------------------------
    <out>/pocket_<T>.pdb              12 A 口袋残基(每靶点一个, 该靶点配体共享)
    <out>/ligands/lig_<T>_<name>.sdf  对接最优位姿
    <out>/index.pkl                   [(pocket, ligand, None, "0.0"), ...]
    <out>/split_by_name.pt            {'train': [...], 'test': [...]}
    <out>/build_report.txt            含治理计数与位姿质量统计

用法
----
    python src/scripts/build_docked_dataset.py \
        --out data/gpcr_dock_v1 --work outputs/dock_build \
        --exhaustiveness 8 --val-ratio 0.15

注意: 划分按**配体对分层随机**(同一靶点内随机), 验证集与训练集共享口袋。
      这是有意为之: 评测靶点(A2A)本身就在训练集中, 度量目标是对同一口袋生成新分子的
      能力(配体空间泛化), 而非跨靶点泛化。报告里会明确写出这一点。
"""
import argparse
import csv
import glob
import os
import pickle
import random
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, DATA, TARGETS, CONFIGS, LIGANDS, OBABEL, ensure_dir  # noqa: E402
from paths import src_on_path  # noqa: E402
src_on_path()

# 与 src/scripts/build_gpcr_v3.py / build_gpcr_dataset.py 保持一致的模型元素词表
SUPPORTED_ELEMENTS = {6, 7, 8, 9, 15, 16, 17}
# 类药性(与 build_gpcr_v3.py 的 medium 档一致)
MW_LO, MW_HI, LOGP_MAX, QED_MIN, RINGS_MIN = 250.0, 600.0, 6.0, 0.30, 2
# 位姿物理可行性: 与 diagnose_dataset.py 的 LIGAND_CLOSE_CONTACT 同界(0.5~1.5 A 记为紧接触)
MIND_MIN = 1.5
POCKET_RADIUS = 12.0

# 靶点 -> (活性配体 CSV, 共晶配体 cid)  靶点名取自 configs/targets.json
TARGET_CSV = {
    "A2A": "A2A_活性配体.csv",
    "B2AR": "B2AR_活性配体.csv",
    "D3": "D3_活性配体.csv",
    "5HT2B": "5HT2B_活性配体.csv",
    # ---- M6 批次（2026-09-28）：按受体扩口袋用的人源 Class A 靶点 ----
    # 这些靶点没有旧的 *_活性配体.csv，只有 fetch_chembl_actives.py 的 *_活性配体_full.csv
    "D2": "D2_活性配体_full.csv",
    "5HT1A": "5HT1A_活性配体_full.csv",
    "5HT2A": "5HT2A_活性配体_full.csv",
    "M1": "M1_活性配体_full.csv",
    "H1": "H1_活性配体_full.csv",
    "OPRM1": "OPRM1_活性配体_full.csv",
    "OPRK1": "OPRK1_活性配体_full.csv",
    "A1": "A1_活性配体_full.csv",
}


def log(msg):
    print(msg, flush=True)


def read_smiles(csv_path):
    out = []
    with open(csv_path, encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            s = (row.get("canonical_smiles") or "").strip()
            if s:
                out.append(s)
    return out


def make_pocket(pdb_path, center, out_pdb):
    """按 12 A 提取口袋残基(与 sample_for_pdb / diagnose_dataset 同口径)。"""
    from utils.protein_ligand import PDBProtein
    prot = PDBProtein(pdb_path)
    residues = prot.query_residues_radius(center, POCKET_RADIUS, criterion="center_of_mass")
    if len(residues) < 8:
        return 0
    with open(out_pdb, "w", encoding="ascii") as h:
        h.write(prot.residues_to_pdb_block(residues, name="GPCR_POCKET"))
    return len(residues)


def pdbqt_to_sdf(pdbqt, sdf, logdir):
    """用 obabel 把 vina 输出位姿转 SDF。不加 --gen3d, 保留对接坐标。"""
    cmd = [OBABEL, "-ipdbqt", pdbqt, "-osdf", "-O", sdf]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True)
    if r.returncode != 0 or not os.path.exists(sdf):
        with open(os.path.join(logdir, os.path.basename(pdbqt) + ".err"), "wb") as f:
            f.write(r.stderr or b"")
        return False
    return True


def validate_ligand(sdf, pocket_pos, mind_min=MIND_MIN):
    """返回 (ok, reason, mind)。检查元素词表/类药性/位姿是否与蛋白重叠。"""
    import numpy as np
    from rdkit import Chem, RDLogger
    from rdkit.Chem import QED, Descriptors, Crippen, rdMolDescriptors
    from scripts.diagnose_dataset import min_distances
    RDLogger.DisableLog("rdApp.*")

    m = Chem.MolFromMolFile(sdf, sanitize=False, removeHs=True)
    if m is None:
        return False, "unreadable", None
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return False, "sanitize_fail", None
    # 注意: 只校验重原子。obabel 从 PDBQT 转出的 SDF 会带显式氢, 而 sanitize=False 时
    # removeHs=True 并不保证去掉 —— 氢(原子序数 1)不在词表里, 会把全部样本误杀。
    bad = sorted({a.GetAtomicNum() for a in m.GetAtoms() if a.GetAtomicNum() != 1}
                 - SUPPORTED_ELEMENTS)
    if bad:
        return False, "element_%s" % bad, None
    mw = Descriptors.MolWt(m)
    lp = Crippen.MolLogP(m)
    q = QED.qed(m)
    nr = rdMolDescriptors.CalcNumRings(m)
    if not (MW_LO <= mw <= MW_HI):
        return False, "mw_%.0f" % mw, None
    if lp > LOGP_MAX:
        return False, "logp_%.2f" % lp, None
    if q < QED_MIN:
        return False, "qed_%.3f" % q, None
    if nr < RINGS_MIN:
        return False, "rings_%d" % nr, None
    conf = m.GetConformer()
    lp_pos = np.asarray([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                          conf.GetAtomPosition(i).z] for i in range(m.GetNumAtoms())],
                        dtype=float)
    _, mind = min_distances(lp_pos, pocket_pos)
    if mind is None:
        return False, "no_pocket", None
    if mind < mind_min:
        return False, "clash_%.3f" % mind, mind
    return True, "ok", mind


def main():
    ap = argparse.ArgumentParser(description="对接构建高密度微调数据集")
    ap.add_argument("--out", default=os.path.join(DATA, "gpcr_dock_v1"))
    ap.add_argument("--work", default=os.path.join(ROOT, "outputs", "dock_build"))
    ap.add_argument("--targets", nargs="+", default=list(TARGET_CSV))
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--parallel", type=int, default=6,
                    help="每个靶点把配体切成 N 块并发对接(docking_pipeline 单进程串行, "
                         "实测约 1.7 个/分钟, 并发可把总时间缩短数倍)")
    ap.add_argument("--val-ratio", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=2021)
    ap.add_argument("--mind-min", type=float, default=MIND_MIN,
                    help="配体-蛋白最近重原子距离下限(A), 低于此值判为穿模并丢弃")
    ap.add_argument("--reuse-dock", action="store_true",
                    help="复用 --work 下已完成的对接结果, 不重新对接")
    args = ap.parse_args()

    import json
    import numpy as np
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from scripts.diagnose_dataset import pocket_geometry

    reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    outdir = ensure_dir(args.out)
    ligdir = ensure_dir(os.path.join(outdir, "ligands"))
    ensure_dir(args.work)

    entries, skipped, per_target, minds_all = [], [], {}, []
    for T in args.targets:
        if T not in reg:
            log("  [跳过] %s 不在 configs/targets.json" % T)
            continue
        t = reg[T]
        pdb = os.path.join(ROOT, t["pdb"])
        center = list(t["center"])
        csv_path = os.path.join(LIGANDS, TARGET_CSV[T])
        smis = read_smiles(csv_path)
        log("=" * 78)
        log("[%s] 活性分子 %d 个 | 受体 %s | 中心 %s"
            % (T, len(smis), os.path.basename(pdb), ",".join("%.2f" % v for v in center)))

        # 1) 口袋
        pname = "pocket_%s.pdb" % T
        nres = make_pocket(pdb, center, os.path.join(outdir, pname))
        if nres == 0:
            log("  口袋残基不足, 跳过该靶点")
            continue

        # 2) 对接
        dockdir = os.path.join(args.work, "dock_" + T)
        smi_file = os.path.join(args.work, "smiles_%s.txt" % T)
        os.makedirs(args.work, exist_ok=True)
        with open(smi_file, "w", encoding="ascii") as f:
            f.write("\n".join(smis) + "\n")
        summary = os.path.join(dockdir, "summary.csv")
        if not (args.reuse_dock and os.path.exists(summary)):
            # 并行分块: docking_pipeline.py 单进程串行跑配体, 实测约 1.7 个/分钟
            # (B2AR 191 个需近 2 小时)。切成 N 块并发跑, 再合并汇总与位姿文件。
            npar = max(1, min(args.parallel, len(smis)))
            if npar == 1:
                cmd = [sys.executable,
                       os.path.join(ROOT, "src", "scripts", "docking_pipeline.py"),
                       "--receptor", pdb,
                       "--center=" + ",".join("%.3f" % v for v in center),
                       "--smiles-file", smi_file, "--out", dockdir,
                       "--exhaustiveness", str(args.exhaustiveness)]
                with open(os.path.join(args.work, "dock_%s.log" % T), "w",
                          encoding="utf-8", errors="ignore") as lf:
                    rc = subprocess.run(cmd, cwd=ROOT, stdout=lf,
                                        stderr=subprocess.STDOUT).returncode
                log("  对接返回码 %d" % rc)
            else:
                chunks = [[] for _ in range(npar)]
                for i, s in enumerate(smis):
                    chunks[i % npar].append(s)
                procs = []
                for ci, ch in enumerate(chunks):
                    if not ch:
                        continue
                    cf = os.path.join(args.work, "smiles_%s_p%d.txt" % (T, ci))
                    with open(cf, "w", encoding="ascii") as f:
                        f.write("\n".join(ch) + "\n")
                    cdir = os.path.join(args.work, "dock_%s_p%d" % (T, ci))
                    cmd = [sys.executable,
                           os.path.join(ROOT, "src", "scripts", "docking_pipeline.py"),
                           "--receptor", pdb,
                           "--center=" + ",".join("%.3f" % v for v in center),
                           "--smiles-file", cf, "--out", cdir,
                           "--exhaustiveness", str(args.exhaustiveness)]
                    lf = open(os.path.join(args.work, "dock_%s_p%d.log" % (T, ci)), "w",
                              encoding="utf-8", errors="ignore")
                    procs.append((subprocess.Popen(cmd, cwd=ROOT, stdout=lf,
                                                   stderr=subprocess.STDOUT), cdir, lf))
                log("  并行 %d 块对接中 ..." % len(procs))
                for p, cdir, lf in procs:
                    p.wait()
                    lf.close()
                # 合并。注意: docking_pipeline 在每个分块内部都把配体命名为
                # lig_000, lig_001 ... —— 若直接合并, 6 个分块会互相覆盖, 最终只剩
                # 一个分块的量(实测 A2A 175 个只剩 30 个)。故统一加分块前缀。
                ensure_dir(os.path.join(dockdir, "docking"))
                rows = []
                for ci, (p, cdir, lf) in enumerate(procs):
                    cs = os.path.join(cdir, "summary.csv")
                    if os.path.exists(cs):
                        with open(cs, encoding="utf-8-sig", newline="") as f:
                            for r in csv.DictReader(f):
                                r = dict(r)
                                r["name"] = "p%d_%s" % (ci, r.get("name", ""))
                                rows.append(r)
                    for pq in glob.glob(os.path.join(cdir, "docking", "*_out.pdbqt")):
                        shutil.move(pq, os.path.join(dockdir, "docking",
                                                     "p%d_%s" % (ci, os.path.basename(pq))))
                with open(summary, "w", encoding="utf-8-sig", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=["name", "smiles", "vina_score", "status"])
                    w.writeheader()
                    w.writerows(rows)
                n_ok = sum(1 for r in rows if r.get("status") == "ok")
                log("  合并完成: %d 条, ok %d" % (len(rows), n_ok))
        if not os.path.exists(summary):
            log("  无 summary.csv, 跳过")
            continue

        # 3) 口袋原子坐标(用于穿模判据)
        pg = pocket_geometry(os.path.join(outdir, pname))
        pocket_pos = pg["pos"] if pg["parse_ok"] else None
        if pocket_pos is None:
            log("  口袋解析失败, 跳过")
            continue

        # 4) 转位姿 + 治理
        logdir = ensure_dir(os.path.join(args.work, "conv_" + T))
        kept = 0
        reasons = {}
        with open(summary, encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                if row.get("status") != "ok":
                    reasons["dock_" + str(row.get("status"))] = \
                        reasons.get("dock_" + str(row.get("status")), 0) + 1
                    continue
                name = row["name"]
                src = os.path.join(dockdir, "docking", name + "_out.pdbqt")
                if not os.path.exists(src):
                    reasons["pose_missing"] = reasons.get("pose_missing", 0) + 1
                    continue
                tmp = os.path.join(logdir, name + ".sdf")
                if not pdbqt_to_sdf(src, tmp, logdir):
                    reasons["convert_fail"] = reasons.get("convert_fail", 0) + 1
                    continue
                ok, why, mind = validate_ligand(tmp, pocket_pos, args.mind_min)
                if not ok:
                    reasons[why.split("_")[0]] = reasons.get(why.split("_")[0], 0) + 1
                    skipped.append(("%s/%s" % (T, name), why))
                    continue
                lname = os.path.join("ligands", "lig_%s_%s.sdf" % (T, name))
                # 写出时剥掉显式氢并 sanitize, 与 build_gpcr_v3 产出的配体 SDF 保持一致
                m_clean = Chem.MolFromMolFile(tmp, sanitize=True, removeHs=True)
                if m_clean is None:
                    reasons["cleanfail"] = reasons.get("cleanfail", 0) + 1
                    continue
                w = Chem.SDWriter(os.path.join(outdir, lname))
                w.write(m_clean)
                w.close()
                entries.append((pname, lname, None, "0.0"))
                minds_all.append(mind)
                kept += 1
        per_target[T] = kept
        log("  保留 %d 个位姿 | 淘汰原因: %s" % (kept, dict(sorted(reasons.items()))))

    if len(entries) < 20:
        raise SystemExit("有效样本过少 (%d)" % len(entries))

    # 5) 划分: 按靶点分层随机(同一靶点内随机选验证集)
    random.seed(args.seed)
    by_pocket = {}
    for e in entries:
        by_pocket.setdefault(e[0], []).append(e)
    train, val = [], []
    for pk, es in sorted(by_pocket.items()):
        es = es[:]
        random.shuffle(es)
        n_val = max(1, int(round(len(es) * args.val_ratio)))
        val.extend(es[:n_val])
        train.extend(es[n_val:])
    split = {"train": [tuple(e[:2]) for e in train],
             "test": [tuple(e[:2]) for e in val]}

    with open(os.path.join(outdir, "index.pkl"), "wb") as f:
        pickle.dump(entries, f)
    import torch
    torch.save(split, os.path.join(outdir, "split_by_name.pt"))

    m = np.asarray(minds_all) if minds_all else np.asarray([0.0])
    with open(os.path.join(outdir, "build_report.txt"), "w", encoding="utf-8") as f:
        f.write("method=autodock_vina_docking_into_experimental_pockets\n")
        f.write("total=%d\ntrain=%d\nval=%d\n" % (len(entries), len(train), len(val)))
        f.write("pockets=%d (%s)\n" % (len(per_target), ",".join(sorted(per_target))))
        f.write("per_pocket=%s\n" % sorted(per_target.items()))
        f.write("exhaustiveness=%d\nmind_min=%.2f\n" % (args.exhaustiveness, args.mind_min))
        f.write("split_by=pair_stratified_by_target (seed=%d, val_ratio=%.3f)\n"
                % (args.seed, args.val_ratio))
        f.write("NOTE=验证集与训练集共享口袋: 评测靶点(A2A)本身在训练集中, 本划分度量的是"
                "对同一口袋生成新分子的能力(配体空间泛化), 不是跨靶点泛化\n")
        f.write("drug_like=MW[%.0f,%.0f] LogP<=%.1f QED>=%.2f rings>=%d\n"
                % (MW_LO, MW_HI, LOGP_MAX, QED_MIN, RINGS_MIN))
        f.write("pose_min_dist: n=%d median=%.3f p10=%.3f min=%.3f\n"
                % (len(minds_all), float(np.median(m)), float(np.percentile(m, 10)),
                   float(m.min())))
        f.write("skipped=%d\n" % len(skipped))
        for s, why in skipped[:40]:
            f.write("  skip %s: %s\n" % (s, why))

    log("=" * 78)
    log("完成: %d 对样本 (train %d / val %d), 口袋 %d 个, 每口袋 %s"
        % (len(entries), len(train), len(val), len(per_target),
           dict(sorted(per_target.items()))))
    log("位姿: 最近重原子距离 中位 %.3f A | P10 %.3f | 最小 %.3f (判据 >= %.2f)"
        % (float(np.median(m)), float(np.percentile(m, 10)), float(m.min()), args.mind_min))
    log("输出: %s" % outdir)


if __name__ == "__main__":
    main()
