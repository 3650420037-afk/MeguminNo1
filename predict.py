# -*- coding: utf-8 -*-
"""predict.py —— 一键生成最终结果文件(候选分子清单)。

流程(训练/设计 → 优化 → 推理/筛选):
    1 生成   在靶点口袋内生成候选分子            (sample_for_pdb, 口袋引导束搜索)
    2 过滤   药化规则 + 3D 构象生成, 入库         (build_library)
    3 打分   类药性/新颖性打分与排序              (本脚本; 指标口径与 src/evaluation 一致)
    4 落盘   results/results.csv + 结构文件 + 运行清单

用法示例
--------
    # 完整一键流程(A2A 靶点)
    python predict.py --target A2A

    # 复用已生成的结果, 只重跑筛选与排序
    python predict.py --skip-design --runs outputs/design/A2A_20260101_120000 --target A2A

    # 小规模自检(几分钟)
    python predict.py --target A2A --num-samples 10 --beam 50 --max-steps 20

输出
----
    results/results.csv              最终候选清单(UTF-8)
    results/structures/*.sdf         候选的三维结构(与清单中的“结构文件”对应)
    results/run_manifest.json        本次运行的参数/种子/环境/耗时(可溯源)
    <outdir>/library/compounds.csv   入库全量分子及其描述符(未截断)

字段说明见 README 的“最终结果文件”一节。
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from paths import (ROOT, CONFIGS, DATA, KNOWN_DRUGS, OUTPUTS, PRETRAINED,  # noqa: E402
                   RESULTS, ensure_dir, missing_hint)

MODEL_NAME = "Pocket2Mol (等变图神经网络) + 口袋引导束搜索"
MODEL_VERSION = "pretrained_Pocket2Mol.pt (ICML 2022 官方权重, 未微调)"
CODE_VERSION = "7-eonmol-competition-submission"

# 靶点 -> 已知活性分子参考集(用于新颖性评估); 缺失则该项留空
KNOWN_SETS = {
    "A2A": "A2A_known_set.txt",
    "B2AR": "beta-2_adrenergic_receptor_known_set.txt",
    "D3": "dopamine_D3_receptor_known_set.txt",
    "5HT2B": "5-HT2B_receptor_known_set.txt",
}

CSV_FIELDS = [
    "候选编号", "所属赛道", "候选SMILES", "结构文件",
    "关键预测指标_类药性综合分", "关键预测指标_QED", "关键预测指标_SA",
    "关键预测指标_新颖性最大Tanimoto",
    "关键预测指标_对接分Vina(kcal/mol)", "关键预测指标_ADMET分级",
    "关键预测指标_选择性SI百分位", "关键预测指标_六维综合Tier",
    "MW", "LogP", "TPSA", "LogS", "HBD", "HBA", "可旋转键", "环数", "Lipinski违反数",
    "ADMET警报", "可测性", "姿态能量差(kcal/mol)",
    "Murcko骨架", "同骨架候选数", "排序", "靶点",
    "模型名称", "模型版本", "代码版本", "随机种子", "运行编号", "备注",
]

# 已有分析结果 -> 清单列。键为库目录下的文件名, 值为 (源列, 目标列)。
# 这些文件由 src/scripts/ 下的分析脚本产出; 缺失时对应列留空, 不影响清单生成。
ENRICH_SPECS = [
    ("admet_screen.csv", {"grade": "关键预测指标_ADMET分级",
                          "alerts": "ADMET警报", "logS": "LogS"}),
    ("selectivity_full.csv", {"SI百分位点": "关键预测指标_选择性SI百分位",
                              "可测性": "可测性"}),
    ("pose_consistency.csv", {"energy_gap": "姿态能量差(kcal/mol)"}),
    ("final_recommendation.csv", {"tier": "关键预测指标_六维综合Tier",
                                  "vina_kcal": "关键预测指标_对接分Vina(kcal/mol)"}),
    ("top_candidates.csv", {"vina_kcal": "关键预测指标_对接分Vina(kcal/mol)"}),
]


def run(cmd, log_path=None):
    """执行子进程。日志写文件(避免管道缓冲问题), 同时回显末尾若干行。"""
    if log_path:
        ensure_dir(os.path.dirname(log_path))
        with open(log_path, "w", encoding="utf-8", errors="ignore") as lf:
            return subprocess.run(cmd, cwd=ROOT, stdout=lf,
                                  stderr=subprocess.STDOUT).returncode
    return subprocess.run(cmd, cwd=ROOT).returncode


def sdf_filename(canon_smiles, scaffold):
    """与 src/scripts/build_library.py 的 sdf_names() 保持一致的命名规则。"""
    mol_h = hashlib.md5(canon_smiles.encode("utf-8")).hexdigest()[:12]
    scaf_h = hashlib.md5(scaffold.encode("utf-8")).hexdigest()[:8] if scaffold else "chain"
    return "mol_%s_%s.sdf" % (scaf_h, mol_h)


def load_known_fps(target):
    """已知活性分子的 Morgan 指纹; 无参考集返回 None。"""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    fn = KNOWN_SETS.get(target)
    if not fn:
        return None, "该靶点无已知活性参考集"
    path = os.path.join(KNOWN_DRUGS, fn)
    if not os.path.exists(path):
        return None, "缺少参考集 %s" % os.path.relpath(path, ROOT)
    fps = []
    for line in open(path, encoding="utf-8", errors="ignore"):
        s = line.strip()
        if not s:
            continue
        m = Chem.MolFromSmiles(s)
        if m is not None:
            fps.append(AllChem.GetMorganFingerprint(m, 2))
    return (fps or None), ("参考集 %d 分子" % len(fps) if fps else "参考集为空")


def load_enrichment(library_dir):
    """载入库目录下已有的分析结果(对接/ADMET/选择性/姿态/六维推荐)。

    返回 (map, note): map 为 {规范化SMILES: {目标列: 值}}, note 说明用到了哪些文件。
    任一文件缺失或缺列都会被跳过, 因此该步骤永远是可选的。
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    out, used = {}, []
    for fn, spec in ENRICH_SPECS:
        p = os.path.join(library_dir, fn)
        if not os.path.exists(p):
            continue
        n, cols = 0, set()
        try:
            for row in csv.DictReader(open(p, encoding="utf-8-sig")):
                smi = (row.get("smiles") or "").strip()
                m = Chem.MolFromSmiles(smi)
                if m is None:
                    continue
                canon = Chem.MolToSmiles(m, isomericSmiles=True)
                d = out.setdefault(canon, {})
                for src, dst in spec.items():
                    v = row.get(src)
                    if v is not None and v != "":
                        # 先到先得: 排在前的文件优先(如 final_recommendation 优先于 top_candidates)
                        d.setdefault(dst, v)
                        cols.add(dst)
                n += 1
        except Exception as e:
            print("    [enrich] %s 读取失败, 跳过: %s" % (fn, e))
            continue
        if n:
            used.append("%s(%d 行, 补 %d 列)" % (fn, n, len(cols)))
    return out, ("; ".join(used) if used else "无(仅用本清单自带指标)")


def score_and_rank(library_dir, target, top, seed, run_id, track, results_dir,
                   rank_by="auto"):
    """读取 compounds.csv, 并入已有分析结果, 打分排序, 写出 results.csv。"""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem import DataStructs
    RDLogger.DisableLog("rdApp.*")

    comp = os.path.join(library_dir, "compounds.csv")
    if not os.path.exists(comp):
        raise SystemExit("[predict] 未找到入库结果 %s" % comp)
    rows = list(csv.DictReader(open(comp, encoding="utf-8-sig")))

    known_fps, known_note = load_known_fps(target)
    print("  新颖性参考: %s" % known_note, flush=True)
    enrich, enrich_note = load_enrichment(library_dir)
    print("  并入已有分析: %s" % enrich_note, flush=True)

    sdf_dir = os.path.join(library_dir, "sdf")
    out_sdf_dir = ensure_dir(os.path.join(results_dir, "structures"))

    recs = []
    for r in rows:
        smi = (r.get("smiles") or "").strip()
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        canon = Chem.MolToSmiles(m, isomericSmiles=True)
        qed = _f(r.get("qed"))
        sa = _f(r.get("sa"))
        mw = _f(r.get("mw"))
        logp = _f(r.get("logp"))
        hbd = _f(r.get("hbd"))
        hba = _f(r.get("hba"))
        rot = rdMolDescriptors.CalcNumRotatableBonds(m)
        # Lipinski 五规则: 与 src/evaluation/scoring_func.py 的 obey_lipinski 同口径
        passed = sum(int(x) for x in [
            (mw is not None and mw < 500),
            (hbd is not None and hbd <= 5),
            (hba is not None and hba <= 10),
            (logp is not None and -2 <= logp <= 5),
            rot <= 10,
        ])
        # 新颖性: 与已知活性分子的最大 Tanimoto 相似度(越低越新)
        nov = ""
        if known_fps:
            fp = AllChem.GetMorganFingerprint(m, 2)
            nov = round(max(DataStructs.BulkTanimotoSimilarity(fp, known_fps)), 3)
        # 类药性综合分 = QED + (1 - SA/10), 与生成阶段 utils/guidance.py 的
        # chem_score 同式, 保证"用来搜索的目标"与"用来排序的指标"一致
        chem = None
        if qed is not None and sa is not None:
            chem = qed + (1.0 - sa / 10.0)
        ex = enrich.get(canon, {})
        recs.append({
            "smiles": canon, "qed": qed, "sa": sa, "mw": mw, "logp": logp,
            "tpsa": _f(r.get("tpsa")), "hbd": hbd, "hba": hba, "rot": rot,
            "rings": _f(r.get("rings")), "scaffold": (r.get("murcko_scaffold") or ""),
            "violations": 5 - passed, "novelty": nov, "chem": chem,
            "pdbqt_ready": r.get("pdbqt_ready"), "extra": ex,
        })

    if not recs:
        raise SystemExit("[predict] 入库分子为空, 无法生成候选清单")

    # 同骨架候选数: 作为"该骨架被反复生成"的稳健性提示, 一并给出
    from collections import Counter
    scaf_cnt = Counter(x["scaffold"] for x in recs)

    # ---- 排序 ----
    vina_col = "关键预测指标_对接分Vina(kcal/mol)"
    tier_col = "关键预测指标_六维综合Tier"
    have_vina = any(x["extra"].get(vina_col) for x in recs)
    have_tier = any(x["extra"].get(tier_col) for x in recs)
    if rank_by == "auto":
        rank_by = "tier" if have_tier else ("vina" if have_vina else "druglikeness")
    if rank_by == "vina" and not have_vina:
        print("  [warn] 无对接分, 退回按类药性综合分排序")
        rank_by = "druglikeness"
    if rank_by == "tier" and not have_tier:
        print("  [warn] 无六维 Tier, 退回按类药性综合分排序")
        rank_by = "druglikeness"

    def key_vina(x):
        v = _f(x["extra"].get(vina_col))
        return (0, v) if v is not None else (1, 0.0)      # 缺分排最后

    def key_tier(x):
        t = str(x["extra"].get(tier_col) or "")
        # Tier1 优先, 其次 Tier2/3...; 无 Tier 者最后
        try:
            n = int("".join(ch for ch in t if ch.isdigit()) or 99)
        except ValueError:
            n = 99
        return (0 if t else 1, n, -(x["chem"] if x["chem"] is not None else -9))

    if rank_by == "vina":
        recs.sort(key=lambda x: (key_vina(x), -(x["chem"] if x["chem"] is not None else -9)))
        ranking_desc = ("关键预测指标_对接分Vina 升序(越负结合越强), 缺失者置底; "
                        "同分按类药性综合分 = QED + (1 − SA/10)")
    elif rank_by == "tier":
        recs.sort(key=key_tier)
        ranking_desc = "已有六维综合 Tier 优先(Tier1→Tier2→…), 组内按类药性综合分"
    else:
        recs.sort(key=lambda x: (-(x["chem"] if x["chem"] is not None else -9),
                                 -(x["qed"] if x["qed"] is not None else -9)))
        ranking_desc = "类药性综合分 = QED + (1 − SA/10) 降序, 同分按 QED"
    picked = recs[:top] if top and top > 0 else recs

    csv_path = os.path.join(results_dir, "results.csv")
    n_struct = 0
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for i, x in enumerate(picked, 1):
            src_sdf = os.path.join(sdf_dir, sdf_filename(x["smiles"], x["scaffold"]))
            struct_rel = ""
            if os.path.exists(src_sdf):
                import shutil
                dst = os.path.join(out_sdf_dir,
                                   "cand%03d_%s" % (i, os.path.basename(src_sdf)))
                if not os.path.exists(dst):
                    shutil.copyfile(src_sdf, dst)
                struct_rel = os.path.relpath(dst, ROOT).replace("\\", "/")
                n_struct += 1
            note = []
            if x["pdbqt_ready"] not in ("True", "true", "1"):
                note.append("3D构象未生成")
            if x["novelty"] != "" and x["novelty"] >= 0.8:
                note.append("与已知活性高度相似(Tanimoto>=0.8)")
            if x["violations"] > 1:
                note.append("Lipinski 违反>1")
            gray = str(x["extra"].get("关键预测指标_ADMET分级") or "").upper()
            if gray.startswith("C"):
                note.append("ADMET 分级 C(风险较高)")
            if rank_by == "vina" and not x["extra"].get(vina_col):
                note.append("无对接分")
            w.writerow({
                "候选编号": "C%04d" % i,
                "所属赛道": track,
                "候选SMILES": x["smiles"],
                "结构文件": struct_rel,
                "关键预测指标_类药性综合分": _s(x["chem"]),
                "关键预测指标_QED": _s(x["qed"]),
                "关键预测指标_SA": _s(x["sa"]),
                "关键预测指标_新颖性最大Tanimoto": _s(x["novelty"]),
                "关键预测指标_对接分Vina(kcal/mol)": x["extra"].get(vina_col, ""),
                "关键预测指标_ADMET分级": x["extra"].get("关键预测指标_ADMET分级", ""),
                "关键预测指标_选择性SI百分位": x["extra"].get("关键预测指标_选择性SI百分位", ""),
                "关键预测指标_六维综合Tier": x["extra"].get(tier_col, ""),
                "MW": _s(x["mw"]), "LogP": _s(x["logp"]), "TPSA": _s(x["tpsa"]),
                "LogS": x["extra"].get("LogS", ""),
                "HBD": _s(x["hbd"]), "HBA": _s(x["hba"]), "可旋转键": x["rot"],
                "环数": _s(x["rings"]), "Lipinski违反数": x["violations"],
                "ADMET警报": x["extra"].get("ADMET警报", ""),
                "可测性": x["extra"].get("可测性", ""),
                "姿态能量差(kcal/mol)": x["extra"].get("姿态能量差(kcal/mol)", ""),
                "Murcko骨架": x["scaffold"], "同骨架候选数": scaf_cnt[x["scaffold"]],
                "排序": i, "靶点": target,
                "模型名称": MODEL_NAME, "模型版本": MODEL_VERSION,
                "代码版本": CODE_VERSION, "随机种子": seed, "运行编号": run_id,
                "备注": "; ".join(note),
            })
    return len(recs), len(picked), n_struct, csv_path, known_note, ranking_desc, enrich_note


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _s(v, nd=3):
    return "" if v is None or v == "" else round(float(v), nd)


def main():
    ap = argparse.ArgumentParser(
        description="一键生成最终候选分子清单",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--target", default="A2A", help="靶点名(configs/targets.json)")
    ap.add_argument("--track", default="赛道 3 · 离子通道/GPCR 小分子药物虚拟筛选 · 子任务 2",
                    help="所属赛道(请按大赛官网赛道实施细则填写)")
    ap.add_argument("--num-samples", type=int, default=50)
    ap.add_argument("--beam", type=int, default=100)
    ap.add_argument("--max-steps", type=int, default=50)
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--lam", type=float, default=3.0, help="引导强度 λ")
    ap.add_argument("--diversity-w", type=float, default=0.5)
    ap.add_argument("--top", type=int, default=100, help="写入清单的候选数; 0 表示全部")
    ap.add_argument("--rank-by", choices=["auto", "druglikeness", "vina", "tier"],
                    default="auto",
                    help="排序依据: auto=有六维Tier则按Tier, 否则有对接分则按对接分, "
                         "再否则按类药性综合分")
    ap.add_argument("--skip-design", action="store_true", help="跳过生成, 复用 --runs")
    ap.add_argument("--runs", nargs="+", help="复用已有的采样会话目录(配合 --skip-design)")
    ap.add_argument("--library", help="直接对已有的化合物库目录(含 compounds.csv 与 sdf/)打分排序, "
                                      "跳过生成与入库两步(用于对既有结果重出清单)")
    ap.add_argument("--outdir", default=os.path.join(OUTPUTS, "predict"))
    ap.add_argument("--results-dir", default=RESULTS)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    miss = missing_hint([("预训练权重", PRETRAINED)])
    if miss and not (args.skip_design or args.library):
        raise SystemExit("[predict] 缺少必需文件:\n" + miss +
                         "\n请按 README 的“模型权重”一节准备。")

    run_id = "%s_%s" % (args.target, time.strftime("%Y%m%d_%H%M%S"))
    t0 = time.time()
    ensure_dir(args.outdir)
    ensure_dir(args.results_dir)
    run_dir = ensure_dir(os.path.join(args.outdir, run_id))

    print("=" * 76)
    print("7-eonmol predict | 靶点=%s | 运行编号=%s" % (args.target, run_id))
    print("=" * 76, flush=True)

    # ---- 1 生成 ----
    if args.library:
        lib_dir = os.path.abspath(args.library)
        if not os.path.exists(os.path.join(lib_dir, "compounds.csv")):
            raise SystemExit("[predict] %s 下没有 compounds.csv" % lib_dir)
        run_dirs, gen_rc, gen_sec = [], 0, 0.0
        print("[1/4] 生成: 跳过(直接使用已有库 %s)" % os.path.relpath(lib_dir, ROOT), flush=True)
    elif args.skip_design:
        if not args.runs:
            raise SystemExit("[predict] --skip-design 需要同时给出 --runs <会话目录>")
        run_dirs = [os.path.abspath(p) for p in args.runs]
        for d in run_dirs:
            if not os.path.isdir(d):
                raise SystemExit("[predict] 采样会话目录不存在: %s" % d)
        gen_rc, gen_sec = 0, 0.0
        print("[1/4] 生成: 跳过(复用 %d 个已有会话)" % len(run_dirs), flush=True)
    else:
        print("[1/4] 生成: 靶点 %s, %d 分子 / 束宽 %d / 种子 %d"
              % (args.target, args.num_samples, args.beam, args.seed), flush=True)
        ts = time.strftime("%Y%m%d_%H%M%S")
        design_out = ensure_dir(os.path.join(run_dir, "design"))
        cmd = [sys.executable, os.path.join(ROOT, "design.py"),
               "--target", args.target, "--num-samples", str(args.num_samples),
               "--beam", str(args.beam), "--max-steps", str(args.max_steps),
               "--seed", str(args.seed), "--lam", str(args.lam),
               "--diversity-w", str(args.diversity_w),
               "--device", args.device, "--outdir", design_out]
        t1 = time.time()
        gen_rc = run(cmd, os.path.join(run_dir, "01_design.log"))
        gen_sec = time.time() - t1
        if gen_rc != 0:
            raise SystemExit("[predict] 生成阶段失败(返回码 %d), 见 %s"
                             % (gen_rc, os.path.join(run_dir, "01_design.log")))
        run_dirs = [os.path.join(design_out, d) for d in sorted(os.listdir(design_out))
                    if os.path.isdir(os.path.join(design_out, d))]
        if not run_dirs:
            raise SystemExit("[predict] 生成阶段未产生会话目录")

    # ---- 2 过滤入库 ----
    if args.library:
        lib_sec = 0.0
        print("[2/4] 过滤入库: 跳过(使用已有库)", flush=True)
    else:
        lib_dir = os.path.join(run_dir, "library")
        print("[2/4] 过滤入库: 药化规则 + 3D 构象生成 ...", flush=True)
        t2 = time.time()
        rc = run([sys.executable, os.path.join(ROOT, "src", "scripts", "build_library.py"),
                  "--runs"] + run_dirs + ["--library", lib_dir],
                 os.path.join(run_dir, "02_library.log"))
        lib_sec = time.time() - t2
        if rc != 0:
            raise SystemExit("[predict] 入库失败(返回码 %d), 见 %s"
                             % (rc, os.path.join(run_dir, "02_library.log")))

    # ---- 3 打分排序 ----
    print("[3/4] 打分排序 ...", flush=True)
    t3 = time.time()
    n_all, n_top, n_struct, csv_path, known_note, ranking_desc, enrich_note = score_and_rank(
        lib_dir, args.target, args.top, args.seed, run_id, args.track, args.results_dir,
        rank_by=args.rank_by)
    score_sec = time.time() - t3

    # ---- 4 运行清单 ----
    import torch
    manifest = {
        "run_id": run_id, "target": args.target, "track": args.track,
        "stage": "predict",
        "params": {"num_samples": args.num_samples, "beam_size": args.beam,
                   "max_steps": args.max_steps, "seed": args.seed, "lam": args.lam,
                   "diversity_w": args.diversity_w, "top": args.top,
                   "skip_design": bool(args.skip_design)},
        "inputs": {"design_runs": [os.path.relpath(d, ROOT) for d in run_dirs],
                   "library_dir": os.path.relpath(lib_dir, ROOT),
                   "checkpoint": os.path.relpath(PRETRAINED, ROOT)},
        "result": {"library_molecules": n_all, "listed_candidates": n_top,
                   "structure_files": n_struct,
                   "results_csv": os.path.relpath(csv_path, ROOT),
                   "ranking": ranking_desc,
                   "enrichment": enrich_note,
                   "novelty_reference": known_note},
        "timing_sec": {"design": round(gen_sec, 1), "library": round(lib_sec, 1),
                       "scoring": round(score_sec, 1),
                       "total": round(time.time() - t0, 1)},
        "env": {"python": sys.version.split()[0], "torch": torch.__version__,
                "cuda": torch.cuda.is_available(),
                "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
                "platform": platform.platform()},
        "model": {"name": MODEL_NAME, "version": MODEL_VERSION,
                  "code_version": CODE_VERSION},
    }
    mpath = os.path.join(args.results_dir, "run_manifest.json")
    with open(mpath, "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=1)

    print("[4/4] 完成")
    print("=" * 76)
    print("入库分子 %d 个 -> 清单候选 %d 个 (含结构文件 %d 个)" % (n_all, n_top, n_struct))
    print("候选清单: %s" % os.path.relpath(csv_path, ROOT))
    print("运行清单: %s" % os.path.relpath(mpath, ROOT))
    print("耗时: 生成 %.1fs / 入库 %.1fs / 打分 %.1fs / 合计 %.1fs"
          % (gen_sec, lib_sec, score_sec, time.time() - t0))
    print("=" * 76)
    return 0


if __name__ == "__main__":
    sys.exit(main())
