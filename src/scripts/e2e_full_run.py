# -*- coding: utf-8 -*-
"""全盘端到端运行测试 (E2E)

依次执行并判定 8 个阶段, 任一阶段失败即记录并继续后续可独立执行的阶段:
  1 环境自检        依赖 + GPU
  2 配置生成        gen_sample_config (含 lam / diversity_w / relax_output)
  3 采样            sample_for_pdb (20/50/20, guided, relax_output=true)
  4 过滤入库        build_library -> compounds.csv / sdf
  5 分子对接        docking_pipeline (对新入库分子)
  6 新颖性          对采样产物的 Top 分子 vs A2A 已知集
  7 训练冒烟        train.py max_iters=2 (multitarget_v2, 验证 lmdb 重建/前向反向/ckpt)
  8 下游分析链      regen_top50 / selectivity / pareto / final_recommendation / roc / admet

用法: python scripts/e2e_full_run.py [--fast]
输出: outputs/e2e_<ts>/  (全过程日志与产物) + 控制台分阶段 PASS/FAIL
"""
import os, sys, time, glob, shutil, subprocess, argparse, csv

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import (CONFIGS, KNOWN_DRUGS, OUTPUTS, PYTHON, RESULTS, ROOT, TARGETS,
                   src_on_path)
src_on_path()
os.chdir(ROOT)
LIB = RESULTS
KNOWN = KNOWN_DRUGS
PDB = os.path.join(TARGETS, "4EIY_A2A受体.pdb")
CENTER = " -0.4,8.5,17.1"

stage_results = []


def run(cmd, log, timeout=None, cwd=ROOT):
    """执行子进程, stdout/stderr 重定向到文件 (避免管道沙箱问题)"""
    with open(log, "w", encoding="utf-8", errors="ignore") as lf:
        r = subprocess.run(cmd, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT, timeout=timeout)
    return r.returncode


def stage(name, ok, detail=""):
    stage_results.append((name, bool(ok), detail))
    print("[%s] %-22s %s" % ("PASS" if ok else "FAIL", name, detail), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="更小规模(10/50/20, 仅连通性验证)")
    args = ap.parse_args()
    # 默认与正式模板一致 (num_samples/beam_size/max_steps = 50/100/50):
    # max_steps 是安全上限而非目标步数, 设小会提前 break 导致 finished 偏少, 不代表链路故障
    ns, beam, steps = (10, 50, 20) if args.fast else (50, 100, 50)

    ts = time.strftime("%Y%m%d_%H%M%S")
    root = os.path.join(OUTPUTS, "e2e_" + ts)
    os.makedirs(root, exist_ok=True)
    print("E2E 开始: %s  (ns=%d beam=%d steps=%d)" % (root, ns, beam, steps), flush=True)

    # ---- 1 环境自检 ----
    code = ("import torch, rdkit, torch_geometric, torch_cluster, torch_scatter, lmdb, yaml, pptx;"
            "print('torch', torch.__version__, 'cuda', torch.cuda.is_available())")
    l1 = os.path.join(root, "01_env.log")
    with open(l1, "w", encoding="utf-8") as f:
        r = subprocess.run([PYTHON, "-c", code], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
    txt1 = open(l1, encoding="utf-8", errors="ignore").read()
    cuda_ok = "cuda True" in txt1 or "cuda True" in txt1.replace("cuda ", "cuda ")
    stage("1 环境自检", r.returncode == 0, txt1.strip().splitlines()[-1][:80] if txt1.strip() else "")

    # ---- 2 配置生成 ----
    cfg = os.path.join(root, "sample.yml")
    l2 = os.path.join(root, "02_gencfg.log")
    rc = run([PYTHON, "src/scripts/gen_sample_config.py", "--template", os.path.join(CONFIGS, "sample_for_pdb_guided_l3.yml"),
              "--out", cfg, "--seed", "2024", "--num-samples", str(ns), "--beam", str(beam),
              "--max-steps", str(steps), "--lam", "3.0", "--diversity-w", "0.5",
              "--guided", "1", "--relax-output", "1"], l2)
    ok2 = False
    if rc == 0 and os.path.exists(cfg):
        import yaml
        c = yaml.safe_load(open(cfg, encoding="utf-8-sig"))
        ok2 = (c["sample"]["num_samples"] == ns and c["sample"]["guided"]["lam"] == 3.0
               and c["sample"].get("relax_output") is True)
    stage("2 配置生成", ok2, "seed=2024 lam=3.0 relax_output=True")

    # ---- 3 采样 (开启构象精修, 顺带验证精修路径) ----
    sampdir = os.path.join(root, "sample_out")
    l3 = os.path.join(root, "03_sample.log")
    rc3 = run([PYTHON, "src/sample_for_pdb.py", "--pdb_path", PDB, "--center", CENTER,
               "--config", cfg, "--outdir", sampdir], l3, timeout=3600)
    sessions = sorted(glob.glob(os.path.join(sampdir, "sample_for_pdb_*")))
    n_smi, n_sdf, has_all = 0, 0, False
    if sessions:
        s = sessions[-1]
        p = os.path.join(s, "SMILES.txt")
        if os.path.exists(p):
            n_smi = len([x for x in open(p, encoding="utf-8", errors="ignore").read().splitlines() if x.strip()])
        n_sdf = len(glob.glob(os.path.join(s, "SDF", "*.sdf")))
        has_all = os.path.exists(os.path.join(s, "samples_all.pt"))
    stage("3 采样", rc3 == 0 and has_all, "finished=%d sdf=%d samples_all=%s" % (n_smi, n_sdf, has_all))

    # ---- 4 过滤入库 ----
    libout = os.path.join(root, "library")
    l4 = os.path.join(root, "04_library.log")
    rc4 = run([PYTHON, "src/scripts/build_library.py", "--runs", sampdir, "--library", libout], l4, timeout=1200)
    comp = os.path.join(libout, "compounds.csv")
    n_comp = 0
    if os.path.exists(comp):
        n_comp = sum(1 for _ in csv.DictReader(open(comp, encoding="utf-8-sig")))
    stage("4 过滤入库", rc4 == 0 and n_comp > 0, "compounds=%d" % n_comp)

    # ---- 5 对接 (对新入库分子) ----
    dock_ok, n_dock, rows = False, 0, []
    if n_comp > 0:
        smif = os.path.join(root, "lib_smiles.txt")
        with open(smif, "w", encoding="utf-8") as f:
            for row in csv.DictReader(open(comp, encoding="utf-8-sig")):
                f.write(row["smiles"] + "\n")
        dockout = os.path.join(root, "docking")
        l5 = os.path.join(root, "05_docking.log")
        rc5 = run([PYTHON, "src/scripts/docking_pipeline.py", "--receptor", PDB, "--center=" + CENTER.strip(),
                   "--smiles-file", smif, "--out", dockout], l5, timeout=3600)
        dsum = os.path.join(dockout, "summary.csv")
        if os.path.exists(dsum):
            rows = [r for r in csv.DictReader(open(dsum, encoding="utf-8-sig")) if r["status"] == "ok"]
            n_dock = len(rows)
            dock_ok = rc5 == 0 and n_dock > 0 and all(r["smiles"] for r in rows)
    stage("5 分子对接", dock_ok, "ok=%d" % n_dock)

    # ---- 6 新颖性 (采样产物 vs A2A 已知集) ----
    nov_ok, nov_detail = False, ""
    if n_dock > 0:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        from rdkit.Chem.Scaffolds import MurckoScaffold
        from rdkit.Chem import AllChem, DataStructs
        known = [s for s in open(os.path.join(KNOWN, "A2A_known_set.txt"), encoding="utf-8").read().splitlines() if s.strip()]
        kf = []
        for s in known:
            m = Chem.MolFromSmiles(s)
            if m:
                kf.append(AllChem.GetMorganFingerprint(m, 2))
        sims = []
        for r in rows[:50]:
            m = Chem.MolFromSmiles(r["smiles"])
            if m is None:
                continue
            sims.append(max(DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprint(m, 2), kf)))
        if sims:
            nov_ok = True
            nov_detail = "n=%d maxTanimoto 中位 %.3f" % (len(sims), sorted(sims)[len(sims)//2])
    stage("6 新颖性验证", nov_ok, nov_detail)

    # ---- 7 训练冒烟 (2 iter) ----
    tr_cfg = os.path.join(root, "train_smoke.yml")
    src = open(os.path.join(CONFIGS, "train_multitarget_v2.yml"), encoding="utf-8-sig").read()
    import re
    src = re.sub(r"max_iters: \d+", "max_iters: 2", src)
    src = re.sub(r"val_freq: \d+", "val_freq: 1", src)
    open(tr_cfg, "w", encoding="utf-8").write(src)
    logdir = os.path.join(root, "train_logs")
    os.makedirs(logdir, exist_ok=True)
    l7 = os.path.join(root, "07_train.log")
    rc7 = run([PYTHON, "train.py", "--config", tr_cfg, "--logdir", logdir], l7, timeout=3600)
    ckpts = glob.glob(os.path.join(logdir, "**", "checkpoints", "*.pt"), recursive=True)
    stage("7 训练冒烟(2 iter)", rc7 == 0 and len(ckpts) > 0, "checkpoints=%d" % len(ckpts))

    # ---- 8 下游分析链 (对全库已有数据) ----
    steps8 = [
        ("regen_top50", [PYTHON, "src/scripts/regen_top50.py", "A2A"]),
        ("selectivity", [PYTHON, "src/scripts/selectivity_analysis.py"]),
        ("pareto", [PYTHON, "src/scripts/pareto_analysis.py"]),
        ("final_rec", [PYTHON, "src/scripts/final_recommendation.py"]),
        ("roc", [PYTHON, "src/scripts/roc_decoys.py"]),
        ("admet", [PYTHON, "src/scripts/admet_screen.py"]),
    ]
    ok8, fails = True, []
    for nm, cmd in steps8:
        lg = os.path.join(root, "08_%s.log" % nm)
        rc = run(cmd, lg, timeout=1800)
        if rc != 0:
            ok8 = False
            fails.append(nm)
    stage("8 下游分析链(6脚本)", ok8, "全部成功" if ok8 else "失败: " + ",".join(fails))

    # ---- 汇总 ----
    print()
    print("=" * 70)
    npass = sum(1 for _, ok, _ in stage_results if ok)
    print("E2E 结果: %d/%d 阶段通过" % (npass, len(stage_results)))
    for nm, ok, d in stage_results:
        print("  [%s] %s %s" % ("PASS" if ok else "FAIL", nm, d))
    print("产物目录:", root)
    if npass != len(stage_results):
        print("失败阶段数:", len(stage_results) - npass)
    sys.exit(0 if npass == len(stage_results) else 1)


if __name__ == "__main__":
    main()
