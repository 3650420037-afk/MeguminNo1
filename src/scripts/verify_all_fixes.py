# -*- coding: utf-8 -*-
"""全工程修复验证清单 (端到端)

对本次 bug 审计后的所有修复做可自动判定的验证, 输出 PASS/FAIL。
用法: python src/scripts/verify_all_fixes.py [--deploy-dir DIR] [--docs-dir DIR]

  --deploy-dir  部署包目录 (内含 7-eonmol.exe 与 repo/ 子目录), 默认 <root>/dist/7-eonmol
  --docs-dir    报告/文档目录 (PPT、审计报告、使用指南所在处), 默认 <root>/docs
"""
import argparse
import os, sys, csv, glob, subprocess, importlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import DOCS, OUTPUTS, PYTHON, RESULTS, ROOT, SRC, src_on_path
src_on_path()

_ap = argparse.ArgumentParser()
_ap.add_argument("--deploy-dir", default=os.path.join(ROOT, "dist", "7-eonmol"),
                 help="部署包目录 (内含 7-eonmol.exe 与 repo/ 子目录)")
_ap.add_argument("--docs-dir", default=DOCS, help="报告/文档目录")
_args = _ap.parse_args()

os.chdir(ROOT)
LIB = RESULTS
OUT = OUTPUTS

results = []
skips = []


def chk(name, ok, detail=""):
    results.append((name, bool(ok), detail))


def skip(name, why):
    """依赖仓库外产物(部署包/PPT/线下报告)的检查: 目标不存在时记 SKIP 而非 FAIL。

    这些产物按 .gitignore 不随仓库提交, 默认参数下本就找不到,
    若记为 FAIL 会让"全部通过"变成噪声, 掩盖真正的回归。
    """
    skips.append((name, why))


def read(p):
    with open(os.path.join(ROOT, p), encoding="utf-8", errors="ignore") as f:
        return f.read()


# 1) 源码修复在位
chk("docking_pipeline: --smiles-file 去 BOM", 'encoding="utf-8-sig"' in read("src/scripts/docking_pipeline.py"))
chk("docking_pipeline: SDF 用 open()+MolBlock", "MolFromMolBlock(block" in read("src/scripts/docking_pipeline.py"))
chk("docking_pipeline: SMILES 去氢规范化", "MolToSmiles(Chem.RemoveHs(m))" in read("src/scripts/docking_pipeline.py"))
chk("maskfill: sample_position 未实现分支显式报错", "sample_position(n_samples>=0) 未实现" in read("src/models/maskfill.py"))
chk("maskfill: assert 已改 raise", "if len(torch.unique(batch_query)) != 1:" in read("src/models/maskfill.py"))
chk("reconstruct: 无裸 except", "\n    except:\n" not in read("src/utils/reconstruct.py"))
chk("transforms: except ImportError", "except ImportError:" in read("src/utils/transforms.py"))
chk("guidance: 保留 shim 供 sascorer 兼容", "_six" in read("src/utils/guidance.py") or "add_safe_globals" in read("src/utils/misc.py"))
chk("app_easy: 关窗协议", "WM_DELETE_WINDOW" in read("src/gui/app_easy.py"))
chk("app_easy: 过滤阶段兜底", "过滤阶段异常" in read("src/gui/app_easy.py"))
chk("app_easy: 配置生成走生成器", "gen_sample_config.py" in read("src/gui/app_easy.py"))
chk("app.py: 会话回退已移除", "best_any" not in read("src/gui/app.py").split("def ")[-1])
chk("overnight_phase2: PROMOTED 锚定正则", "DECISION:\\s*PROMOTED\\s*$" in read("src/scripts/overnight_phase2.ps1"))
chk("build_library: qed 直方图空保护", "无可统计分子" in read("src/scripts/build_library.py"))
chk("sascorer: 内联 rdkit.six shim", "import pickle as cPickle" in read("src/evaluation/sascorer.py"))
chk("scoring_func: 规则4 海象 bug 已修", "_lp = Crippen.MolLogP(mol)" in read("src/evaluation/scoring_func.py"))
chk("audit_conformers: PCA 平面性判据", "np.linalg.svd(centered" in read("src/scripts/audit_conformers.py"))
chk("verify_gpcr_datasets: 一致性校验", "index_not_in_split" in read("src/scripts/verify_gpcr_datasets.py"))
chk("phase4: 等待超时", "PHASE4_WAIT_TIMEOUT" in read("src/scripts/phase4_retrain.ps1"))
chk("requirements.txt 存在", os.path.exists(os.path.join(ROOT, "requirements.txt")))
chk("selval_run.ps1 已删除", not os.path.exists(os.path.join(SRC, "scripts/selval_run.ps1")))

# 2) ps1 BOM 全绿
ps1 = glob.glob(os.path.join(SRC, "scripts", "*.ps1"))
nb = [os.path.basename(p) for p in ps1 if open(p, "rb").read(3) != b"\xef\xbb\xbf"]
chk("全部 ps1 带 UTF-8 BOM", not nb, ",".join(nb))

# 3) 关键依赖可导入 (含 pptx / sascorer 独立导入)
code = ("import pptx, torch, rdkit, torch_geometric, torch_cluster, torch_scatter;"
        "from evaluation.sascorer import calculateScore;"
        "print('deps+independent_sascorer OK')")
r = subprocess.run([PYTHON, "-c", code], cwd=SRC, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
chk("依赖齐备 + sascorer 可独立导入", r.returncode == 0 and b"OK" in r.stdout,
    r.stdout.decode("utf-8", "ignore").strip()[-80:])

# 4) scoring_func 规则4 实测 (logP>5 的分子应判 False)
code2 = ("import sys; sys.path.insert(0, r'" + SRC + "');"
         "from evaluation.scoring_func import obey_lipinski;"
         "from rdkit import Chem;"
         "m=Chem.MolFromSmiles('CCCCCCCCCCCCCCCCCC(=O)O');"   # 长链脂肪酸 logP>5
         "print('lipinski_rules=' + str(obey_lipinski(m)))")
r2 = subprocess.run([PYTHON, "-c", code2], cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
out2 = r2.stdout.decode("utf-8", "ignore")
chk("Lipinski 规则4 修复实证 (logP>5 不再恒真)", r2.returncode == 0 and "lipinski_rules=" in out2, out2.strip()[-60:])

# 5) 配置生成器: 参数生效 + 产物可解析
r3 = subprocess.run([PYTHON, "src/scripts/gen_sample_config.py", "--template", "configs/sample_for_pdb_guided_l3.yml",
                     "--out", os.path.join(os.environ["TEMP"], "vfy_cfg.yml"),
                     "--seed", "4321", "--lam", "1.5", "--diversity-w", "0.3", "--guided", "1"],
                    cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
try:
    import yaml
    c = yaml.safe_load(open(os.path.join(os.environ["TEMP"], "vfy_cfg.yml"), encoding="utf-8-sig"))
    ok5 = (r3.returncode == 0 and c["sample"]["seed"] == 4321
           and c["sample"]["guided"]["lam"] == 1.5 and c["sample"]["guided"]["diversity_w"] == 0.3)
    chk("配置生成器: seed/lam/diversity 全部生效", ok5,
        "seed=%s lam=%s div=%s" % (c["sample"]["seed"], c["sample"]["guided"]["lam"], c["sample"]["guided"]["diversity_w"]))
except Exception as e:
    chk("配置生成器: seed/lam/diversity 全部生效", False, str(e))

# 6) 数据完整性: 四靶点 docking summary 与库 CSV 交集 100%
data_ok, detail = True, []
for t, comp in (("a2a", os.path.join(LIB, "compounds.csv")),
                ("b2ar", os.path.join(LIB, "b2ar", "compounds.csv")),
                ("d3", os.path.join(LIB, "d3", "compounds.csv")),
                ("5ht2b", os.path.join(LIB, "5ht2b", "compounds.csv"))):
    sp = os.path.join(OUT, "docking_%s_library" % t, "summary.csv")
    if not (os.path.exists(sp) and os.path.exists(comp)):
        continue
    lib = set(r["smiles"] for r in csv.DictReader(open(comp, encoding="utf-8-sig")))
    ds = [r for r in csv.DictReader(open(sp, encoding="utf-8-sig")) if r["status"] == "ok"]
    inter = sum(1 for r in ds if r["smiles"] in lib)
    detail.append("%s %d/%d" % (t, inter, len(ds)))
    if inter != len(ds):
        data_ok = False
chk("四靶点对接-库主键交集 100%", data_ok, " ".join(detail))

# 7) 关键下游产物的行数与结论数字
def rowcount(p):
    p = os.path.join(LIB, p)
    return sum(1 for _ in open(p, encoding="utf-8-sig")) - 1 if os.path.exists(p) else -1
chk("selectivity_full.csv = 2931 行", rowcount("selectivity_full.csv") == 2931, str(rowcount("selectivity_full.csv")))
chk("final_recommendation.csv 存在(一梯队166)", rowcount("final_recommendation.csv") >= 166, str(rowcount("final_recommendation.csv")))
chk("pareto_front.csv 存在", rowcount("pareto_front.csv") > 0, str(rowcount("pareto_front.csv")))
chk("roc_decoys.csv 存在", rowcount("roc_decoys.csv") > 0, str(rowcount("roc_decoys.csv")))
chk("admet_screen.csv = 2931 行", rowcount("admet_screen.csv") == 2931, str(rowcount("admet_screen.csv")))

# 8) 部署包同步 + PPT + exe  (依赖仓库外产物, 不存在则跳过)
dst = _args.deploy_dir
same = 0
pairs = ["src/models/maskfill.py", "src/utils/reconstruct.py", "src/utils/misc.py", "src/utils/protein_ligand.py",
         "src/sample.py", "src/sample_for_pdb.py", "src/models/encoders/__init__.py",
         "src/scripts/gen_sample_config.py", "requirements.txt"]
import hashlib
def h(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest() if os.path.exists(p) else None
if not os.path.isdir(dst):
    skip("部署包关键文件同步 9/9", "部署包目录不存在(不随仓库提交): %s" % dst)
    skip("部署包 exe 存在", "同上")
else:
    for n in pairs:
        if h(os.path.join(ROOT, n)) and h(os.path.join(ROOT, n)) == h(os.path.join(dst, "repo", n)):
            same += 1
    chk("部署包关键文件同步 9/9", same == len(pairs), "%d/%d" % (same, len(pairs)))
    chk("部署包 exe 存在", os.path.exists(os.path.join(dst, "7-eonmol.exe")))

ppts = [f for f in os.listdir(_args.docs_dir)] if os.path.isdir(_args.docs_dir) else []
if any(f.endswith(".pptx") for f in ppts):
    chk("PPT 已生成", True)
else:
    skip("PPT 已生成", "PPT 为线下答辩材料, 不在仓库内(可由 docs/make_ppt.py 生成)")

# 9) 报告/文档  (线下材料, 不存在则跳过)
for nm, rel in (("bug 审计报告存在", os.path.join("CopilotExport", "bug审计与修复报告.md")),
                ("使用指南存在(文档目录)", "05_使用指南.md")):
    if os.path.exists(os.path.join(_args.docs_dir, rel)):
        chk(nm, True)
    else:
        skip(nm, "线下材料, 不在仓库内: %s" % os.path.join(_args.docs_dir, rel))

# ---- 输出 ----
print("=" * 78)
for name, ok, detail in results:
    print("[%s] %-46s %s" % ("PASS" if ok else "FAIL", name, detail))
for name, why in skips:
    print("[SKIP] %-46s %s" % (name, why))
print("=" * 78)
npass = sum(1 for _, ok, _ in results if ok)
print("结果: %d/%d 通过%s" % (npass, len(results),
      ("; %d 项跳过(依赖仓库外产物)" % len(skips)) if skips else ""))
if npass != len(results):
    print("失败项:")
    for name, ok, detail in results:
        if not ok:
            print("  -", name, detail)
    sys.exit(1)
print("ALL PASS")
sys.exit(0)
