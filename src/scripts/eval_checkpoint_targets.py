# -*- coding: utf-8 -*-
"""一个检查点 x 多靶点 x 多种子 的批量评测编排器(只编排, 自己不加载模型)。

为什么需要它
------------
`select_ckpt_by_generation.py` 一次只看一个靶点, `ab_docking_compare.py` 一次也只看
一个靶点。要给一个候选权重出一张"多靶点体检表", 只能手工反复拼命令 —— 既容易把
采样参数(阈值/规模/种子)写岔, 也容易漏掉某一轴。本脚本把这件事固定成一条可复现
的流水线, 供后续批量评测(例如 12 个靶点全跑)使用。

流程(对每个 靶点 x 种子)
------------------------
a) `src/scripts/select_ckpt_by_generation.py`(**子进程**)采样, 产物保留到
   `<keep-dir>/<靶点>/<种子>/`, 读它的 `--out` CSV 取生成级口径:
   n_finished / qed_med / sa_med / n_scaffold / clash_rate / mind_med;
b) `src/scripts/arm_gen_metrics.py`(**子进程**)算尺寸与 strict 通过率:
   mw_med / ha_med / strict_rate(该脚本只读 SMILES, 纯 CPU);
c) 除非 `--skip-dock`, 再调 `src/scripts/ab_docking_compare.py`(**子进程**)做
   **与官方权重的同种子成对 A/B**: baseline = `models/pretrained_Pocket2Mol.pt`,
   第二臂 = `--ckpt`, 产出该 (靶点, 种子) 的对接 CSV。
   拿到官方臂的采样产物后再跑一次 (b)(仍是纯 CPU 读 CSV), 于是多出
   strict 占官方比例、以及"尺寸差"两个对照量。

产物与列的口径
--------------
`--out` CSV 每行 = 靶点/种子/完成数/QED/SA/骨架/穿模/MW中位/重原子中位/strict率;
有对接时再附 Vina中位 / Δ vs 官方 / 重原子 / LE / 是否退化。
注意 `ha_med` 是**生成口径**的重原子中位, `dock_ha_med` 是**对接口径**(真正对接
成功那批分子的重原子中位), 两者不可互换。

目录布局(每个 (靶点, 种子) 一个目录, 便于单独复核)
--------------------------------------------------
    <keep-dir>/<靶点>/<种子>/
        select_gen_metrics.csv          a) 的产物指标(子脚本 --out)
        arm_gen_metrics.csv             b) 的尺寸/strict
        ab_docking.csv                  c) 的成对对接结果
        arm_gen_metrics_vs_official.csv b') 与官方臂的尺寸/strict 对照
        eval_meta.json                  该目录产物对应的采样参数(复用核对用)
        select.log / arm_gen_metrics.log / ab_docking.log
        <检查点名>/samples/<会话>/      采样产物(SMILES.txt / SDF/)
        ab_dock/                       A/B 两臂各自的采样与对接产物

安全与纪律(写进代码, 不是口号)
------------------------------
1. **绝不并发采样**: (靶点, 种子) 严格串行; `--parallel` 只透传给
   `ab_docking_compare.py` 作为 **Vina 对接**的分块并发。理由: 采样要占 GPU, 评测机
   上往往还有训练在跑, 并发采样会 OOM —— 省下的时间远不够赔一次重跑。
2. 任何子进程返回码非 0 **只记录、不中断**(`run()` 连启动异常都兜住), 该行写进
   `errors` 列后继续下一个 (靶点, 种子)。
3. **不重跑已有产物**: `--keep-dir` 下已有该 (靶点, 种子) 的指标 CSV 就直接复用;
   采样产物在而指标 CSV 丢了, 也只在 **CPU 上离线复算**(口径与子脚本同源, 见
   `offline_gen_row`), 不再占 GPU。复用时会核对 `eval_meta.json` 里的采样参数,
   不一致就打警告并在行内留痕。
4. **单种子标注 PROVISIONAL**: `--seeds` 少于 2 个时所有行 status=PROVISIONAL,
   不得据此晋级(依据: `ab_docking_compare.py` 实测种子间 Vina 中位噪声约 0.7
   kcal/mol, 远大于同臂重复运行的 ~0.14)。
5. **尺寸差 > 15% 警告**: 与官方臂的重原子/MW 中位相差超过 15% 时打印警告 ——
   两臂分子大小差得多时, 原始 Vina 分的差异主要来自尺寸而非亲和力, 必须以 LE 复核。
6. **不调用 git**, 不修改任何既有文件: 本脚本只写 `--keep-dir` 与 `--out` 下的文件。
7. 退出码: 0 = 至少有一行可用; 1 = 全部 (靶点, 种子) 都不可用。

开销(先算账再跑)
----------------
每个 (靶点, 种子) 最多 3 次采样(a 步 1 次 + A/B 两臂各 1 次)+ 1 轮对接。
只看生成级就加 `--skip-dock`, 采样次数从 3 降到 1。

用法
----
    # 全量: 4 靶点 x 2 种子, 生成 + 与官方权重的成对对接
    python src/scripts/eval_checkpoint_targets.py \\
        --ckpt models/7-eonmol_ft_gpcr_v3.pt \\
        --targets A2A B2AR D3 5HT2B --num-samples 50 --beam 50 --max-steps 40 \\
        --seeds 2024 2025 --parallel 6

    # 只做生成级(不跑对接, 采样次数 3 -> 1):
    python src/scripts/eval_checkpoint_targets.py --ckpt <ckpt> --skip-dock

    # 先看会跑哪些子命令(不执行任何子进程, 不建目录、不写任何文件):
    python src/scripts/eval_checkpoint_targets.py --ckpt <ckpt> --dry-run

⚠ 本脚本会通过子进程调用 `src/sample_for_pdb.py`(**需要 CUDA**), 因此**不要在 GPU
   正在训练时执行**。`--help` 与 `--dry-run` 都不加载 torch/rdkit, 可安全运行。
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, CONFIGS, OUTPUTS, ensure_dir  # noqa: E402

# ---------------------------------------------------------------- 子脚本
SELECT = os.path.join(ROOT, "src", "scripts", "select_ckpt_by_generation.py")
ARM = os.path.join(ROOT, "src", "scripts", "arm_gen_metrics.py")
AB = os.path.join(ROOT, "src", "scripts", "ab_docking_compare.py")

# A/B 的 baseline 固定为官方权重(与 ab_docking_compare.py 的成对口径一致)
BASELINE_NAME = "pretrained_Pocket2Mol.pt"
BASELINE = os.path.join(ROOT, "models", BASELINE_NAME)

# 判据容差: 与 ab_docking_compare.py / arm_gen_metrics.py 保持同源, 不在这里另立标准。
TOL_N = -5                # 完成/对接分子数允许的最大退步
TOL_QED = -0.05
TOL_SA = 0.30
TOL_SCAF = 0.8
TOL_VINA = 0.20           # Vina 越低越强, 允许 0.2 kcal/mol 噪声
SIZE_TOL = 0.15           # 尺寸差 >15% 必须警告
STRICT_RATIO_MIN = 0.90   # 同 arm_gen_metrics.STRICT_RATIO_MIN
SIZE_MW_MIN = 340.0       # 同 arm_gen_metrics.SIZE_MW_MIN(预登记下限)
SIZE_HA_MIN = 25.0        # 同 arm_gen_metrics.SIZE_HA_MIN(预登记下限)
CLASH_MAX = 0.5           # 同 select_ckpt_by_generation.CLASH_MAX

SELECT_COLS = ["n_finished", "qed_med", "sa_med", "n_scaffold", "clash_rate", "mind_med"]
ARM_COLS = ["n_finished", "mw_med", "ha_med", "strict_rate"]
AB_COLS = ["n_smiles", "n_docked", "vina_med"]

# 汇总 CSV 的固定列序(某列缺失也照写, 便于增量落盘以及多检查点结果拼接)
SUMMARY_COLS = [
    "target", "seed", "ckpt", "status",
    "n_finished", "qed_med", "sa_med", "n_scaffold", "clash_rate", "mind_med",
    "mw_med", "ha_med", "strict_rate",
    "dock_vina_med", "dock_delta_vs_official", "dock_ha_med", "dock_le_med",
    "dock_n_smiles", "dock_n_docked", "dock_frac_le10",
    "ha_ratio_vs_official", "mw_ratio_vs_official", "strict_ratio_vs_official",
    "degraded", "reasons",
    "gen_source", "gen_rc", "arm_rc", "dock_rc", "errors",
]


# ---------------------------------------------------------------- 小工具
def log(msg):
    print(msg, flush=True)


def rel(p):
    try:
        return os.path.relpath(p, ROOT)
    except Exception:  # noqa: BLE001
        return p


def stem(path):
    """检查点路径 -> 产物目录名。

    必须与两个子脚本一致: select 用 os.path.splitext(basename)[0],
    ab_docking_compare.gen_one 用 basename[:-3](仅 *.pt)。对 *.pt 二者等价。
    """
    b = os.path.basename(path or "")
    return b[:-3] if b.endswith(".pt") else os.path.splitext(b)[0]


def resolve_path(p):
    """检查点路径: 绝对路径原样; 相对路径先按当前目录, 再按仓库根(两种都常见)。"""
    if os.path.isabs(p):
        return p
    a = os.path.abspath(p)
    if os.path.exists(a):
        return a
    b = os.path.join(ROOT, p)
    return b if os.path.exists(b) else a


def root_join(p):
    """输出类路径: 相对路径一律按仓库根解析(与 ab_docking_compare.py 同约定)。"""
    return p if os.path.isabs(p) else os.path.join(ROOT, p)


def fnum(row, key):
    """数值列 -> float; 空/None/NA 一律 None(不抛异常)。"""
    if not row:
        return None
    v = row.get(key)
    if v is None or v == "" or v == "None":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def inum(row, key):
    """计数列 -> int(整数值时); 否则退回 float/None。只影响 CSV 观感。"""
    v = fnum(row, key)
    if v is None:
        return None
    return int(v) if float(v).is_integer() else v


def read_csv(path, need_cols, tag=""):
    """读指标 CSV(utf-8-sig)。列不全/空/读失败 -> None, 即"没有可用产物"。

    纪律 3 的第一道闸: 只有通过列名校验的 CSV 才允许被当成"已有产物"复用, 免得把
    别家脚本的异形 CSV(例如旧版 ab_docking 缺 ha_med/le_med)当成本次口径。
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        with open(path, encoding="utf-8-sig", newline="") as f:
            rows = [r for r in csv.DictReader(f)]
    except Exception as e:  # noqa: BLE001
        log("  [警告] %s 读取失败(%r), 视为无产物" % (rel(path), e))
        return None
    rows = [r for r in rows if r]
    if not rows:
        return None
    missing = [c for c in need_cols if c not in rows[0]]
    if missing:
        log("  [警告] %s 缺列 %s(%s) -> 视为无产物, 需重算" % (rel(path), missing, tag))
        return None
    return rows


def run(cmd, log_path, tag):
    """执行子进程, stdout/stderr 追加进 log_path; **永不抛异常**, 返回 returncode。

    纪律 2: 非 0 只记录、不中断 —— "一个靶点失败就整批崩掉"会白烧掉其它靶点已经
    等到的 GPU 时间。启动异常(解释器/脚本找不到)返回 -1。
    """
    ensure_dir(os.path.dirname(log_path))
    t0 = time.time()
    try:
        with open(log_path, "a", encoding="utf-8", errors="ignore") as lf:
            lf.write("\n$ %s\n" % " ".join(cmd))
            lf.flush()
            rc = subprocess.run(cmd, cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT).returncode
    except Exception as e:  # noqa: BLE001
        log("  [警告] %s 子进程启动失败: %r" % (tag, e))
        try:
            with open(log_path, "a", encoding="utf-8", errors="ignore") as lf:
                lf.write("\n[启动失败] %r\n" % (e,))
        except Exception:  # noqa: BLE001
            pass
        return -1
    if rc != 0:
        log("  [警告] %s 子进程返回码 %d (已记录, 继续); 日志 %s" % (tag, rc, rel(log_path)))
    else:
        log("  [ok] %s 完成 (%.0f s)" % (tag, time.time() - t0))
    return rc


def latest_session(parent):
    """<parent>/samples/<session>/ 中最新的采样会话目录(会话名带时间戳)。"""
    cands = sorted(d for d in glob.glob(os.path.join(parent, "samples", "*"))
                   if os.path.isdir(d))
    return cands[-1] if cands else None


def load_meta(base):
    try:
        with open(os.path.join(base, "eval_meta.json"), encoding="utf-8") as f:
            return json.load(f)
    except Exception:  # noqa: BLE001
        return None


def save_meta(base, d):
    if not d:
        return
    try:
        with open(os.path.join(base, "eval_meta.json"), "w", encoding="utf-8") as f:
            json.dump(d, f, ensure_ascii=False, indent=1)
    except Exception as e:  # noqa: BLE001
        log("  [提示] eval_meta.json 写不了: %r" % (e,))


def write_rows(path, rows):
    """按 rows[0] 的键写 CSV(用于把离线复算结果落盘成子脚本同款 schema)。"""
    ensure_dir(os.path.dirname(path))
    cols = list(rows[0].keys())
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_summary(rows, out):
    """整表重写(行数很少); 每完成一个 (靶点,种子) 就写一次 -> 中途被打断也留痕。"""
    ensure_dir(os.path.dirname(out))
    tmp = out + ".tmp"
    with open(tmp, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, out)


# ---------------------------------------------------------------- 离线复算(不占 GPU)
def offline_gen_row(sess, target, art, ckpt_name):
    """从**既有**采样产物离线复算生成级口径, 返回单行 dict; 失败返回 None。

    只在"产物在、指标 CSV 丢了"时兜底, 目的是**不因为丢一个 CSV 重烧一次 GPU**。
    口径刻意与子脚本同源, 不另立标准:
      - QED/SA/骨架/MW/重原子/strict: `sweep_frontier_threshold.metrics()`
        —— 与 arm_gen_metrics.py 用的是同一个真源实现;
      - 穿模: `diagnose_dataset.min_distances()` + `utils.protein_ligand.PDBProtein`,
        口袋半径 12 A / center_of_mass, 最近重原子距离 < 0.5 A 记穿模
        —— 与 select_ckpt_by_generation.evaluate_ckpt 的实体块同参数。
    """
    try:
        import numpy as np
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        if ROOT not in sys.path:
            sys.path.insert(0, ROOT)
        scripts_dir = os.path.join(ROOT, "src", "scripts")
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        from sweep_frontier_threshold import metrics          # noqa: E402
        from diagnose_dataset import min_distances            # noqa: E402
        from utils.protein_ligand import PDBProtein           # noqa: E402
    except Exception as e:  # noqa: BLE001
        log("  [警告] 离线复算依赖导入失败(%r), 只能回头重跑采样" % (e,))
        return None

    res = {"ckpt": ckpt_name, "returncode": 0, "n_finished": 0, "n_unique": 0,
           "qed_med": None, "sa_med": None, "n_scaffold": 0,
           "clash_rate": None, "mind_med": None, "mind_min": None,
           "mw_med": None, "ha_med": None, "strict_rate": None,
           "metric_source": "offline_from_session"}
    p = os.path.join(sess, "SMILES.txt")
    if not os.path.exists(p):
        log("  [警告] %s 下没有 SMILES.txt, 无法离线复算" % rel(sess))
        return None
    smis = sorted({s.strip() for s in open(p, encoding="utf-8", errors="ignore") if s.strip()})
    mols = [m for m in (Chem.MolFromSmiles(s) for s in smis) if m is not None]
    res["n_finished"] = len(smis)
    res["n_unique"] = len({Chem.MolToSmiles(m) for m in mols})
    m = metrics(smis) if smis else None
    if m:
        res["qed_med"] = round(m["qed"], 4)
        res["sa_med"] = round(m["sa"], 4)
        res["n_scaffold"] = int(m["n_scaf"])
        res["mw_med"] = round(m["mw"], 1)
        res["ha_med"] = round(m["ha"], 1)
        res["strict_rate"] = round(m["pass_pct"] / 100.0, 4)

    sdf_dir = os.path.join(sess, "SDF")
    sdfs = sorted(glob.glob(os.path.join(sdf_dir, "*.sdf"))) if os.path.isdir(sdf_dir) else []
    if sdfs:
        try:
            reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
            t = reg[target]
            prot = PDBProtein(os.path.join(ROOT, t["pdb"]))
            residues = prot.query_residues_radius(t["center"], 12.0,
                                                 criterion="center_of_mass")
            pk = os.path.join(art, "pocket.pdb")
            with open(pk, "w", encoding="ascii") as h:
                h.write(prot.residues_to_pdb_block(residues, name="POCKET"))
            ppos = np.asarray(PDBProtein(pk).pos, dtype=float)
            minds = []
            for f in sdfs[:200]:
                mol = Chem.MolFromMolFile(f, sanitize=False, removeHs=True)
                if mol is None or mol.GetNumConformers() == 0:
                    continue
                c = mol.GetConformer()
                lp = np.asarray([[c.GetAtomPosition(i).x, c.GetAtomPosition(i).y,
                                  c.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())],
                                dtype=float)
                _, mind = min_distances(lp, ppos)
                if mind is not None:
                    minds.append(mind)
            if minds:
                a = np.asarray(minds)
                res["clash_rate"] = round(float((a < CLASH_MAX).mean()), 4)
                res["mind_med"] = round(float(np.median(a)), 4)
                res["mind_min"] = round(float(a.min()), 4)
        except Exception as e:  # noqa: BLE001
            log("  [警告] 离线复算穿模失败(%r), 该列留空" % (e,))
    return res


# ---------------------------------------------------------------- 子命令构造
def build_select_cmd(ckpt, target, seed, args, keep_dir, out):
    cmd = [sys.executable, SELECT, "--ckpts", ckpt, "--target", target,
           "--num-samples", str(args.num_samples), "--beam", str(args.beam),
           "--max-steps", str(args.max_steps), "--seed", str(seed),
           "--rank-by", "both", "--keep-dir", keep_dir, "--out", out]
    if args.frontier_threshold is not None:
        cmd += ["--frontier-threshold", str(args.frontier_threshold)]
    return cmd


def build_arm_cmd(arms, baseline_name, out):
    cmd = [sys.executable, ARM, "--arms"] + ["%s=%s" % (n, p) for n, p in arms]
    if baseline_name:
        cmd += ["--baseline", baseline_name]
    return cmd + ["--out", out]


def build_ab_cmd(ckpt, target, seed, args, work, out):
    cmd = [sys.executable, AB, "--ckpts", BASELINE, ckpt, "--target", target,
           "--n", str(args.num_samples), "--beam", str(args.beam),
           "--max-steps", str(args.max_steps), "--seeds", str(seed),
           "--parallel", str(args.parallel), "--work", work, "--out", out]
    if args.frontier_threshold is not None:
        cmd += ["--frontier-threshold", str(args.frontier_threshold)]
    return cmd


def ab_names(ckpts):
    """复刻 ab_docking_compare.main() 的产物目录命名(含重名加 _1 的规则)。

    必须一致, 否则会去错目录找官方臂的采样产物。
    """
    names, seen = [], {}
    for c in ckpts:
        b = os.path.basename(c)
        if b in seen:
            seen[b] += 1
            b = "%s_%d" % (b, seen[b])
        else:
            seen[b] = 0
        names.append(b)
    return names


def find_arm_parent(ab_work, arm_name, seed, ckpt):
    """ab_docking_compare 为某臂保留产物的"检查点目录"路径(找不到返回 None)。"""
    p = os.path.join(ab_work, "%s_s%d" % (arm_name, seed), stem(ckpt))
    if os.path.isdir(p):
        return p
    cands = [d for d in glob.glob(os.path.join(ab_work, "*_s%d" % seed, "*"))
             if os.path.isdir(d) and os.path.basename(d) == stem(ckpt)]
    return cands[0] if cands else None


def pick_arm(rows, key):
    """按 ckpt/arm 列取指定臂的行: 先精确匹配, 再按 stem 匹配(arm_gen_metrics 存 stem)。"""
    if not rows or not key:
        return None
    for r in rows:
        if r.get("ckpt") == key or r.get("arm") == key:
            return r
    for r in rows:
        for col in ("ckpt", "arm"):
            v = r.get(col)
            if v and stem(v) == stem(key):
                return r
    return None


# ---------------------------------------------------------------- 三个步骤
def step_generate(ckpt, target, seed, args, base):
    """(a) 生成级: 复用已有指标 CSV -> 复用 A/B 的 select 产物 -> 离线复算 -> 才真跑 select。

    返回 (rows, rc, source, note)。**任何情况下都不重复占用 GPU**。
    """
    art = os.path.join(base, stem(ckpt))
    out = os.path.join(base, "select_gen_metrics.csv")
    ckpt_name = os.path.basename(ckpt)

    rows = read_csv(out, SELECT_COLS, "select")
    if rows:
        return rows, 0, "reused-csv", ""

    # A/B 自己生成的 select 产物就在本 (靶点,种子) 目录下, 也是真品 -> 可复用(省一次采样)
    for name in ab_names([BASELINE, ckpt])[1:]:
        cand = os.path.join(base, "ab_dock", "%s_s%d" % (name, seed), "gen_metrics.csv")
        rows = read_csv(cand, SELECT_COLS, "ab-select")
        if rows:
            log("  [复用] 生成级口径取自 A/B 的 select 产物: %s" % rel(cand))
            return rows, 0, "reused-ab-csv", ""

    sess = latest_session(art)
    if sess:
        log("  [复用] 已有采样产物 %s -> CPU 离线复算, 不再占 GPU" % rel(sess))
        row = offline_gen_row(sess, target, art, ckpt_name)
        if row:
            if not args.dry_run:
                write_rows(out, [row])
            return [row], 0, "reused-session(offline)", ""

    cmd = build_select_cmd(ckpt, target, seed, args, base, out)
    if args.dry_run:
        log("  [dry-run] " + " ".join(cmd))
        return [], 0, "dry-run", ""

    rc = run(cmd, os.path.join(base, "select.log"), "生成(select)")
    rows = read_csv(out, SELECT_COLS, "select")
    if rows:
        return rows, rc, "fresh", ""
    # 子进程可能"采样成功但 CSV 没写成/被截断" -> 从产物离线复算, 绝不白烧 GPU
    sess = latest_session(art)
    if sess:
        row = offline_gen_row(sess, target, art, ckpt_name)
        if row:
            write_rows(out, [row])
            return [row], rc, "post-run-session(offline)", "采样成功但指标CSV缺失, 已离线复算"
    return [], rc, "failed", "生成级无产物"


def step_arm_metrics(arms, baseline_name, out, log_path, tag, args):
    """(b) 尺寸/strict 口径。纯 CPU 读 SMILES, 重跑它不算"重跑产物"。"""
    rows = read_csv(out, ARM_COLS, "arm")
    if rows:
        return rows, 0, "reused"
    if args.dry_run:
        log("  [dry-run] " + " ".join(build_arm_cmd(arms, baseline_name, out)))
        return [], 0, "dry-run"
    arms = [(n, p) for n, p in arms if p and os.path.isdir(p)]
    if not arms:
        log("  [警告] arm_gen_metrics 无可用产物目录, 跳过(%s)" % tag)
        return [], 0, "no-arm"
    cmd = build_arm_cmd(arms, baseline_name, out)
    rc = run(cmd, log_path, tag)
    return read_csv(out, ARM_COLS, "arm") or [], rc, "fresh"


def step_dock(ckpt, target, seed, args, base):
    """(c) 与官方权重的同种子成对 A/B(生成 + Vina 对接)。"""
    out = os.path.join(base, "ab_docking.csv")
    rows = read_csv(out, AB_COLS, "ab")
    if rows:
        return rows, 0, "reused-csv", ""
    work = os.path.join(base, "ab_dock")
    cmd = build_ab_cmd(ckpt, target, seed, args, work, out)
    if args.dry_run:
        log("  [dry-run] " + " ".join(cmd))
        return [], 0, "dry-run", ""
    rc = run(cmd, os.path.join(base, "ab_docking.log"), "对接 A/B(ab_docking_compare)")
    rows = read_csv(out, AB_COLS, "ab")
    if not rows:
        return [], rc, "failed", "对接无产物"
    return rows, rc, "fresh", ""


# ---------------------------------------------------------------- 判定
def derive(gen, arm_row, base_arm_row, ab_row, base_ab_row, ref_ha):
    """派生对照量与退化判定, 返回 dict(并入汇总行)。

    判据与容差全部来自 ab_docking_compare.py / arm_gen_metrics.py(同源, 不另立标准)。
    """
    d = {"ha_ratio_vs_official": None, "mw_ratio_vs_official": None,
         "strict_ratio_vs_official": None, "dock_vina_med": None,
         "dock_delta_vs_official": None, "dock_ha_med": None, "dock_le_med": None,
         "dock_n_smiles": None, "dock_n_docked": None, "dock_frac_le10": None,
         "degraded": "NA", "reasons": ""}
    reasons = []

    if ab_row:
        d["dock_vina_med"] = fnum(ab_row, "vina_med")
        d["dock_ha_med"] = fnum(ab_row, "ha_med")
        d["dock_le_med"] = fnum(ab_row, "le_med")
        d["dock_n_smiles"] = inum(ab_row, "n_smiles")
        d["dock_n_docked"] = inum(ab_row, "n_docked")
        d["dock_frac_le10"] = fnum(ab_row, "frac_le10")

    # ---- 相对官方臂的尺寸/strict/生成级硬判据 ----
    if arm_row and base_arm_row:
        ha, ha_b = fnum(arm_row, "ha_med"), fnum(base_arm_row, "ha_med")
        mw, mw_b = fnum(arm_row, "mw_med"), fnum(base_arm_row, "mw_med")
        sr, sr_b = fnum(arm_row, "strict_rate"), fnum(base_arm_row, "strict_rate")
        if ha and ha_b:
            d["ha_ratio_vs_official"] = round(ha / ha_b, 3)
        if mw and mw_b:
            d["mw_ratio_vs_official"] = round(mw / mw_b, 3)
        if sr is not None and sr_b:
            d["strict_ratio_vs_official"] = round(sr / sr_b, 3)
            if d["strict_ratio_vs_official"] < STRICT_RATIO_MIN:
                reasons.append("strict通过率仅为官方%.0f%%"
                               % (100 * d["strict_ratio_vs_official"]))
        for label, ratio in (("重原子", d["ha_ratio_vs_official"]),
                             ("MW", d["mw_ratio_vs_official"])):
            if ratio is not None and abs(ratio - 1.0) > SIZE_TOL:
                log("  [警告] **尺寸差 >15%%**: %s中位 我方/官方 = %.2f (%+.0f%%) —— 原始 Vina "
                    "的差异很可能主要来自分子大小而非亲和力, 必须以 LE 复核, 单看原始分会误导。"
                    % (label, ratio, 100 * (ratio - 1)))
                reasons.append("%s差官方%+.0f%%" % (label, 100 * (ratio - 1)))
        n, n_b = fnum(arm_row, "n_finished"), fnum(base_arm_row, "n_finished")
        if n is not None and n_b is not None and (n - n_b) < TOL_N:
            reasons.append("生成数%+d(停止策略可能被破坏)" % (n - n_b))
        q, q_b = fnum(arm_row, "qed_med"), fnum(base_arm_row, "qed_med")
        if q is not None and q_b is not None and (q - q_b) < TOL_QED:
            reasons.append("QED中位%+.3f" % (q - q_b))
        s, s_b = fnum(arm_row, "sa_med"), fnum(base_arm_row, "sa_med")
        if s is not None and s_b is not None and (s - s_b) > TOL_SA:
            reasons.append("SA中位%+.2f" % (s - s_b))
        sc, sc_b = fnum(arm_row, "n_scaffold"), fnum(base_arm_row, "n_scaffold")
        if sc is not None and sc_b and sc / sc_b < TOL_SCAF:
            reasons.append("骨架数仅为官方%.0f%%" % (100 * sc / sc_b))

    # ---- 绝对口径(arm_gen_metrics 的预登记下限), 没有官方臂时仍可用 ----
    mw_now = fnum(arm_row, "mw_med") if arm_row else fnum(gen, "mw_med")
    ha_now = fnum(arm_row, "ha_med") if arm_row else fnum(gen, "ha_med")
    if mw_now is not None and mw_now < SIZE_MW_MIN:
        reasons.append("MW中位 %.0f < 预登记下限 %.0f" % (mw_now, SIZE_MW_MIN))
    if ha_now is not None and ha_now < SIZE_HA_MIN:
        reasons.append("重原子中位 %.1f < 预登记下限 %.1f" % (ha_now, SIZE_HA_MIN))
    if not (arm_row and base_arm_row) and ref_ha and ha_now:
        r = ha_now / float(ref_ha)
        if abs(r - 1.0) > SIZE_TOL:
            log("  [警告] **尺寸差 >15%%**: 重原子中位 %.1f vs 该靶点天然配体 %s (%.2f 倍) —— "
                "无官方臂对照, 此处仅作参考。" % (ha_now, ref_ha, r))
            reasons.append("重原子为该靶点天然配体的%.2f倍" % r)

    # ---- 对接级(同种子成对) ----
    if ab_row and base_ab_row:
        v, v_b = fnum(ab_row, "vina_med"), fnum(base_ab_row, "vina_med")
        if v is not None and v_b is not None:
            d["dock_delta_vs_official"] = round(v - v_b, 3)
            if d["dock_delta_vs_official"] > TOL_VINA:
                reasons.append("Vina中位%+.2f(更弱)" % d["dock_delta_vs_official"])
        for key, label in (("n_smiles", "A/B生成数"), ("n_docked", "对接数")):
            a, b = fnum(ab_row, key), fnum(base_ab_row, key)
            if a is not None and b is not None and (a - b) < TOL_N:
                reasons.append("%s%+d" % (label, a - b))
        # 尺寸不可比时, LE 是唯一还站得住的亲和力口径 -> 明确提示
        if (d["ha_ratio_vs_official"] is not None
                and abs(d["ha_ratio_vs_official"] - 1.0) > SIZE_TOL
                and d["dock_delta_vs_official"] is not None):
            le, le_b = fnum(ab_row, "le_med"), fnum(base_ab_row, "le_med")
            if le is not None and le_b is not None:
                log("  [提示] 尺寸不可比, 请以 LE 复核: LE %.3f vs 官方 %.3f (%+.3f, 越大越好)"
                    % (le, le_b, le - le_b))

    if reasons:
        d["degraded"] = "是"
    elif arm_row or base_arm_row or ab_row or gen:
        d["degraded"] = "否"
    d["reasons"] = "; ".join(reasons)
    return d


# ---------------------------------------------------------------- 汇总行
def make_row(target, seed, ckpt_name, gen, arm_row, base_arm_row, ab_row, base_ab_row,
             gen_src, multi_seed, errors=None, gen_rc=None, arm_rc=None,
             dock_rc=None, ref_ha=None):
    """拼一行汇总。缺失一律空字符串(不写 "None"), 便于 pandas/Excel 直接读。"""
    gen = gen or {}
    arm = arm_row or {}
    d = derive(gen, arm_row, base_arm_row, ab_row, base_ab_row, ref_ha)

    def pick(*vals):
        for v in vals:
            if v not in (None, ""):
                return v
        return ""

    status = "PROVISIONAL" if not multi_seed else "MULTI_SEED"
    if errors:
        status += "/PARTIAL"
    return {
        "target": target, "seed": seed, "ckpt": ckpt_name, "status": status,
        "n_finished": pick(inum(gen, "n_finished"), inum(arm, "n_finished")),
        "qed_med": pick(gen.get("qed_med"), arm.get("qed_med")),
        "sa_med": pick(gen.get("sa_med"), arm.get("sa_med")),
        "n_scaffold": pick(inum(gen, "n_scaffold"), inum(arm, "n_scaffold")),
        "clash_rate": pick(gen.get("clash_rate")),
        "mind_med": pick(gen.get("mind_med")),
        "mw_med": pick(arm.get("mw_med"), gen.get("mw_med")),
        "ha_med": pick(arm.get("ha_med"), gen.get("ha_med")),
        "strict_rate": pick(arm.get("strict_rate"), gen.get("strict_rate")),
        "dock_vina_med": pick(d["dock_vina_med"]),
        "dock_delta_vs_official": pick(d["dock_delta_vs_official"]),
        "dock_ha_med": pick(d["dock_ha_med"]),
        "dock_le_med": pick(d["dock_le_med"]),
        "dock_n_smiles": pick(d["dock_n_smiles"]),
        "dock_n_docked": pick(d["dock_n_docked"]),
        "dock_frac_le10": pick(d["dock_frac_le10"]),
        "ha_ratio_vs_official": pick(d["ha_ratio_vs_official"]),
        "mw_ratio_vs_official": pick(d["mw_ratio_vs_official"]),
        "strict_ratio_vs_official": pick(d["strict_ratio_vs_official"]),
        "degraded": d["degraded"], "reasons": d["reasons"],
        "gen_source": gen_src,
        "gen_rc": "" if gen_rc is None else gen_rc,
        "arm_rc": "" if arm_rc is None else arm_rc,
        "dock_rc": "" if dock_rc is None else dock_rc,
        "errors": "; ".join(errors) if errors else "",
    }


def check_meta(base, target, ckpt_name, meta_gen, gen_src, dry_run):
    """核对复用产物的采样参数; 返回需写进 errors 列的提示(可能为空)。"""
    if dry_run:
        return ""
    old = load_meta(base)
    if old is None:
        save_meta(base, {"target": target, "ckpt": ckpt_name, "gen": meta_gen})
        return ("NO_META(该目录无 eval_meta.json, 复用产物的采样参数无法核对)"
                if gen_src.startswith("reused") else "")
    prev = old.get("ckpt")
    if prev and prev != ckpt_name:
        log("  [警告] 该目录上次评测的是 %s, 本次是 %s -> 产物可能不属于本次检查点"
            % (prev, ckpt_name))
    if gen_src.startswith("reused") and old.get("gen") and old["gen"] != meta_gen:
        log("  [警告] 复用产物的采样参数与本次不一致: 产物=%s 本次=%s (不重跑, 但结论需注明口径)"
            % (old["gen"], meta_gen))
        return "PARAM_MISMATCH(产物 %s != 本次 %s)" % (old["gen"], meta_gen)
    return ""


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser(
        description="一次评估一个检查点在多个靶点上的表现(编排 生成/尺寸/对接 三个子脚本)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="安全: (靶点,种子) 严格串行, 不并发采样(避免与训练抢显存); 子进程非 0 只记录不中断; "
               "已有产物一律复用, 不重跑; 单种子结果标 PROVISIONAL; 尺寸差 >15% 打警告; "
               "本脚本不调用 git。\n"
               "相对路径按仓库根解析。a/b/c 三步分别调用 src/scripts/ 下的 "
               "select_ckpt_by_generation.py / arm_gen_metrics.py / ab_docking_compare.py。")
    ap.add_argument("--ckpt", required=True, help="被评测的检查点路径(.pt)")
    ap.add_argument("--targets", nargs="+", default=["A2A", "B2AR", "D3", "5HT2B"],
                    help="靶点名(须在 configs/targets.json 里), 默认 A2A B2AR D3 5HT2B")
    ap.add_argument("--num-samples", type=int, default=50,
                    help="每个 (靶点,种子) 的采样分子数, 默认 50")
    ap.add_argument("--beam", type=int, default=50, help="束宽 beam_size, 默认 50")
    ap.add_argument("--max-steps", type=int, default=40, help="最大生成步数 max_steps, 默认 40")
    ap.add_argument("--seeds", nargs="+", type=int, default=[2024, 2025],
                    help="随机种子(>=2 个才不标 PROVISIONAL), 默认 2024 2025")
    ap.add_argument("--parallel", type=int, default=6,
                    help="**只透传给 A/B 的 Vina 对接**分块并发数, 默认 6; 采样永远串行")
    ap.add_argument("--frontier-threshold", type=float, default=None,
                    help="frontier 判定阈值(直接控制分子大小), 默认 None = 不干预, 由子脚本"
                         "使用各自默认 0.0; 若该权重自带推荐阈值会打印口径提示")
    ap.add_argument("--keep-dir", default=os.path.join(OUTPUTS, "eval_ckpt"),
                    help="产物根目录, 默认 outputs/eval_ckpt; 每个 (靶点,种子) 一个子目录, "
                         "已有产物直接复用")
    ap.add_argument("--out", default=os.path.join(OUTPUTS, "eval_ckpt_summary.csv"),
                    help="汇总 CSV, 默认 outputs/eval_ckpt_summary.csv")
    ap.add_argument("--skip-dock", action="store_true",
                    help="只做生成级(不跑 ab_docking_compare): 采样次数 3 -> 1, 汇总无对接列")
    ap.add_argument("--dry-run", action="store_true",
                    help="额外开关: 只打印将要执行的子命令, 不执行子进程、不建目录、不写文件")
    args = ap.parse_args()

    try:  # 中文输出被重定向到 cp936 文件时不至于 UnicodeEncodeError
        sys.stdout.reconfigure(errors="replace")
    except Exception:  # noqa: BLE001
        pass

    args.ckpt = resolve_path(args.ckpt)
    args.keep_dir = root_join(args.keep_dir)
    args.out = root_join(args.out)

    log("=" * 96)
    log("检查点多靶点评测 | ckpt=%s" % rel(args.ckpt))
    log("靶点=%s | 种子=%s | 规模 %d/%d/%d | frontier_threshold=%s | 对接=%s"
        % (args.targets, args.seeds, args.num_samples, args.beam, args.max_steps,
           "不干预(子脚本默认 0.0)" if args.frontier_threshold is None
           else args.frontier_threshold,
           "跳过(--skip-dock)" if args.skip_dock else "开(与官方权重成对)"))
    log("产物根 %s | 汇总 %s" % (rel(args.keep_dir), rel(args.out)))
    if args.dry_run:
        log("** --dry-run: 只打印子命令, 不执行、不建目录、不写文件 **")
    log("=" * 96)

    if not os.path.exists(args.ckpt):
        raise SystemExit("检查点不存在: %s" % args.ckpt)
    ckpt_name = os.path.basename(args.ckpt)

    # 口径提示(只读侧车表, 不做 torch.load: 评测机上 GPU 可能在训练, 不必做权重读放大)
    if args.frontier_threshold is None:
        try:
            from paths import recommended_frontier_threshold as _rec
            v, src = _rec(args.ckpt, peek_ckpt=False)
            if abs(v) > 1e-9:
                log("[口径提示] %s 的推荐 frontier_threshold 是 %+.2f(来源 %s); 本次不干预, "
                    "子脚本会用默认 0.0。要对齐出货口径请显式传 --frontier-threshold %+.2f"
                    % (ckpt_name, v, src, v))
        except Exception:  # noqa: BLE001
            pass

    dock_on = not args.skip_dock
    if dock_on and not os.path.exists(BASELINE):
        log("[警告] 找不到官方权重 %s -> 无法做同种子成对 A/B, 本轮退化为只做生成级"
            % rel(BASELINE))
        dock_on = False
    if not dock_on and args.frontier_threshold is not None:
        log("[提示] frontier_threshold=%s 会写进生成配置(影响分子大小); 但本轮无官方臂对照, "
            "尺寸只能与该靶点天然配体比。" % args.frontier_threshold)

    try:
        reg = json.load(open(os.path.join(CONFIGS, "targets.json"), encoding="utf-8"))
    except Exception as e:  # noqa: BLE001
        raise SystemExit("读不了 configs/targets.json: %r" % (e,))
    unknown = [t for t in args.targets if t not in reg]
    if unknown:
        log("[警告] 未知靶点 %s(不在 configs/targets.json), 各记一行错误后跳过" % unknown)

    n_pairs = len([t for t in args.targets if t in reg]) * len(args.seeds)
    log("计划: %d 个 (靶点,种子), 每对最多 %s 次采样 + 1 轮对接; 采样严格串行"
        % (n_pairs, "3(本次 + A/B 两臂)" if dock_on else "1(仅本次)"))
    log("-" * 96)

    multi_seed = len(args.seeds) >= 2
    if not multi_seed:
        log("[警告] 只有 %d 个种子 -> 全部结果标 PROVISIONAL(未证实), 不得据此晋级。"
            % len(args.seeds))
        log("        依据: 种子间 Vina 中位噪声实测约 0.7 kcal/mol, 远大于同臂重复的 ~0.14。")
        log("-" * 96)

    meta_gen = dict(num_samples=args.num_samples, beam=args.beam,
                    max_steps=args.max_steps, frontier_threshold=args.frontier_threshold)
    names = ab_names([BASELINE, args.ckpt])
    if names[0] == names[1]:
        log("[警告] --ckpt 与官方权重是同一路径 -> A/B 两臂相同, Δ 必然接近 0, 该轮无判据意义。")

    rows = []
    for target in args.targets:
        for seed in args.seeds:
            base = os.path.join(args.keep_dir, target, str(seed))
            log("[%s | 种子 %d]" % (target, seed))
            if not args.dry_run:      # dry-run 不建目录
                ensure_dir(base)

            if target in unknown:
                rows.append(make_row(target, seed, ckpt_name, None, None, None, None, None,
                                     "unknown-target", multi_seed,
                                     errors=["靶点不在 configs/targets.json"]))
                continue

            # ---- a) 生成级 ----
            gen_rows, gen_rc, gen_src, gen_note = step_generate(args.ckpt, target, seed,
                                                               args, base)
            if gen_note:
                log("  [提示] %s" % gen_note)
            errors = [gen_note] if gen_note else []
            if gen_rc != 0:
                errors.append("select rc=%d" % gen_rc)
            gen = gen_rows[0] if gen_rows else None
            if gen is None and not args.dry_run:
                log("  [失败] 生成级无产物, 记一行后继续下一个 (靶点,种子)")
                rows.append(make_row(target, seed, ckpt_name, None, None, None, None, None,
                                     gen_src, multi_seed, errors=errors, gen_rc=gen_rc))
                write_summary(rows, args.out)
                log("-" * 96)
                continue

            note = check_meta(base, target, ckpt_name, meta_gen, gen_src, args.dry_run)
            if note:
                errors.append(note)

            # ---- b) 尺寸/strict ----
            art = os.path.join(base, stem(args.ckpt))
            arm_rows, arm_rc, _ = step_arm_metrics(
                [(stem(args.ckpt), art)], None, os.path.join(base, "arm_gen_metrics.csv"),
                os.path.join(base, "arm_gen_metrics.log"), "尺寸/strict(arm_gen_metrics)", args)
            if arm_rc != 0:
                errors.append("arm_gen_metrics rc=%d" % arm_rc)
            arm_row = pick_arm(arm_rows, stem(args.ckpt)) or pick_arm(arm_rows, ckpt_name)

            # ---- c) 与官方权重的同种子成对 A/B ----
            ab_rows, ab_row, base_ab_row, base_arm_row, dock_rc = [], None, None, None, None
            if dock_on:
                ab_rows, dock_rc, _src, dock_note = step_dock(args.ckpt, target, seed, args, base)
                if dock_note:
                    log("  [提示] %s" % dock_note)
                    errors.append(dock_note)
                if dock_rc:
                    errors.append("ab_docking_compare rc=%d" % dock_rc)
                ab_row = pick_arm(ab_rows, names[1])
                base_ab_row = pick_arm(ab_rows, names[0])
                if base_ab_row is None and ab_rows:
                    errors.append("A/B CSV 里找不到官方臂 %s" % names[0])
                # 官方臂产物到手后补一次尺寸/strict 对照(纯 CPU, 不占 GPU)
                if base_ab_row is not None or args.dry_run:
                    p_official = find_arm_parent(os.path.join(base, "ab_dock"),
                                                 names[0], seed, BASELINE)
                    if p_official is None and args.dry_run:
                        # dry-run: 产物目录还不存在, 用预期的路径把子命令打出来
                        p_official = os.path.join(base, "ab_dock",
                                                  "%s_s%d" % (names[0], seed), stem(BASELINE))
                    if p_official:
                        vs_rows, vs_rc, _ = step_arm_metrics(
                            [("official", p_official), (stem(args.ckpt), art)], "official",
                            os.path.join(base, "arm_gen_metrics_vs_official.csv"),
                            os.path.join(base, "arm_gen_metrics.log"),
                            "对照官方臂(arm_gen_metrics)", args)
                        if vs_rc != 0:
                            errors.append("arm_gen_metrics(vs官方) rc=%d" % vs_rc)
                        base_arm_row = pick_arm(vs_rows, "official")
                        if arm_row is None:
                            arm_row = pick_arm(vs_rows, stem(args.ckpt))
                    else:
                        log("  [提示] 找不到官方臂的采样产物目录 -> 本轮无尺寸/strict 对照")
                        errors.append("未找到官方臂产物目录")

            ref_ha = (reg.get(target) or {}).get("n_ligand_heavy")
            row = make_row(target, seed, ckpt_name, gen, arm_row, base_arm_row, ab_row,
                           base_ab_row, gen_src, multi_seed, errors=errors,
                           gen_rc=gen_rc, arm_rc=arm_rc, dock_rc=dock_rc, ref_ha=ref_ha)
            rows.append(row)
            if not args.dry_run:
                write_summary(rows, args.out)
            log("        完成 %s | QED %s | SA %s | 骨架 %s | 穿模 %s | MW %s | 重原子 %s | "
                "strict %s | Vina %s (Δ官方 %s) | LE %s | 退化 %s%s"
                % (row["n_finished"], row["qed_med"], row["sa_med"], row["n_scaffold"],
                   row["clash_rate"], row["mw_med"], row["ha_med"], row["strict_rate"],
                   row["dock_vina_med"], row["dock_delta_vs_official"], row["dock_le_med"],
                   row["degraded"], ("  [" + str(row["reasons"]) + "]") if row["reasons"] else ""))
            log("-" * 96)

    # ---- 收尾 ----
    log("=" * 96)
    if args.dry_run:
        log("--dry-run 结束: 未执行任何子进程, 未建目录, 未写文件。")
        return 0
    usable = [r for r in rows if r.get("n_finished") not in (None, "")]
    bad = [r for r in rows if r.get("errors")]
    if not dock_on:
        log("[提示] 本轮只做生成级(--skip-dock 或官方权重缺失): 汇总无对接列, "
            "不得据此对结合强度下任何结论。")
    log("汇总: %d 行(可用 %d) | 有告警/失败 %d 行 | 对接 %s | 单种子=%s"
        % (len(rows), len(usable), len(bad),
           "开" if dock_on else "跳过", "否" if multi_seed else "是"))
    if not multi_seed and rows:
        log("** 单种子 -> 全部行 status=PROVISIONAL, 不得用来宣布改进 **")
    deg = [r for r in rows if r.get("degraded") == "是"]
    if deg:
        log("** 相对官方退化 %d 行 **:" % len(deg))
        for r in deg:
            log("   %-6s 种子 %-5s %s" % (r["target"], r["seed"], r["reasons"]))
    if bad:
        log("需要人看的行(errors 列已在 CSV 里留痕):")
        for r in bad:
            log("   %-6s 种子 %-5s %s" % (r["target"], r["seed"], r["errors"]))
    log("列口径: ha_med=生成口径重原子中位, dock_ha_med=对接口径重原子中位, "
        "dock_le_med=-Vina/重原子(越大越好)")
    log("明细: %s (utf-8-sig)" % rel(args.out))
    log("=" * 96)
    return 0 if usable else 1


if __name__ == "__main__":
    sys.exit(main())
