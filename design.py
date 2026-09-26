# -*- coding: utf-8 -*-
"""design.py —— 一键分子设计/生成入口。

在指定 GPCR 靶点的结合口袋内, 用 Pocket2Mol 的等变图神经网络 + 口袋引导束搜索
生成候选分子, 输出 3D 结构(SDF)与 SMILES。

用法示例
--------
    # 用登记表里的靶点(推荐)
    python design.py --target A2A

    # 指定规模与引导强度
    python design.py --target D3 --num-samples 100 --beam 100 --lam 3.0 --seed 2024

    # 用自己的受体结构 + 口袋中心(埃)
    python design.py --pdb data/targets/4EIY_A2A受体.pdb --center "-0.42,8.53,17.13"

输出
----
    <outdir>/<session>/SDF/*.sdf     每个候选的三维结构
    <outdir>/<session>/SMILES.txt    候选 SMILES 清单
    <outdir>/<session>/samples_*.pt  采样过程快照
    <outdir>/<session>/design_manifest.json  本次运行的参数/种子/版本(可溯源)
"""
import argparse
import json
import os
import platform
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from paths import (ROOT, CONFIGS, OUTPUTS, PRETRAINED, DEFAULT_CKPT, TARGETS,  # noqa: E402
                   ensure_dir, missing_hint)

DEFAULT_TEMPLATE = os.path.join(CONFIGS, "sample_for_pdb_guided_l3.yml")


def load_targets():
    p = os.path.join(CONFIGS, "targets.json")
    if not os.path.exists(p):
        raise SystemExit(
            "[design] 缺少靶点登记表 %s\n请先运行: python src/scripts/build_target_registry.py" % p)
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def resolve_target(args, reg):
    """返回 (pdb_path, center_list, target_name, reference_smiles)。"""
    if args.target:
        if args.target not in reg:
            raise SystemExit("[design] 未知靶点 %s; 可选: %s"
                             % (args.target, ", ".join(sorted(reg))))
        t = reg[args.target]
        pdb = t["pdb"] if os.path.isabs(t["pdb"]) else os.path.join(ROOT, t["pdb"])
        return pdb, list(t["center"]), args.target, t.get("ligand_smiles", "")
    if not args.pdb:
        raise SystemExit("[design] 需要 --target, 或同时给出 --pdb 与 --center")
    if not args.center:
        raise SystemExit("[design] 给出 --pdb 时必须同时给出 --center \"x,y,z\"")
    if not os.path.exists(args.pdb):
        raise SystemExit("[design] 受体文件不存在: %s" % args.pdb)
    return args.pdb, args.center, args.name or "custom", ""


def build_config(args, center, out_cfg):
    """在模板配置基础上生成本次运行配置(不改动仓库里的模板)。"""
    import yaml
    with open(args.config, encoding="utf-8-sig") as f:
        cfg = yaml.safe_load(f)
    cfg.setdefault("model", {})
    cfg["model"]["checkpoint"] = args.ckpt or DEFAULT_CKPT
    s = cfg.setdefault("sample", {})
    s["seed"] = args.seed
    s["num_samples"] = args.num_samples
    s["beam_size"] = args.beam
    s["max_steps"] = args.max_steps
    g = s.setdefault("guided", {})
    g["enabled"] = bool(args.lam and args.lam > 0)
    g["qed_w"] = args.qed_w
    g["sa_w"] = args.sa_w
    g["lam"] = args.lam
    g["diversity_w"] = args.diversity_w
    s["relax_output"] = bool(args.relax)
    ensure_dir(os.path.dirname(out_cfg))
    with open(out_cfg, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)
    return cfg


def main():
    ap = argparse.ArgumentParser(
        description="一键分子设计/生成 (Pocket2Mol + 口袋引导束搜索)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--target", help="靶点名, 取自 configs/targets.json (如 A2A/B2AR/D3/5HT2B/A1AR)")
    ap.add_argument("--pdb", help="自定义受体 PDB 路径(与 --center 同时使用)")
    ap.add_argument("--center", type=lambda s: [float(x) for x in s.split(",")],
                    help="口袋中心 \"x,y,z\" (埃)")
    ap.add_argument("--name", help="自定义靶点的显示名")
    ap.add_argument("--num-samples", type=int, default=50, help="目标生成分子数")
    ap.add_argument("--beam", type=int, default=100, help="束宽")
    ap.add_argument("--max-steps", type=int, default=50, help="步数上限(安全上限, 非目标值)")
    ap.add_argument("--seed", type=int, default=2024, help="随机种子(固定以保证可复现)")
    ap.add_argument("--lam", type=float, default=3.0, help="引导强度 λ; 0 表示关闭引导")
    ap.add_argument("--qed-w", type=float, default=1.0, help="QED 权重")
    ap.add_argument("--sa-w", type=float, default=1.0, help="SA 权重")
    ap.add_argument("--diversity-w", type=float, default=0.5, help="多样性惩罚权重")
    ap.add_argument("--relax", type=int, default=0, choices=[0, 1],
                    help="是否对输出做构象精修(保持姿态)")
    ap.add_argument("--config", default=DEFAULT_TEMPLATE, help="采样配置模板")
    ap.add_argument("--ckpt", default=None,
                    help="自定义模型权重路径; 默认用本项目自训并通过两轮 A/B 晋级的 "
                         "models/7-eonmol_ft_gpcr_v2.pt(不存在时回退官方权重)。"
                         "想跑官方基线请显式传 models/pretrained_Pocket2Mol.pt")
    ap.add_argument("--outdir", default=os.path.join(OUTPUTS, "design"), help="输出根目录")
    ap.add_argument("--device", default="cuda", help="cuda 或 cpu")
    args = ap.parse_args()

    reg = load_targets()
    pdb, center, tname, ref_smi = resolve_target(args, reg)
    ensure_dir(args.outdir)
    ts = time.strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(args.outdir, "%s_%s" % (tname, ts))
    ensure_dir(run_dir)
    cfg_path = os.path.join(run_dir, "design_config.yml")

    miss = missing_hint([("预训练权重", PRETRAINED), ("采样配置模板", args.config),
                         ("受体结构", pdb)])
    if miss:
        raise SystemExit("[design] 缺少必需文件:\n" + miss +
                         "\n请按 README 的“模型权重”一节准备。")

    build_config(args, center, cfg_path)

    print("=" * 74)
    print("7-eonmol design  |  靶点=%s  受体=%s" % (tname, os.path.basename(pdb)))
    print("  口袋中心 = %s" % ", ".join("%.2f" % v for v in center))
    print("  规模 = %d 分子 / 束宽 %d / 步数上限 %d / 种子 %d" % (
        args.num_samples, args.beam, args.max_steps, args.seed))
    print("  引导 = %s (λ=%.2f, 多样性w=%.2f, 构象精修=%s)" % (
        "开" if args.lam > 0 else "关", args.lam, args.diversity_w, bool(args.relax)))
    print("  输出 = %s" % run_dir)
    print("=" * 74, flush=True)

    # 用 --center=x,y,z 的等号形式传递: 若写成 "--center" " -0.4,8.5,17.1" 这种空格
    # 分隔形式, argparse 会把以 "-" 开头的值误判为选项(这正是旧文档里"前导空格不能省"
    # 的原因)。等号形式对负数首值同样安全, 无需前导空格。
    cmd = [sys.executable, os.path.join(ROOT, "src", "sample_for_pdb.py"),
           "--pdb_path", pdb,
           "--center=" + ",".join("%.2f" % v for v in center),
           "--config", cfg_path,
           "--device", args.device,
           "--outdir", run_dir]
    t0 = time.time()
    rc = subprocess.run(cmd, cwd=ROOT).returncode
    el = time.time() - t0

    # sample_for_pdb 会在 run_dir 下自建 <配置名>_<时间戳> 会话目录
    sessions = sorted(d for d in (os.path.join(run_dir, x) for x in os.listdir(run_dir))
                      if os.path.isdir(d))
    sess = sessions[-1] if sessions else run_dir
    n_sdf = n_smi = 0
    if sess:
        sdf_dir = os.path.join(sess, "SDF")
        if os.path.isdir(sdf_dir):
            n_sdf = len([x for x in os.listdir(sdf_dir) if x.endswith(".sdf")])
        smi = os.path.join(sess, "SMILES.txt")
        if os.path.exists(smi):
            n_smi = len([x for x in open(smi, encoding="utf-8", errors="ignore")
                         if x.strip()])

    import torch
    manifest = {
        "stage": "design", "target": tname, "receptor_pdb": os.path.relpath(pdb, ROOT),
        "center": center, "cocrystal_ligand_smiles": ref_smi,
        "params": {"num_samples": args.num_samples, "beam_size": args.beam,
                   "max_steps": args.max_steps, "seed": args.seed, "lam": args.lam,
                   "qed_w": args.qed_w, "sa_w": args.sa_w,
                   "diversity_w": args.diversity_w, "relax_output": bool(args.relax)},
        "result": {"session_dir": os.path.relpath(sess, ROOT) if sess else None,
                   "n_sdf": n_sdf, "n_smiles": n_smi, "returncode": rc,
                   "elapsed_sec": round(el, 1)},
        "env": {"python": sys.version.split()[0], "torch": torch.__version__,
                "cuda": torch.cuda.is_available(),
                "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
                "platform": platform.platform()},
        "timestamp": ts,
    }
    with open(os.path.join(run_dir, "design_manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=1)

    print("=" * 74)
    print("完成: SDF=%d, SMILES=%d, 用时 %.1f 秒, 返回码 %d" % (n_sdf, n_smi, el, rc))
    if sess:
        print("会话目录: %s" % sess)
        print("下一步(打分排序并生成候选清单):")
        print("    python predict.py --skip-design --runs \"%s\" --target %s" % (sess, tname))
    print("=" * 74)
    return 0 if rc == 0 else rc


if __name__ == "__main__":
    sys.exit(main())
