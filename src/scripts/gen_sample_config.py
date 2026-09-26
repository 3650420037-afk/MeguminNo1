# -*- coding: utf-8 -*-
"""采样配置生成器 (根治"文本末尾追加"导致的 YAML 非法问题)

背景: 旧实现(GUI app_easy.py、overnight_*.ps1)在模板文件**末尾追加** `    diversity_w: x`
缩进行。当模板最后一个 `guided:` 块之后又新增了其它同级键(如 `relax_output`)时,
追加会产生非法 YAML -> load_config 抛 ScannerError -> 采样/训练秒崩。
本脚本用 PyYAML 做结构化读改写, 不再依赖模板的文本布局, 并顺带支持 lam 覆盖。

用法:
  python scripts/gen_sample_config.py --template configs/sample_for_pdb_guided_l3.yml \
      --out <输出yml> [--seed N] [--num-samples N] [--beam N] [--max-steps N] \
      [--lam F] [--diversity-w F] [--guided 0|1] [--relax-output 0|1]

约定: 未提供的参数保持模板原值; 写出的 yml 为 UTF-8 无 BOM。
退出码: 0 成功 / 2 参数或模板错误 (调用方应检查)。
"""
import argparse
import sys

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML 不可用", file=sys.stderr)
    sys.exit(2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--num-samples", type=int)
    ap.add_argument("--beam", type=int)
    ap.add_argument("--max-steps", type=int)
    ap.add_argument("--lam", type=float)
    ap.add_argument("--diversity-w", type=float)
    ap.add_argument("--guided", type=int, choices=[0, 1])
    ap.add_argument("--relax-output", type=int, choices=[0, 1])
    args = ap.parse_args()

    try:
        with open(args.template, encoding="utf-8-sig") as f:
            cfg = yaml.safe_load(f) or {}
    except Exception as e:
        print("ERROR: 读取模板失败: %s" % e, file=sys.stderr)
        sys.exit(2)

    sample = cfg.setdefault("sample", {})
    guided = sample.setdefault("guided", {})

    if args.seed is not None:
        sample["seed"] = args.seed
    if args.num_samples is not None:
        sample["num_samples"] = args.num_samples
    if args.beam is not None:
        sample["beam_size"] = args.beam
    if args.max_steps is not None:
        sample["max_steps"] = args.max_steps

    if args.guided is not None:
        guided["enabled"] = bool(args.guided)
    if args.lam is not None:
        guided["lam"] = args.lam
    if args.diversity_w is not None:
        guided["diversity_w"] = args.diversity_w
    if args.relax_output is not None:
        sample["relax_output"] = bool(args.relax_output)

    try:
        with open(args.out, "w", encoding="utf-8") as f:   # 无 BOM
            yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False, default_flow_style=False)
    except Exception as e:
        print("ERROR: 写出配置失败: %s" % e, file=sys.stderr)
        sys.exit(2)

    # 自检: 立即回读解析, 确保产物一定可被 load_config 使用
    try:
        with open(args.out, encoding="utf-8-sig") as f:
            yaml.safe_load(f)
    except Exception as e:
        print("ERROR: 生成结果自检不通过: %s" % e, file=sys.stderr)
        sys.exit(2)
    print("OK %s" % args.out)


if __name__ == "__main__":
    main()
