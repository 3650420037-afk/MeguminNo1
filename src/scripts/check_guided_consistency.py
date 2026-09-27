# -*- coding: utf-8 -*-
"""护栏：保证「引导参数」在模板、入口脚本、评测脚本三处完全一致。

为什么需要它
------------
历史缺陷（2026-09 发现）：`diversity_w` 在 `src/sample_for_pdb.py` 里的**代码默认是 0.0**
（键缺失即关闭），而 `design.py` / `predict.py` 会**注入** 0.5，模板
`configs/sample_for_pdb_guided_l3.yml` 当时**没有这个键**。后果有两层：

1. 「直接跑模板」与「走入口脚本」行为不同 —— 文档宣称的多样性惩罚会静默失效；
2. 更严重：检查点筛选 `src/scripts/select_ckpt_by_generation.py` 沿用模板，
   于是 `outputs/ab_*` 五组 A/B **全部在惩罚=0 下评测**，而出货药库在 0.5 下生成，
   **评测目标 != 出货目标**。

本脚本把那三处的取值读出来直接比对，任何一处漂移都会 FAIL。
CI/夜间流程应把本脚本作为前置门禁。

用法
----
    python src/scripts/check_guided_consistency.py
    退出码 0 = 一致；1 = 不一致（并打印每处的实际取值）
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, CONFIGS  # noqa: E402

TEMPLATE = os.path.join(CONFIGS, "sample_for_pdb_guided_l3.yml")
GUIDED_KEYS = ["enabled", "qed_w", "sa_w", "lam", "diversity_w"]


def _read(p):
    return open(p, encoding="utf-8-sig", errors="ignore").read()


def template_guided():
    """从模板解析 guided 块（只用标准库，避免依赖 pyyaml 的版本差异）。

    注意：不能用 `(?=^\\s*\\w+:)` 这种「遇到下一个键就结束」的写法 —— guided 的
    **第一个子键**本身就是一个 `^\\s*\\w+:`，会让块体长度为 0。正确判据是
    「缩进回到 <= guided 所在行的缩进」才算离开该块。
    """
    out = {}
    in_guided = False
    base = None
    for line in _read(TEMPLATE).splitlines():
        if not line.strip() or line.strip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        if re.match(r"\s*guided\s*:\s*$", line):
            in_guided, base = True, indent
            continue
        if in_guided:
            if indent <= base:
                in_guided = False
                continue
            mm = re.match(r"\s*([A-Za-z_]\w*)\s*:\s*(\S+)", line)
            if mm:
                out[mm.group(1)] = mm.group(2)
    return out


def py_default(path, key):
    """读某脚本里 `--<key>` 参数的 default=... 或模块级常量。"""
    txt = _read(path)
    m = re.search(r'"--%s"[^)]*?default=([0-9.]+)' % key.replace("_", "-"), txt, re.S)
    if m:
        return m.group(1)
    m = re.search(r"^DIVERSITY_W\s*=\s*([0-9.]+)", txt, re.M)
    return m.group(1) if m else None


def code_default(key):
    """读 sample_for_pdb.py 里 guided_cfg.get(key, X) 的 X —— 即键缺失时的兜底值。"""
    txt = _read(os.path.join(ROOT, "src", "sample_for_pdb.py"))
    m = re.search(r"guided_cfg\.get\(\s*'%s'\s*,\s*([0-9.]+)" % key, txt)
    return m.group(1) if m else None


def main():
    tpl = template_guided()
    design = os.path.join(ROOT, "design.py")
    predict = os.path.join(ROOT, "predict.py")
    selck = os.path.join(ROOT, "src", "scripts", "select_ckpt_by_generation.py")

    rows, bad = [], []

    # 1) 模板必须显式写出每一个引导键
    for k in GUIDED_KEYS:
        if k not in tpl:
            bad.append("模板 %s 缺少 guided.%s" % (os.path.basename(TEMPLATE), k))
        rows.append(("模板", "guided.%s" % k, tpl.get(k, "<缺失>")))

    # 2) design.py / predict.py 的 --diversity-w 默认值必须等于模板值
    tpl_div = tpl.get("diversity_w")
    for name, p in (("design.py", design), ("predict.py", predict)):
        v = py_default(p, "diversity_w")
        rows.append((name, "--diversity-w default", v))
        if v is None:
            bad.append("%s 找不到 --diversity-w 默认值" % name)
        elif tpl_div is not None and float(v) != float(tpl_div):
            bad.append("%s 默认 %s != 模板 %s" % (name, v, tpl_div))

    # 3) 评测脚本的 DIVERSITY_W 必须等于模板值
    v = py_default(selck, "diversity-w")
    rows.append(("select_ckpt_by_generation.py", "DIVERSITY_W", v))
    if v is None:
        bad.append("select_ckpt_by_generation.py 没有 DIVERSITY_W 常量（评测目标会随模板漂移）")
    elif tpl_div is not None and float(v) != float(tpl_div):
        bad.append("select_ckpt_by_generation.py DIVERSITY_W=%s != 模板 %s" % (v, tpl_div))

    # 4) 评测脚本必须显式覆盖 diversity_w，而不能沿用模板
    if not re.search(r'g\["diversity_w"\]\s*=', _read(selck)):
        bad.append("select_ckpt_by_generation.py 未显式注入 diversity_w（会再次沿用模板）")

    # 5) 代码兜底值应当是 0.0 并且模板已显式给出 —— 显式优于隐式
    cd = code_default("diversity_w")
    rows.append(("src/sample_for_pdb.py", "get('diversity_w', X) 兜底", cd))

    w = max(len(a) for a, _, _ in rows) if rows else 10
    print("=" * 62)
    print("引导参数一致性检查")
    print("=" * 62)
    for a, b, c in rows:
        print("  %-*s | %-28s = %s" % (w, a, b, c))
    print("-" * 62)
    if bad:
        print("结果: FAIL (%d 项)" % len(bad))
        for b in bad:
            print("  ✗ %s" % b)
        return 1
    print("结果: PASS —— 模板 / design.py / predict.py / 评测脚本 四处取值一致")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
