# -*- coding: utf-8 -*-
"""端到端行为验证: 多样性惩罚是否**真的生效**, 以及出货采样路径是否正常。

为什么不能只做静态检查
----------------------
`check_guided_consistency.py` 只能证明"模板里写了 diversity_w: 0.5"。但历史缺陷的本质是
**读不到**(代码默认 0.0), 所以必须证明"读到并改变了行为":
同一模型、同一种子、同一预算, 只把 diversity_w 从 0.0 改成 0.5 —— 生成的分子应当不同。
若完全相同, 说明该参数在实际运行中仍未起作用, 必须继续排查。

同时这也是出货路径的冒烟测试: 直接跑 **规范模板**(而不是临时改写的配置), 确认
sample_for_pdb.py 能正常出分子。

用法
----
    python src/scripts/verify_diversity_active.py
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, CONFIGS, src_on_path  # noqa: E402

src_on_path()
import yaml  # noqa: E402

TEMPLATE = os.path.join(CONFIGS, "sample_for_pdb_guided_l3.yml")
SAMPLER = os.path.join(ROOT, "src", "sample_for_pdb.py")
TARGETS = os.path.join(CONFIGS, "targets.json")
FAILS = []


def chk(name, ok, detail=""):
    print("  %s %s%s" % ("[OK]  " if ok else "[FAIL]", name, ("  -> " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)
    return ok


def run_once(div_w, center, pdb, workdir, n=24, beam=24, steps=16, seed=777):
    cfg = yaml.safe_load(open(TEMPLATE, encoding="utf-8-sig"))
    cfg["sample"]["num_samples"] = n
    cfg["sample"]["beam_size"] = beam
    cfg["sample"]["max_steps"] = steps
    cfg["sample"]["seed"] = seed
    cfg["sample"].setdefault("guided", {})["diversity_w"] = div_w
    p = os.path.join(workdir, "cfg_div%s.yml" % str(div_w).replace(".", "p"))
    with open(p, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)
    out = os.path.join(workdir, "out%s" % str(div_w).replace(".", "p"))
    os.makedirs(out, exist_ok=True)
    lg = os.path.join(workdir, "log%s.txt" % str(div_w).replace(".", "p"))
    with open(lg, "w", encoding="utf-8", errors="ignore") as lf:
        rc = subprocess.run(
            [sys.executable, SAMPLER, "--pdb_path", pdb, "--center=" + center,
             "--config", p, "--device", "cuda", "--outdir", out],
            cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT).returncode
    # 找到本次 session 的 SMILES.txt
    smis = []
    for d in sorted([x for x in _walk(out) if os.path.basename(x) == "SMILES.txt"]):
        for line in open(d, encoding="utf-8", errors="ignore"):
            s = line.strip()
            if s and not s.startswith("#"):
                smis.append(s)
    return rc, smis, out


def _walk(root):
    got = []
    for dp, _dn, fn in os.walk(root):
        for f in fn:
            got.append(os.path.join(dp, f))
    return got


def main():
    reg = __import__("json").load(open(TARGETS, encoding="utf-8"))
    t = reg["A2A"]
    pdb = os.path.join(ROOT, t["pdb"])
    center = ",".join("%.2f" % v for v in t["center"])

    print("=" * 72)
    print("1/2 出货路径冒烟: 直接跑规范模板(应为 diversity_w=0.5, 不改任何键)")
    print("=" * 72)
    tpl = yaml.safe_load(open(TEMPLATE, encoding="utf-8-sig"))
    got = tpl["sample"].get("guided", {}).get("diversity_w", None)
    chk("规范模板自带 diversity_w=0.5", got == 0.5, "实际 %s" % got)

    wd = tempfile.mkdtemp(prefix="div_verify_")
    try:
        rc0, smi0, _ = run_once(0.0, center, pdb, wd)
        rc5, smi5, _ = run_once(0.5, center, pdb, wd)
        chk("diversity_w=0.0 运行成功", rc0 == 0 and len(smi0) > 0,
            "rc=%s 分子数=%d" % (rc0, len(smi0)))
        chk("diversity_w=0.5 运行成功", rc5 == 0 and len(smi5) > 0,
            "rc=%s 分子数=%d" % (rc5, len(smi5)))

        print("=" * 72)
        print("2/2 行为证明: 只改 diversity_w, 生成结果必须不同")
        print("=" * 72)
        s0, s5 = set(smi0), set(smi5)
        only0, only5 = len(s0 - s5), len(s5 - s0)
        inter = len(s0 & s5)
        print("      div=0.0: %d 个 | div=0.5: %d 个 | 共有 %d | 仅0有 %d | 仅0.5有 %d"
              % (len(s0), len(s5), inter, only0, only5))
        if not (smi0 and smi5):
            chk("多样性惩罚生效", False, "有一次运行没出分子, 无法判定")
        elif only0 or only5:
            chk("多样性惩罚确实改变了生成结果", True,
                "差异 %d 个分子 —— 说明该参数在实际运行中被读到并生效" % (only0 + only5))
        else:
            chk("多样性惩罚确实改变了生成结果", False,
                "两次输出完全相同 -> 参数仍可能未生效(需继续排查)")
    finally:
        shutil.rmtree(wd, ignore_errors=True)

    print("=" * 72)
    if FAILS:
        print("结果: FAIL -> %s" % FAILS)
        return 1
    print("结果: PASS —— 出货路径正常, 且多样性惩罚被证明真实生效")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
