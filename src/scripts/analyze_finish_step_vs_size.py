# -*- coding: utf-8 -*-
"""**S1 前提的直接检验**：早完成的分子是否就是小分子？

为什么先做这个（跑之前验证，比跑完再解释强）：
S1（§5.46）的假设是"分子在前几步就落后 = 掉进小模式，且再也回不来"。若"完成步数"与
"最终尺寸"无关，那 S1 的干预时机假设就不成立，不该花 GPU 去跑它。
数据来源：既有采样会话的 `sample.log` —— 里面既有每步一行 `[Pool] Queue/Finished/Failed`，
也有完成时的 `Success: <SMILES>`。**按 Success 之前出现的 [Pool] 行数即可还原"第几步完成"**。

用法：
    python src/scripts/analyze_finish_step_vs_size.py <sample.log> [<sample.log> ...]
"""
import argparse
import os
import re
import statistics as st
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

POOL = re.compile(r"\[Pool\] Queue (\d+) \| Finished (\d+) \| Failed (\d+)")
SUCC = re.compile(r"Success: (\S+)")


def parse(path):
    """返回 [(step, smiles, heavy_atoms)]，step 为完成时的步序号（1 起）。"""
    step = 0
    out = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            m = POOL.search(line)
            if m:
                step += 1
                continue
            m2 = SUCC.search(line)
            if m2:
                smi = m2.group(1)
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                out.append((step, smi, mol.GetNumHeavyAtoms()))
    return out


def spearman(xs, ys):
    def rank(v):
        idx = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(idx):
            j = i
            while j + 1 < len(idx) and v[idx[j + 1]] == v[idx[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[idx[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    args = ap.parse_args()
    pooled_x, pooled_y = [], []
    print("%-52s %5s %7s %7s %8s %s" % ("log", "n", "HA中位", "首/末HA", "rho", "最早3个"))
    for p in args.logs:
        if not os.path.exists(p):
            print("跳过(不存在): %s" % p)
            continue
        recs = parse(p)
        if len(recs) < 5:
            print("跳过(有效分子<5): %s" % p)
            continue
        steps = [r[0] for r in recs]
        has = [r[2] for r in recs]
        rho = spearman(steps, has)
        name = os.path.relpath(p, os.path.dirname(os.path.dirname(os.path.dirname(p))))
        print("%-52s %5d %7.0f %3d/%-3d %+8.3f %s"
              % (name[:52], len(recs), st.median(has), min(has), max(has), rho,
                 ", ".join("%d(HA%d)" % (s, h) for s, _, h in recs[:3])))
        pooled_x += steps
        pooled_y += has
    if pooled_x:
        print("\n池化（全部日志）: n=%d | Spearman(完成步, 重原子数) = **%+.3f**"
              % (len(pooled_x), spearman(pooled_x, pooled_y)))
        # 早完成 vs 晚完成两组的尺寸对比（按步数中位切）
        med = st.median(pooled_x)
        early = [h for s, h in zip(pooled_x, pooled_y) if s <= med]
        late = [h for s, h in zip(pooled_x, pooled_y) if s > med]
        print("早完成组（步≤%.0f）HA 中位 %.0f（n=%d） vs 晚完成组 HA 中位 %.0f（n=%d）→ 差 %+.0f"
              % (med, st.median(early), len(early), st.median(late), len(late),
                 st.median(early) - st.median(late)))
        print("\n判读：rho 显著为正且早完成组明显更小 → **S1 的前提成立**（落后=小模式，值得按轨迹干预）；")
        print("      rho ≈ 0 或为负 → S1 干预时机假设不成立，应先改设计再花 GPU。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
