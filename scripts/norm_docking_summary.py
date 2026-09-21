# -*- coding: utf-8 -*-
"""规范化已有 docking summary.csv 的 smiles 列 (修复历史污染, 无需重跑对接)

修复两类污染 (与 docking_pipeline.py 的历史实现缺陷对应):
  1) U+FEFF BOM 前缀: --smiles-file 旧写法用 utf-8 读取, strip() 无法去除 BOM,
     导致每个库的首行 smiles 带不可见字符, 与 compounds.csv 主键不一致而静默丢样
  2) 显式氢 SMILES: --sdf-dir 分支旧写法 MolToSmiles(含氢分子), 输出
     [H]OC([H])([H])... 形式, 与库规范 SMILES 无法关联

处理: smiles 去 BOM; 再用 RDKit canonical 规范化 (RemoveHs) —— 仅当结果变化时替换。
安全: 原文件先备份为 <name>.bak (已存在则不覆盖备份); 分数/状态等其它列原样保留。

用法: python scripts/norm_docking_summary.py [--dry]
"""
import os, sys, csv, glob, shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

OUTROOT = r"D:\MMModel\Pocket2Mol\outputs"
LIB = r"D:\MMModel\化合物库"


def norm_smiles(s):
    if not s:
        return s, "empty"
    s2 = s.lstrip("\ufeff")
    tag = "bom" if s2 != s else ""
    try:
        m = Chem.MolFromSmiles(s2)
        if m is None:
            return s2, (tag + "+unparsable") or "unparsable"
        can = Chem.MolToSmiles(Chem.RemoveHs(m))
        if can != s2:
            tag = (tag + "+canon") if tag else "canon"
        return can, tag
    except Exception:
        return s2, (tag + "+err") or "err"


def main():
    dry = "--dry" in sys.argv
    files = []
    for pat in ("docking_*/summary.csv", "*_library/summary.csv", "*/summary.csv"):
        files += glob.glob(os.path.join(OUTROOT, pat))
    files = sorted(set(files))
    print("待检查 summary.csv: %d 份  (dry=%s)" % (len(files), dry))
    tot_fix = 0
    for f in files:
        try:
            rows = list(csv.DictReader(open(f, encoding="utf-8-sig")))
        except Exception as e:
            print("  跳过(读取失败) %s: %s" % (f, e)); continue
        if not rows or "smiles" not in rows[0]:
            continue
        fixed = 0
        stats = {}
        for r in rows:
            new, tag = norm_smiles(r["smiles"])
            if new != r["smiles"]:
                fixed += 1
                stats[tag] = stats.get(tag, 0) + 1
                r["smiles"] = new
        if fixed == 0:
            continue
        rel = os.path.relpath(f, OUTROOT)
        print("  %-58s 修复 %d 行  %s" % (rel, fixed, stats))
        if not dry:
            bak = f + ".bak"
            if not os.path.exists(bak):
                shutil.copyfile(f, bak)
            with open(f, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader(); w.writerows(rows)
        tot_fix += fixed
    print("合计修复行数: %d" % tot_fix)

    # 与库 CSV 的主键交集变化 (验证丢样是否恢复)
    print()
    print("=== 与库 compounds.csv 的交集核对 ===")
    for t, comp in (("A2A", os.path.join(LIB, "compounds.csv")),
                    ("b2ar", os.path.join(LIB, "b2ar", "compounds.csv")),
                    ("d3", os.path.join(LIB, "d3", "compounds.csv")),
                    ("5ht2b", os.path.join(LIB, "5ht2b", "compounds.csv"))):
        if not os.path.exists(comp):
            continue
        lib = set(r["smiles"] for r in csv.DictReader(open(comp, encoding="utf-8")))
        cands = glob.glob(os.path.join(OUTROOT, "docking_%s_library" % t, "summary.csv")) or \
                glob.glob(os.path.join(OUTROOT, "docking_%s_library" % t.lower(), "summary.csv"))
        if not cands:
            continue
        dsum = [r for r in csv.DictReader(open(cands[0], encoding="utf-8-sig")) if r["status"] == "ok"]
        inter = sum(1 for r in dsum if r["smiles"] in lib)
        print("  %-6s 库 %4d | 对接 %4d | 交集 %4d %s" % (
            t, len(lib), len(dsum), inter, "(完整)" if inter == len(dsum) else "(仍有 %d 未匹配)" % (len(dsum) - inter)))


if __name__ == "__main__":
    main()
