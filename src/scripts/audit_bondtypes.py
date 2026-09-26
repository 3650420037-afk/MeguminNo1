# -*- coding: utf-8 -*-
"""数据集键类型洁净性审计 (问题5 的取证工具)

背景: Pocket2Mol 的 parse_sdf_file 手工解析 SDF 文本键块第 4 列数字
(1/单 2/双 3/三 4/芳香) 并直接映射为训练用 bond_type。若数据含芳香标志(4),
会导致: ①tri_edge 特征模板 [-1,0,1,2,3] 覆盖不到而退化为全零行;
        ②valence 特征按键级加和 +12 污染。
本脚本按"与 parse_sdf_file 一致的口径"审计数据集的真实键值分布。

重要方法论陷阱: RDKit MolFromMolBlock(sanitize=True) 会**自动芳香化**交替单双键,
GetBondType() 大量返回 AROMATIC (本数据实测 52.7%) — 这**不能**代表文件里实际存储的
键值。必须解析 SDF 文本键块原始数字 (本脚本做法)。

用法: python scripts/audit_bondtypes.py [数据集目录名...]
"""
import os, sys, pickle
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import DATA, src_on_path
src_on_path()

DEFAULT = ["gpcr_a2a", "gpcr_b2ar", "gpcr_d3", "gpcr_5ht2b", "gpcr_multitarget_v2"]


def sdf_bond_types(path):
    """按 SDF 文本键块第 7-9 列原始数字统计 (与 utils/protein_ligand.parse_sdf_file 口径一致)"""
    try:
        lines = open(path, encoding="utf-8", errors="ignore").read().splitlines()
        na, nb = int(lines[3][0:3]), int(lines[3][3:6])
        c = Counter()
        for bl in lines[4 + na:4 + na + nb]:
            c[int(bl[6:9])] += 1
        return c
    except Exception:
        return None


def main(datasets):
    total = Counter()
    n_skipped = 0        # 解析失败计数: 不能静默跳过, 否则"洁净"结论的分母会失真
    print("=== SDF 文本原始键值分布（parse_sdf_file 口径）===")
    for ds in datasets:
        d = os.path.join(DATA, ds)
        idxp = os.path.join(d, "index.pkl")
        if not os.path.exists(idxp):
            continue
        idx = pickle.load(open(idxp, "rb"))
        c, n = Counter(), 0
        for pocket, lig, _, _ in idx:
            r = sdf_bond_types(os.path.join(d, lig))
            if r is None:
                n_skipped += 1          # 显式计数, 避免"洁净"结论分母漏项
                print("    WARN 解析失败, 已跳过: %s" % lig)
                continue
            n += 1
            c.update(r)
        total.update(c)
        print("%-22s 分子 %4d | 单 %5d 双 %5d 三 %4d 芳香(4) %d" % (
            ds, n, c.get(1, 0), c.get(2, 0), c.get(3, 0), c.get(4, 0)))
    s = sum(total.values())
    print()
    if s == 0:
        print("汇总: 无有效键数据 (全部解析失败或数据集为空), 无法判定洁净性")
        return
    print("汇总 %d 键: 单 %d (%.1f%%) | 双 %d (%.1f%%) | 三 %d (%.1f%%) | 芳香(4) %d" % (
        s, total.get(1, 0), 100.0 * total.get(1, 0) / s,
        total.get(2, 0), 100.0 * total.get(2, 0) / s,
        total.get(3, 0), 100.0 * total.get(3, 0) / s, total.get(4, 0)))
    if n_skipped:
        print("注意: 有 %d 个 SDF 解析失败被跳过, 洁净性结论未覆盖这些文件" % n_skipped)
    print("结论:", "数据洁净（无芳香标志，训练特征无污染）" if total.get(4, 0) == 0
          else "存在芳香键 %d 条 -> 需 kekulize 规范化后再训练" % total.get(4, 0))


if __name__ == "__main__":
    main(sys.argv[1:] or DEFAULT)
