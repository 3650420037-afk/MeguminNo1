# -*- coding: utf-8 -*-
"""数据集芳香键规范化工具 (条件性缺陷的第一层防御: 入口规范化)

用途: 当 scripts/audit_bondtypes.py 检出数据集 SDF 含芳香键标志(4) 时, 用本脚本把
芳香键 kekulize 成交替单双键, 使 Pocket2Mol 下游特征 (tri_edge 模板 / valence 加和)
与生成路径假设一致。

背景: parse_sdf_file 手工读 SDF 键块第 4 列 (1单/2双/3三/4芳香)。键值 4 会被映射为
AROMATIC(数值12), 而 tri_edge 特征模板只认 [-1,0,1,2,3]、valence 特征按键级直接加和
(12 会被当 12 价) -> 静默特征污染。本项目自建数据集实测零芳香键, 此工具为未来
引入外部数据 (如 CrossDocked) 时的预案。

用法:
  python scripts/kekulize_dataset.py --dataset data/gpcr_a2a            # 干跑(仅统计)
  python scripts/kekulize_dataset.py --dataset data/gpcr_a2a --apply    # 原地规范化(自动备份)

行为:
  - 默认 dry-run, 只报告命中数
  - --apply 时: 原文件备份为 <name>.sdf.bak (已存在则跳过备份), 再写回 kekulize 后的 MolBlock
  - 单个分子 kekulize 失败 (如无法 kekulize 的芳香体系) -> 跳过并记录, 不改动该文件
  - 写文件全程用 Python open() (规避 Windows 中文路径下 RDKit C++ 文件 API 的失败)
"""
import os, sys, argparse, pickle, shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")


def count_aromatic_bonds(path):
    """按文本键块统计芳香标志(4)的条数; 失败返回 None"""
    try:
        lines = open(path, encoding="utf-8", errors="ignore").read().splitlines()
        na, nb = int(lines[3][0:3]), int(lines[3][3:6])
        return sum(1 for bl in lines[4 + na:4 + na + nb] if int(bl[6:9]) == 4)
    except Exception:
        return None


def kekulize_file(path, apply=False):
    """返回 (状态, 芳香键数)。状态: clean/kept/fixed/skip/kekulize_fail/write_fail"""
    n_arom = count_aromatic_bonds(path)
    if n_arom is None:
        return "skip", 0
    if n_arom == 0:
        return "clean", 0
    block = open(path, encoding="utf-8", errors="ignore").read()
    try:
        m = Chem.MolFromMolBlock(block, removeHs=False, sanitize=False)
        if m is None:
            return "kekulize_fail", n_arom
        m.UpdatePropertyCache(strict=False)
        Chem.Kekulize(m, clearAromaticFlags=True)
        new_block = Chem.MolToMolBlock(m, kekulize=True)
        if count_aromatic_bonds_text(new_block) != 0:
            return "kekulize_fail", n_arom
    except Exception:
        return "kekulize_fail", n_arom
    if not apply:
        return "kept", n_arom
    try:
        bak = path + ".bak"
        if not os.path.exists(bak):
            shutil.copyfile(path, bak)
        open(path, "w", encoding="utf-8").write(new_block)
        return "fixed", n_arom
    except Exception:
        return "write_fail", n_arom


def count_aromatic_bonds_text(block):
    """对 MolBlock 文本统计芳香标志(4)"""
    try:
        lines = block.splitlines()
        na, nb = int(lines[3][0:3]), int(lines[3][3:6])
        return sum(1 for bl in lines[4 + na:4 + na + nb] if int(bl[6:9]) == 4)
    except Exception:
        return -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="数据集目录 (含 index.pkl 与 ligands/)")
    ap.add_argument("--apply", action="store_true", help="原地规范化 (自动 .bak 备份); 默认仅干跑")
    args = ap.parse_args()

    d = args.dataset
    idx = pickle.load(open(os.path.join(d, "index.pkl"), "rb"))
    from collections import Counter
    stats = Counter()
    total_arom = 0
    for _pocket, lig, _a, _b in idx:
        st, n = kekulize_file(os.path.join(d, lig), apply=args.apply)
        stats[st] += 1
        total_arom += n
    print("数据集: %s | 分子 %d" % (d, len(idx)))
    print("状态分布:", dict(stats))
    print("芳香键总数: %d" % total_arom)
    if total_arom == 0:
        print("结论: 数据洁净, 无需规范化")
    elif not args.apply:
        print("结论: 检出芳香键 -> 加 --apply 执行规范化 (会备份 .bak)")
    else:
        print("结论: 已规范化 fixed=%d, 失败=%d" % (stats.get("fixed", 0), stats.get("kekulize_fail", 0)))


if __name__ == "__main__":
    main()
