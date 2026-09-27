# -*- coding: utf-8 -*-
"""实测"建数据集"这一步的速率, 把估算里最大的不确定项钉死。

为什么单独测它: 15,000 个位姿要逐个 PDBQT->SDF(obabel) + RDKit 清洗 + 穿模校验,
我此前只能给 30-90 分钟的区间。用**仓库自己的转换函数**跑一批真实位姿,
按实测每分子耗时外推, 才能给出可用数字。

用法: python src/scripts/measure_dataset_build_rate.py --target A2A --n 40
"""
import argparse
import glob
import os
import statistics as st
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, src_on_path  # noqa: E402
src_on_path()
import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "bdd", os.path.join(ROOT, "src", "scripts", "build_docked_dataset.py"))
bdd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bdd)          # 只加载函数, 不执行 main(有 __main__ 保护)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="A2A")
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--total", type=int, default=15000, help="外推用的总位姿数")
    args = ap.parse_args()

    poses = sorted(glob.glob(os.path.join(ROOT, "outputs", "mass_dock",
                                         args.target, "docking", "*_out.pdbqt")))[:args.n]
    if not poses:
        raise SystemExit("找不到位姿: outputs/mass_dock/%s/docking" % args.target)
    tmpd = os.path.join(ROOT, "outputs", "_rate_test")
    os.makedirs(tmpd, exist_ok=True)

    ok, tms = 0, []
    for i, p in enumerate(poses):
        out = os.path.join(tmpd, "t%03d.sdf" % i)
        t0 = time.perf_counter()
        good = bdd.pdbqt_to_sdf(p, out, tmpd)
        dt = time.perf_counter() - t0
        if good:
            ok += 1
            tms.append(dt)
    print("实测 %s: %d/%d 转换成功" % (args.target, ok, len(poses)))
    if not tms:
        raise SystemExit("全部转换失败, 无法估速率")
    med = st.median(tms)
    print("  每分子转换耗时: 中位 %.3f s  均值 %.3f s  最大 %.3f s"
          % (med, sum(tms) / len(tms), max(tms)))
    for tot in (args.total, 7249, 15000):
        sec = med * tot
        print("  外推 %6d 个位姿 -> %5.1f 分钟 (中位速率)" % (tot, sec / 60.0))
    print()
    print("注: 这只是**转换**耗时; 建库还要做穿模校验(RDKit 几何, 很快)与写 SDF,")
    print("    以及 `adapt_mass_dock` 的硬链接(秒级)。故实际建库时间≈上表 + 2-3 分钟。")
    for f in glob.glob(os.path.join(tmpd, "*")):
        os.remove(f)
    os.rmdir(tmpd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
