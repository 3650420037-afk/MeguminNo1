# -*- coding: utf-8 -*-
"""把 mass_dock_actives.py 的产物布局转换成 build_docked_dataset.py 期望的布局。

为什么要这个转换
----------------
两个脚本各自都正确, 但落盘布局不同, 不能直接对接:

    mass_dock_actives.py  产出:  outputs/mass_dock/<T>/summary_merged.csv
                                  outputs/mass_dock/<T>/docking/*.pdbqt
    build_docked_dataset.py 期望: <work>/dock_<T>/summary.csv
                                  <work>/dock_<T>/docking/*.pdbqt

比起改 build_docked_dataset.py(那是已被验证能产出 gpcr_dock_v1 的代码, 改动风险大),
这里用一次性适配器做布局转换, 并在转换时**逐个校验**「每个 ok 行都有对应位姿文件」,
缺失就报错退出 —— 否则建出来的数据集会静默少样本。

为什么不需要传 SMILES
---------------------
build_docked_dataset.py 只在**发起对接**时用 SMILES 列表; 复用已有对接结果
(--reuse-dock)时, 训练对完全由 summary 的 name + docking/<name>_out.pdbqt 位姿
构造(配体从位姿 SDF 重建)。因此本适配器只搬 name 与位姿文件即可。

位姿文件用**硬链接**而非复制(同盘符), 省空间; 失败则回退为复制。

用法
----
    python src/scripts/adapt_mass_dock.py --src outputs/mass_dock \
                                          --dst outputs/mass_dock_build
    python src/scripts/build_docked_dataset.py --work outputs/mass_dock_build \
        --reuse-dock --out data/gpcr_mass_v1
"""
import argparse
import csv
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, ensure_dir  # noqa: E402

FIELDS = ["name", "smiles", "vina_score", "status"]


def log(m):
    print(m, flush=True)


def link_or_copy(src, dst):
    if os.path.exists(dst):
        return "exists"
    try:
        os.link(src, dst)
        return "link"
    except OSError:
        shutil.copyfile(src, dst)
        return "copy"


def adapt(src_root, dst_root, targets=None, strict=True):
    src_root = src_root if os.path.isabs(src_root) else os.path.join(ROOT, src_root)
    dst_root = dst_root if os.path.isabs(dst_root) else os.path.join(ROOT, dst_root)
    ensure_dir(dst_root)

    tdirs = sorted(d for d in glob.glob(os.path.join(src_root, "*")) if os.path.isdir(d))
    if targets:
        want = set(targets)
        tdirs = [d for d in tdirs if os.path.basename(d) in want]
    if not tdirs:
        raise SystemExit("在 %s 下没找到任何靶点目录" % src_root)

    problems, grand = [], {"rows": 0, "ok": 0, "linked": 0, "missing": 0}
    for td in tdirs:
        T = os.path.basename(td)
        summ_in = os.path.join(td, "summary_merged.csv")
        if not os.path.exists(summ_in):
            problems.append("%s: 缺 summary_merged.csv(对接可能还没跑完)" % T)
            continue

        dock_in = os.path.join(td, "docking")
        dock_out = ensure_dir(os.path.join(dst_root, "dock_" + T, "docking"))
        summ_out = os.path.join(dst_root, "dock_" + T, "summary.csv")
        ensure_dir(os.path.dirname(summ_out))

        rows, n_ok, n_link, n_missing = [], 0, 0, 0
        with open(summ_in, encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                row = dict(row)
                rows.append(row)
                name = (row.get("name") or "").strip()
                if row.get("status") != "ok":
                    continue
                n_ok += 1
                psrc = os.path.join(dock_in, name + "_out.pdbqt")
                if not os.path.exists(psrc):
                    # ok 却没有位姿 = 数据不完整, 必须暴露而不是静默跳过
                    n_missing += 1
                    if n_missing <= 5:
                        problems.append("%s/%s: status=ok 但缺位姿文件" % (T, name))
                    continue
                link_or_copy(psrc, os.path.join(dock_out, name + "_out.pdbqt"))
                n_link += 1

        with open(summ_out, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

        log("  [%s] %d 行 (ok %d), 位姿就位 %d, 缺失 %d"
            % (T, len(rows), n_ok, n_link, n_missing))
        grand["rows"] += len(rows)
        grand["ok"] += n_ok
        grand["linked"] += n_link
        grand["missing"] += n_missing

    log("-" * 66)
    log("合计: %d 行, ok %d, 位姿 %d, 缺失 %d" % (grand["rows"], grand["ok"],
                                            grand["linked"], grand["missing"]))
    if problems:
        log("问题:")
        for p in problems[:20]:
            log("  ✗ %s" % p)
        if strict:
            raise SystemExit("适配失败: 存在 %d 处问题(加 --lenient 可只告警)" % len(problems))
    else:
        log("校验通过: 所有 status=ok 的行都有对应位姿文件")
    log("产物: %s/dock_<靶点>/{summary.csv,docking/}" % os.path.relpath(dst_root, ROOT))
    return 0 if not problems else 1


def main():
    ap = argparse.ArgumentParser(description="mass_dock 产物 -> build_docked_dataset 布局")
    ap.add_argument("--src", default=os.path.join(ROOT, "outputs", "mass_dock"))
    ap.add_argument("--dst", default=os.path.join(ROOT, "outputs", "mass_dock_build"))
    ap.add_argument("--targets", nargs="+", default=None)
    ap.add_argument("--lenient", action="store_true", help="有问题也只告警, 不失败")
    args = ap.parse_args()
    return adapt(args.src, args.dst, args.targets, strict=not args.lenient)


if __name__ == "__main__":
    raise SystemExit(main())
