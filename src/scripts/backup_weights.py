# -*- coding: utf-8 -*-
"""给"权重 + 配置 + 路径解析"打一个带 md5 的可回退备份点。

为什么需要它
------------
用户确认的路线里明确要求：**每次结构改动前先备份当前权重与配置以便回退**。
手工 copy 容易漏（这次就漏过一次——只备份了 .pt，没备份 configs 与 paths.py，
而 `DEFAULT_CKPT` 就在 paths.py 里），所以固化成脚本。

备份内容
--------
- `models/*.pt`（全部权重：官方基线 + 各版自训）
- `models/README.md`
- `configs/**`（全部训练/采样配置）
- `src/paths.py`（`DEFAULT_CKPT` 的解析处）
- `MANIFEST.json`：每个文件的字节数与 md5 —— 回退后可比对是否与原状一致

输出位置
--------
默认 `other/weight_backups/<时间戳>/`。`other/` 按 .gitignore **不随仓库提交**，
因此这是**本机回退点**，不污染交付物，也不会撑大仓库。
（权重单文件 15–45 MB，多次备份会占盘；脚本会报告新备份的体积与 `other/` 下
累计备份数，便于清理。）

用法
----
    python src/scripts/backup_weights.py            # 打一个备份点
    python src/scripts/backup_weights.py --list     # 只列出已有备份点
"""
import argparse
import glob
import hashlib
import json
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT  # noqa: E402

DEST_ROOT = os.path.join(ROOT, "other", "weight_backups")


def human(n):
    return "%.1f MB" % (n / 1e6) if n >= 1e6 else "%.1f KB" % (n / 1e3)


def main():
    ap = argparse.ArgumentParser(description="权重/配置可回退备份点")
    ap.add_argument("--dest-root", default=DEST_ROOT)
    ap.add_argument("--list", action="store_true", help="只列出已有备份点")
    ap.add_argument("--tag", default="", help="附加到目录名的说明")
    args = ap.parse_args()

    if args.list:
        ds = sorted(glob.glob(os.path.join(args.dest_root, "*")))
        if not ds:
            print("尚无备份点: %s" % args.dest_root)
            return 0
        print("已有 %d 个备份点:" % len(ds))
        for d in ds:
            if not os.path.isdir(d):
                continue
            n = sum(len(f) for _, _, f in os.walk(d))
            sz = sum(os.path.getsize(os.path.join(r, f))
                     for r, _, fs in os.walk(d) for f in fs)
            print("  %-24s %3d 个文件  %s" % (os.path.basename(d), n, human(sz)))
        return 0

    stamp = time.strftime("%Y%m%d_%H%M%S") + (("_" + args.tag) if args.tag else "")
    dst = os.path.join(args.dest_root, stamp)
    os.makedirs(dst, exist_ok=True)

    made = []

    def put(src, rel):
        """copy2 到备份点的 rel 位置; 支持 glob。"""
        for p in sorted(glob.glob(os.path.join(ROOT, src))):
            if not os.path.isfile(p):
                continue
            t = os.path.join(dst, rel, os.path.basename(p))
            os.makedirs(os.path.dirname(t), exist_ok=True)
            shutil.copy2(p, t)
            made.append(t)

    put("models/*.pt", "")
    put("models/README.md", "")
    # configs 保持目录结构(内含子目录时也不丢)
    for p in sorted(glob.glob(os.path.join(ROOT, "configs", "**", "*"), recursive=True)):
        if os.path.isfile(p):
            rel = os.path.relpath(p, os.path.join(ROOT, "configs"))
            t = os.path.join(dst, "configs", rel)
            os.makedirs(os.path.dirname(t), exist_ok=True)
            shutil.copy2(p, t)
            made.append(t)
    put("src/paths.py", "")

    manifest = {}
    for p in made:
        manifest[os.path.relpath(p, dst)] = {
            "bytes": os.path.getsize(p),
            "md5": hashlib.md5(open(p, "rb").read()).hexdigest(),
        }
    with open(os.path.join(dst, "MANIFEST.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=1)

    weights = {k: v for k, v in manifest.items() if k.endswith(".pt")}
    total = sum(v["bytes"] for v in manifest.values())
    print("备份点: %s" % os.path.relpath(dst, ROOT))
    print("  文件 %d 个（权重 %d 个），合计 %s" % (len(manifest), len(weights), human(total)))
    for k in sorted(weights):
        print("    %-40s %s" % (k, human(weights[k]["bytes"])))
    print("  校验清单: MANIFEST.json（含每个文件的 md5，回退后可逐项比对）")
    others = [d for d in glob.glob(os.path.join(args.dest_root, "*")) if os.path.isdir(d)]
    if len(others) > 1:
        allsz = sum(os.path.getsize(os.path.join(r, f))
                    for r, _, fs in os.walk(args.dest_root) for f in fs)
        print("  ⚠ 备份点累计 %d 个，占 %s —— 权重较大，请定期清理旧备份点"
              % (len(others), human(allsz)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
