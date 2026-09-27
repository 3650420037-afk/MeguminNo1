# -*- coding: utf-8 -*-
"""把仓库里非 UTF-8 的文本日志就地规范化为 UTF-8。

为什么需要它
------------
`src/utils/misc.py` 的 `get_logger` 原先没有给 `FileHandler` 指定 `encoding`,
在中文 Windows 上会按 locale 默认编码(GBK)写文件。结果是**同一份日志里英文行是 ASCII、
中文行是 GBK** —— 这类文件不是合法 UTF-8: grep / read / CI / 编辑器都会看到乱码,
甚至直接报 "line is not valid UTF-8"(实测差点因此漏看关键诊断行, 见
`docs/任务状态与记忆.md` §4-⑥)。代码侧已修(显式 encoding='utf-8'), 但**历史日志仍是 GBK**,
作为交付物留在仓库里就是不达标文件。

本脚本只做编码规范化, **不改变任何语义内容**: 逐行保留, 只把字节序换掉。
已经是合法 UTF-8 的文件**原样跳过**(幂等, 可重复运行)。

安全措施
--------
- 默认跳过最近 `--min-age` 秒内被修改过的文件 —— 正在被写作进程追加的日志不能改,
  否则会出现"前半 UTF-8、后半 GBK"的更糟状态;  等进程结束再跑一次即可。
- `--dry-run` 只报告不写入。

用法
----
    python src/scripts/normalize_log_encoding.py --dry-run
    python src/scripts/normalize_log_encoding.py
"""
import argparse
import glob
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT  # noqa: E402

# 候选源编码: GBK 是中文 Windows 的 locale 默认值(实测病根); cp936 与 GBK 等价,
# 这里再兜一个 latin-1(永不失败)作为最后手段, 但会单独标记出来以免静默损坏内容。
FALLBACK_ENCODINGS = ["gbk", "cp936", "big5", "latin-1"]


def decode_any(raw):
    """返回 (text, encoding)。先试 UTF-8, 再依次试其它编码。"""
    try:
        return raw.decode("utf-8"), "utf-8"
    except UnicodeDecodeError:
        pass
    for enc in FALLBACK_ENCODINGS:
        try:
            return raw.decode(enc), enc
        except UnicodeDecodeError:
            continue
    return None, None


def iter_targets(patterns):
    seen = set()
    for pat in patterns:
        for p in glob.glob(os.path.join(ROOT, pat), recursive=True):
            if os.path.isfile(p) and p not in seen:
                seen.add(p)
                yield p


def main():
    ap = argparse.ArgumentParser(description="把非 UTF-8 文本日志规范化为 UTF-8")
    ap.add_argument("--patterns", nargs="+",
                    default=["logs/**/log.txt", "logs/**/*.jsonl", "outputs/**/sample.log"])
    ap.add_argument("--min-age", type=float, default=180.0,
                    help="跳过最近 N 秒内被修改的文件(默认 180, 避免改到正在写入的日志)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    conv, skipped_utf8, skipped_fresh, latin = [], 0, [], []
    for p in iter_targets(args.patterns):
        age = time.time() - os.path.getmtime(p)
        raw = open(p, "rb").read()
        if not raw:
            continue
        try:
            raw.decode("utf-8")
            skipped_utf8 += 1
            continue
        except UnicodeDecodeError:
            pass
        if age < args.min_age:
            skipped_fresh.append((p, age))
            continue
        text, enc = decode_any(raw)
        if text is None:
            print("  [跳过] 无法解码: %s" % os.path.relpath(p, ROOT))
            continue
        if enc == "latin-1":
            # latin-1 永不失败, 但可能把真正的 GBK 解成乱码 —— 必须让人看见
            latin.append(p)
        rel = os.path.relpath(p, ROOT)
        print("  %-58s %s -> utf-8" % (rel, enc))
        if not args.dry_run:
            with open(p, "w", encoding="utf-8", newline="") as f:
                f.write(text)
        conv.append((rel, enc))

    print("-" * 72)
    print("转换 %d 个 | 已是 UTF-8 跳过 %d 个 | 正在写入而跳过 %d 个"
          % (len(conv), skipped_utf8, len(skipped_fresh)))
    for p, age in skipped_fresh:
        print("  [正在写入, %.0f 秒前] %s" % (age, os.path.relpath(p, ROOT)))
    if latin:
        print("⚠ 以下文件用 latin-1 兜底解码, 请人工确认没有损坏:")
        for p in latin:
            print("   ", os.path.relpath(p, ROOT))
    if args.dry_run:
        print("(dry-run, 未写入)")
    # 若仍有"正在写入"的文件, 提示稍后重跑 —— 否则仓库里会残留非 UTF-8 文件
    return 0 if not skipped_fresh else 2


if __name__ == "__main__":
    raise SystemExit(main())
