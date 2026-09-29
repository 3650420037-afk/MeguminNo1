# -*- coding: utf-8 -*-
"""CPU 单测：A′ 采样器 v1.2 的**内部超时**行为（超时必须标 invalid 而不是挂死/崩掉）。

背景：§5.39 那次"thr −0.50 得 n=0"实际是外层 900 s 超时把批次杀掉，被误记成 INVALID。
v1.2 起 `sample_once(..., timeout=)` 用内部超时：超时 → 杀子进程 → rc=124 → 该次标 invalid、
阈值不变、同阈值重试；并在 run.log 里留下 `[TIMEOUT]` 证据。

跑法：
    python tests/test_adaptive_timeout.py
"""
import importlib.util
import os
import subprocess
import sys

try:  # 无人值守/重定向下 stdout 可能是 GBK 管道：非 GBK 字符会直接崩（本项目已真实踩过）
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

spec = importlib.util.spec_from_file_location(
    "ats", os.path.join(ROOT, "src", "scripts", "adaptive_threshold_sample.py"))
ats = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ats)


def test_timeout_marks_invalid_and_logs():
    calls = {}

    def fake_run(cmd, **kw):
        calls["timeout"] = kw.get("timeout")
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=kw.get("timeout"))

    real = subprocess.run
    subprocess.run = fake_run
    try:
        with tempfile.TemporaryDirectory() as td:
            smis, sess, rc = ats.sample_once("/nonexistent.pt", "A2A", 5, 50, 40,
                                             2024, -0.25, td, timeout=17.0)
            log = open(os.path.join(td, "run.log"), encoding="utf-8").read()
    finally:
        subprocess.run = real

    assert calls.get("timeout") == 17.0, "timeout 没有传到 subprocess.run: %s" % calls
    assert rc == 124, "超时应返回 rc=124，实际 %s" % rc
    assert smis == [] and sess is None, (smis, sess)
    assert "[TIMEOUT]" in log, "run.log 里没有留下 TIMEOUT 证据"
    print("PASS 超时 → rc=124 / 无产物 / run.log 留痕，且 timeout 透传到子进程")


def test_no_timeout_passes_none():
    """默认（timeout=None 或 0→None）时不设超时，行为与旧版一致。"""
    calls = {}

    class R:
        returncode = 0

    def fake_run(cmd, **kw):
        calls["timeout"] = kw.get("timeout")
        return R()

    real = subprocess.run
    subprocess.run = fake_run
    try:
        with tempfile.TemporaryDirectory() as td:
            ats.sample_once("/nonexistent.pt", "A2A", 5, 50, 40, 2024, 0.0, td, timeout=None)
    finally:
        subprocess.run = real
    assert calls.get("timeout") is None, calls
    print("PASS 不传 timeout 时 subprocess.run 收到 None（不设超时，向后兼容）")


if __name__ == "__main__":
    test_timeout_marks_invalid_and_logs()
    test_no_timeout_passes_none()
    print("ALL PASS (2/2)")
