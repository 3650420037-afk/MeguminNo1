# -*- coding: utf-8 -*-
"""CPU 单测：轨迹内尺寸调度 S1（`size_schedule_threshold`）的判定逻辑。

关键不变量：**调度关闭（缺省 / target_ha<=0）时必须恒返回原阈值** ——
这是"默认行为与历史逐位一致、可一键回退"的前提。

跑法：
    python tests/test_size_schedule.py
"""
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

spec = importlib.util.spec_from_file_location(
    "sfp", os.path.join(ROOT, "src", "sample_for_pdb.py"))
sfp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sfp)
f = sfp.size_schedule_threshold


def test_off_is_identity():
    for sched in (None, {}, {'target_ha': 0}, {'relax': 0.5}):
        for n_ha, step in ((3, 1), (20, 20), (40, 39)):
            thr, act = f(0.0, n_ha, step, 40, sched)
            assert thr == 0.0 and act == 'off', (sched, thr, act)
    thr, act = f(-0.5, 10, 10, 40, {})
    assert thr == -0.5 and act == 'off'
    print("PASS 调度关闭时恒为原阈值（可一键回退的前提）")


def test_warmup_no_intervention():
    sched = {'target_ha': 28, 'warmup_frac': 0.25, 'relax': 0.35, 'slack': 4.0}
    # step 5/40 = 0.125 < 0.25 → 即使严重落后也不动
    thr, act = f(0.0, 1, 5, 40, sched)
    assert thr == 0.0 and act == 'warmup', (thr, act)
    print("PASS warmup 期间不干预（即使严重落后）")


def test_relax_when_behind():
    sched = {'target_ha': 28, 'warmup_frac': 0.25, 'relax': 0.35, 'slack': 4.0}
    # step 20/40 → prog 0.5 → expect 14；落后到 5（< 14-4=10）→ 下调
    thr, act = f(0.0, 5, 20, 40, sched)
    assert abs(thr + 0.35) < 1e-9 and act == 'relax', (thr, act)
    # 刚好在容差边界：n=10 → 不小于 expect-slack → 不动
    thr2, act2 = f(0.0, 10, 20, 40, sched)
    assert abs(thr2) < 1e-9 and act2 == 'keep', (thr2, act2)
    print("PASS 落后超容差才下调；恰在边界不动")


def test_tighten_only_when_enabled():
    base = {'target_ha': 28, 'warmup_frac': 0.25, 'relax': 0.35, 'slack': 4.0}
    # 超前：n=30 vs expect 14 → 默认 tighten=0 → 不动
    thr, act = f(0.0, 30, 20, 40, base)
    assert abs(thr) < 1e-9 and act == 'keep', (thr, act)
    # 开启 tighten=0.25 → 上调
    sched = dict(base, tighten=0.25)
    thr2, act2 = f(0.0, 30, 20, 40, sched)
    assert abs(thr2 - 0.25) < 1e-9 and act2 == 'tighten', (thr2, act2)
    print("PASS tighten 默认关闭；开启后按幅度上调")


def test_expected_scales_with_progress():
    """目标随步数线性外推：同一个小分子，早步不罚、晚步才罚。"""
    sched = {'target_ha': 40, 'warmup_frac': 0.0, 'relax': 0.5, 'slack': 2.0}
    early = f(0.0, 5, 5, 40, sched)     # prog .125 → expect 5 → keep
    late = f(0.0, 5, 35, 40, sched)     # prog .875 → expect 35 → relax
    assert early[1] == 'keep' and late[1] == 'relax', (early, late)
    assert abs(late[0] + 0.5) < 1e-9
    print("PASS 期望值随进度线性外推（早步不罚、晚步才罚）")


if __name__ == "__main__":
    test_off_is_identity()
    test_warmup_no_intervention()
    test_relax_when_behind()
    test_tighten_only_when_enabled()
    test_expected_scales_with_progress()
    print("ALL PASS (5/5)")
