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
    # 明显偏小（5 个原子 << 24）→ 下调
    thr, act = f(0.0, 5, 20, 40, sched)
    assert abs(thr + 0.35) < 1e-9 and act == 'relax', (thr, act)
    # 恰在容差边界：n=24 → 不小于 target-slack → 不动
    thr2, act2 = f(0.0, 24, 20, 40, sched)
    assert abs(thr2) < 1e-9 and act2 == 'keep', (thr2, act2)
    print("PASS 偏小超容差才下调；恰在边界不动")


def test_tighten_only_when_enabled():
    base = {'target_ha': 28, 'warmup_frac': 0.25, 'relax': 0.35, 'slack': 4.0}
    # 偏大：n=40 > 32（target+slack）→ 默认 tighten=0 → 不动
    thr, act = f(0.0, 40, 20, 40, base)
    assert abs(thr) < 1e-9 and act == 'keep', (thr, act)
    # 开启 tighten=0.25 → 上调
    sched = dict(base, tighten=0.25)
    thr2, act2 = f(0.0, 40, 20, 40, sched)
    assert abs(thr2 - 0.25) < 1e-9 and act2 == 'tighten', (thr2, act2)
    print("PASS tighten 默认关闭；开启后按幅度上调")


def test_absolute_size_thermostat():
    """修正后的机制：**绝对尺寸恒温器**（不再按步数线性外推）。

    依据（实测 §5.47）：8 份日志 410 个分子 `重原子数 − 步数 ≡ 0` —— 每步恰好加 1 个重原子，
    真实轨迹是 1×step，永远高于"target×step/max_steps"的线性期望 → 旧规则在正常 40 步运行里
    几乎永不触发。故改为按**绝对尺寸**判断。
    """
    sched = {'target_ha': 28, 'warmup_frac': 0.25, 'relax': 0.35, 'slack': 4.0}
    # 正常轨迹（1 原子/步）：step 30 → n=30 > 28-4=24 → 不干预
    keep = f(0.0, 30, 30, 40, sched)
    assert keep == (0.0, 'keep'), keep
    # 掉进小模式：step 30 却只有 15 个原子（远小于 24）→ 放松阈值
    relax = f(0.0, 15, 30, 40, sched)
    assert abs(relax[0] + 0.35) < 1e-9 and relax[1] == 'relax', relax
    # 偏小但还在 warmup 内（step 5 < 0.25*40=10）→ 不干预
    warm = f(0.0, 3, 5, 40, sched)
    assert warm == (0.0, 'warmup'), warm
    # 略小于目标但在 slack 内（n=25 vs 24 阈值）→ 不干预
    edge = f(0.0, 25, 30, 40, sched)
    assert edge == (0.0, 'keep'), edge
    print("PASS 绝对尺寸恒温器：偏小才放松、正常轨迹不干预、warmup 内不动、边界不误触")


def test_real_trajectory_never_relaxes_when_on_target():
    """回归固化：对'每步 1 原子'的真实轨迹（n == step），只要按目标尺寸生长就不会被干预。"""
    sched = {'target_ha': 28, 'warmup_frac': 0.25, 'relax': 0.35, 'slack': 4.0}
    for step in range(1, 41):
        n = step                      # 实测：每步恰好 1 个重原子
        thr, act = f(0.0, n, step, 40, sched)
        if step < 10:
            assert act == 'warmup', (step, act)
        elif n >= 24:
            assert act == 'keep', (step, act)
        else:
            assert act == 'relax', (step, act)
    print("PASS 真实 1 原子/步 轨迹下：≥24 原子后不再被干预（只在真偏小时放松）")


if __name__ == "__main__":
    test_off_is_identity()
    test_warmup_no_intervention()
    test_relax_when_behind()
    test_tighten_only_when_enabled()
    test_absolute_size_thermostat()
    test_real_trajectory_never_relaxes_when_on_target()
    print("ALL PASS (6/6)")
