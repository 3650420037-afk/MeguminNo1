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

try:  # 无人值守/重定向下 stdout 可能是 GBK 管道：非 GBK 字符会直接崩（本项目已真实踩过）
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

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


def test_floor_s2_hard_size_floor():
    """S2 尺寸下限：前 floor_steps 步一律用 floor_thr（方向由构造保证）。

    依据（§5.47）：每步恰好 1 个重原子 → "前 N 步强制继续长" ≡ "分子至少 N 个重原子"。
    这比阈值微调可靠：实测阈值→尺寸**并不单调**（v10 在 thr −0.25 时 HA 反而 25/19，§5.48）。
    """
    sched = {'floor_steps': 24, 'floor_thr': -2.0, 'target_ha': 28,
             'warmup_frac': 0.25, 'relax': 0.35, 'slack': 4.0}
    # 第 1..23 步：无论大小，一律用 floor_thr
    for step in (1, 10, 23):
        thr, act = f(0.0, step, step, 40, sched)
        assert act == 'floor' and abs(thr + 2.0) < 1e-9, (step, thr, act)
    # 第 24 步起：回到恒温器逻辑（24 ≥ 28-4 → keep）
    thr2, act2 = f(0.0, 24, 24, 40, sched)
    assert act2 == 'keep' and abs(thr2) < 1e-9, (thr2, act2)
    # 只开 S2（target_ha=0，恒温器关闭）也要生效
    sched2 = {'floor_steps': 5, 'floor_thr': -1.5}
    thr3, act3 = f(0.0, 2, 2, 40, sched2)
    assert act3 == 'floor' and abs(thr3 + 1.5) < 1e-9, (thr3, act3)
    thr4, act4 = f(0.0, 6, 6, 40, sched2)
    assert act4 == 'off' and abs(thr4) < 1e-9, (thr4, act4)
    print("PASS S2 下限期：前 N 步强制 floor_thr；之后回到恒温器/关闭；可单独使用")


def test_floor_takes_precedence_over_thermostat():
    """下限期优先于恒温器（避免两个机制在同一步互相覆盖）。"""
    sched = {'floor_steps': 10, 'floor_thr': -2.0, 'target_ha': 28,
             'warmup_frac': 0.25, 'relax': 0.35, 'slack': 4.0}
    thr, act = f(0.0, 3, 3, 40, sched)   # 又小又在限期内 → 应判 floor
    assert act == 'floor' and abs(thr + 2.0) < 1e-9, (thr, act)
    print("PASS 下限期优先于恒温器")


def _run_sampler_with_schedule(sched, max_steps=4, timeout=300):
    """用给定 schedule 直接跑采样器子进程，返回 (rc, stdout+stderr)。"""
    import io as _io
    import subprocess
    import tempfile
    import yaml
    cfg = yaml.safe_load(_io.open(os.path.join(ROOT, "configs", "sample_for_pdb_guided_l3.yml"),
                                  encoding="utf-8"))
    cfg.setdefault("sample", {})
    cfg["sample"]["max_steps"] = max_steps
    cfg["sample"]["num_samples"] = 1
    cfg["sample"]["beam_size"] = 2
    cfg["sample"].setdefault("threshold", {})["frontier_threshold_schedule"] = sched
    cfg.setdefault("model", {})["checkpoint"] = os.path.join(ROOT, "models",
                                                             "7-eonmol_ft_gpcr_v10.pt")
    p = os.path.join(tempfile.mkdtemp(prefix="sched_bad_"), "cfg.yml")
    with _io.open(p, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, allow_unicode=True)
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="")
    r = subprocess.run([sys.executable, os.path.join(ROOT, "src", "sample_for_pdb.py"),
                        "--pdb_path", os.path.join(ROOT, "data", "targets", "4EIY_A2A受体.pdb"),
                        "--center=0,0,0", "--config", p, "--device", "cpu",
                        "--outdir", os.path.join(os.path.dirname(p), "out")],
                       cwd=ROOT, capture_output=True, text=True, errors="ignore",
                       timeout=timeout, env=env)
    return r.returncode, (r.stdout or "") + (r.stderr or "")


def test_floor_ge_max_steps_fails_fast():
    """floor_steps >= max_steps 必须**快速失败**，而不是静默产出 0 个分子（独立审查发现的 D5）。"""
    rc, out = _run_sampler_with_schedule(
        {"target_ha": 0, "slack": 4.0, "relax": 0.35, "tighten": 0.0,
         "warmup_frac": 0.25, "floor_steps": 4.0, "floor_thr": -2.0, "ceiling_ha": 0.0},
        max_steps=4)
    assert rc != 0, "应当以非零码退出"
    assert "floor_steps" in out and "max_steps" in out, out[-300:]
    print("PASS S4 floor_steps>=max_steps 快速失败（rc=%d，含明确错误信息）" % rc)


def test_tiny_ceiling_fails_fast():
    """ceiling_ha<=2 等于禁掉生长 → 也要快速失败。"""
    rc, out = _run_sampler_with_schedule(
        {"target_ha": 0, "slack": 4.0, "relax": 0.35, "tighten": 0.0,
         "warmup_frac": 0.25, "floor_steps": 0.0, "floor_thr": -2.0, "ceiling_ha": 2.0},
        max_steps=6)
    assert rc != 0, "应当以非零码退出"
    assert "ceiling_ha" in out, out[-300:]
    print("PASS S5 ceiling_ha 过小快速失败（rc=%d）" % rc)


if __name__ == "__main__":
    test_off_is_identity()
    test_warmup_no_intervention()
    test_relax_when_behind()
    test_tighten_only_when_enabled()
    test_absolute_size_thermostat()
    test_real_trajectory_never_relaxes_when_on_target()
    test_floor_s2_hard_size_floor()
    test_floor_takes_precedence_over_thermostat()
    test_floor_ge_max_steps_fails_fast()
    test_tiny_ceiling_fails_fast()
    print("ALL PASS (10/10)")
