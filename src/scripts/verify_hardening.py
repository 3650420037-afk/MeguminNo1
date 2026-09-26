# -*- coding: utf-8 -*-
"""加固项验证: BOM读取 / 芳香键告警 / kekulize工具 / 编译"""
import os, sys, io, warnings, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import PYTHON, ROOT, src_on_path
src_on_path()
os.chdir(ROOT)
ok = []
bad = []

# 1) BOM 读取加固: 写出带 BOM 的 yml, load_config 应能正确读到 model 键
from utils.misc import load_config
tmp = os.path.join(tempfile.gettempdir(), "bom_test.yml")
with open(tmp, "wb") as f:
    f.write("\ufeff".encode("utf-8") + b"model:\n    checkpoint: ./ckpt/pretrained_Pocket2Mol.pt\nsample:\n  seed: 1\n")
try:
    c = load_config(tmp)
    if c.model.checkpoint.endswith("pretrained_Pocket2Mol.pt"):
        ok.append("BOM 读取加固: load_config 正确解析带 BOM 配置 (config.model.checkpoint 可达)")
    else:
        bad.append("BOM 读取加固: 解析结果异常 %r" % c)
except Exception as e:
    bad.append("BOM 读取加固: 仍失败 %s" % e)

# 2) 芳香键告警: 构造含键值 4 的最小 SDF 文本, parse_sdf_file 应告警
from utils.protein_ligand import parse_sdf_file
sdf_arom = """benzene
  test

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0
  2  3  4  0
  3  4  4  0
  4  5  4  0
  5  6  4  0
  6  1  4  0
M  END
"""
ptmp = os.path.join(tempfile.gettempdir(), "arom_test.sdf")
open(ptmp, "w", encoding="utf-8").write(sdf_arom)
try:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        d = parse_sdf_file(ptmp)
        hit = any("芳香键" in str(x.message) for x in w)
    if hit:
        ok.append("芳香键告警: 含键值4的 SDF 正确触发告警 (数据防御可见)")
    else:
        bad.append("芳香键告警: 未触发告警 (键类型 %s)" % set(d["bond_type"].tolist()))
except Exception as e:
    bad.append("芳香键告警: 异常 %s" % e)

# 3) 洁净数据不误报: 用真实数据集跑 kekulize 工具干跑
r = subprocess.run([PYTHON, "src/scripts/kekulize_dataset.py", "--dataset", "data/gpcr_a2a"],
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
out = r.stdout.decode("utf-8", errors="ignore")
if "结论: 数据洁净" in out:
    ok.append("kekulize 工具: 真实数据集干跑报'洁净'(无芳香键, 不误报)")
else:
    bad.append("kekulize 工具: 干跑输出异常 -> %s" % out.strip().splitlines()[-1:])

# 4) 编译全部改动文件
files = ["src/utils/misc.py", "src/utils/protein_ligand.py", "src/models/common.py",
         "src/models/encoders/__init__.py", "src/scripts/kekulize_dataset.py",
         "src/scripts/build_gpcr_dataset.py", "src/gui/app_easy.py"]
r = subprocess.run([PYTHON, "-m", "py_compile"] + files, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
if r.returncode == 0:
    ok.append("编译: 7 个改动文件全部通过")
else:
    bad.append("编译失败: %s" % r.stdout.decode("utf-8", errors="ignore"))

print("=" * 60)
for x in ok:
    print("[PASS]", x)
for x in bad:
    print("[FAIL]", x)
print("=" * 60)
print("结果:", "ALL PASS" if not bad else "%d 项失败" % len(bad))
