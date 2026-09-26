# Pocket2Mol 采样 GUI 使用说明

## v2.0 Pro 版（推荐）

面向"外行友好 + 可直接分发"的新版界面，源码 `src/gui/app_pro.py`，打包产物
**`src\gui\dist\7-eonmol.exe`**（单文件、无控制台、Windows x64、9.9 MB）。

**启动方式**

- **快速启动版（推荐日常使用）**：双击桌面快捷方式 `7-eonmol (快速启动).lnk`，
  或 `src\gui\dist\fast\7-eonmol\7-eonmol.exe` —— 免解压目录版，实测启动 1.3–2.7 秒
- **单文件版（推荐分发）**：`src\gui\dist\7-eonmol.exe` —— 单文件免安装，实测启动约 2.6 秒
- 源码运行（仓库根目录下）：`python src\gui\app_pro.py [--tab 0-3] [--theme dark]`
  （`python` 用 Pocket2Mol 环境的解释器，即 `src/paths.py` 的 `PYTHON`）

**启动耗时优化记录**（原 onefile 版冷启动 22 秒 → 2.6 秒）

| 措施 | 效果 |
|---|---|
| 环境自检改为后台线程执行（原实现同步等待 `import torch` 子进程 5–15 秒，把窗口显示一起卡住） | 主要收益 |
| 排除无关模块（numpy/torch/rdkit/PIL/ssl/hashlib/sqlite3/setuptools 等） | 单文件 9.9→7.7 MB |
| 裁剪 Tcl/Tk 冗余数据（tzdata 时区、msgs 多语言、多余编码，共约 800 个文件） | 快速版 983→171 个文件 |
| 额外提供 onedir 快速启动版（无解压步骤） | 再次启动 1.3 秒 |

> 首次在新机器上运行或刚复制完成时，Windows Defender 会扫描新文件，可能额外增加数秒；
> 把程序目录加入杀软白名单可消除该影响。

**界面结构（4 个标签页 + 顶栏按钮）**

| 标签页 | 作用 |
|---|---|
| 生成分子 | 选靶点卡片 → 开始生成 → 进度/状态 → 候选表排序筛选 → 2D 结构预览 → 导出 SDF/CSV |
| 候选库浏览 | 读取 `results/`（`paths.RESULTS`）下任意 CSV（26 个数据文件），排序、看结构、导出视图 |
| 分析结果 | 汇总 `results/`、`docs/` 与输出目录下的报告与 CSV（md/pptx/pdf/对接 summary），双击用系统程序打开并预览前 60 行 |
| 环境设置 | 修改仓库/解释器/靶点目录等路径（写 `gui_config.json`，立即生效）+ 一键环境自检 |

顶栏右侧：`?` 打开使用说明窗口，`☾/☀` 切换浅色 / 深色主题，环境状态灯常驻显示。
主流程刻意保持最简：**选靶点 → 点「开始生成分子」→ 看结果**，界面上不出现分步引导文字。

**相对旧版 (app_easy.py / app.py) 的改进**

- 现代观感：圆角自绘按钮、卡片式靶点选择、分区标题、树表隔行色、Canvas 圆角图标
- **浅色 / 深色双主题**，右上角一键切换（含全部 ttk 样式与自绘控件）
- 高 DPI 适配：源码运行按 DPI 放大窗口并同步 Tk 缩放；打包版交由系统统一缩放，避免字体二次放大
- 自绘控件宽度用字体实测值自适应（`ui_kit.measure`），任何 DPI/字体下都不裁字
- 结果表：点击表头排序、类药分阈值滑块、关键字搜索、统计摘要（候选数 / QED 中位）、导出 CSV
- 靶点卡片直接显示 PDB 代码、口袋中心与适应症；界面无"第 N 步"式引导文字
- 帮助改为顶栏 `?` 弹窗（原先占用一个主标签页，属于层级误用）
- `--tab`/`--theme` 启动参数
- 关闭窗口时仍在采样会二次确认并终止子进程（不产生孤儿进程占 GPU）

**打包（可复现）**

```powershell
powershell -ExecutionPolicy Bypass -File src\gui\build_exe.ps1
```

脚本依次：生成图标 → 归档旧 exe 到 `src\gui\dist\legacy\` → PyInstaller onefile 打包
（`--add-binary tcl86t.dll/tk86t.dll`，conda Tk 8.6）→ 校验产物大小。
解释器取自 `EONMOL_PYTHON`（缺省用当前激活 conda 环境的 `python.exe`），
`Tcl/Tk` 取该环境的 `Library\bin`。
`.ps1` 必须带 UTF-8 BOM，否则 Windows PowerShell 5.1 按 GBK 读会语法报错。

**开发自检**

```powershell
python -m py_compile src\gui\app_pro.py src\gui\ui_kit.py
```

---

## 旧版说明（`src/gui/app.py`，保留）

基于 **Python 标准库 Tkinter** 的本地图形界面（零额外依赖），用于驱动后端
`sample_for_pdb.py` 完成口袋分子生成：选靶点 → 配参数 → 跑采样 → 看进度 →
浏览 SMILES → 导出 SDF。

- GUI 文件：`src/gui/app.py`
- 运行解释器：Pocket2Mol 环境的 python（`src/paths.py` 的 `PYTHON`，
  可用环境变量 `EONMOL_PYTHON` 覆盖）（GUI 与后端共用该环境）
- 后端工作目录：仓库根（`src/paths.py` 的 `ROOT`，GUI 启动子进程时自动设为 `cwd`）

## 一、启动方法

在命令行（cmd / PowerShell）执行：

```bat
python src\gui\app.py
```

工作目录任意均可；GUI 内部固定以仓库根（`paths.ROOT`）作为后端工作目录。

## 二、界面与操作流程

1. **靶点与口袋**
   - `靶点` 下拉框：内置 4 个靶点，选中即自动填入 PDB 路径与口袋中心：

     | 靶点 | PDB 文件（`data/targets/`，即 `paths.TARGETS`） | 口袋中心 (x,y,z) | 备注 |
     |---|---|---|---|
     | A2A 腺苷受体 | `4EIY_A2A受体.pdb` | `-0.4, 8.5, 17.1` | 已确认 |
     | β2-AR | `2RH1_β2肾上腺素受体.pdb` | `-29.5, 9.2, 6.9` | 已确认 |
     | D3 多巴胺受体 | `3PBL_D3多巴胺受体.pdb` | `0.1, -14.8, 10.4` | **待确认** |
     | 5-HT2B | `4IB4_5HT2B受体.pdb` | `22.5, 18.3, 11.7` | **待确认** |

     **关于两个“待确认”中心**：4EIY / 2RH1 的已知中心经核实恰为其参考配体
     （ZMA / CAU，链 A）的原子质心；用同一方法对 3PBL 的配体 ETQ（链 A，残基
     1200）与 4IB4 的配体 ERM（链 A，残基 2001）计算质心得到上表数值（分别约
     `0.09, -14.83, 10.43` 与 `22.45, 18.28, 11.73`）。方法一致性已验证，但数值
     仍建议跑样前人工核对；如需修改，选中靶点后直接改中心三个输入框即可。
   - `PDB 文件` + `浏览...`：自定义 PDB；选中后手工填写中心 x/y/z 三个输入框。
     程序启动命令中的 `--center` 值首字符空格由 GUI 自动添加（项目记忆约定，
     无需手动输入）。
2. **采样参数**（写入临时 yml 的 `sample:` 段）
   - `num_samples`（目标分子数，默认 100）、`beam_size`（束宽，300）、
     `max_steps`（最大步数，50）、`seed`（随机种子，2020）——默认值取自模板。
   - `guided 引导`：勾选后写 `sample.guided: {enabled: true, lam, diversity_w,
     qed_w: 1.0, sa_w: 1.0}`，启用后端 `utils/guidance` 的引导束搜索（QED/SA
     化学分数 + Tanimoto 多样性惩罚）；`lam` 平衡 logp 与化学分数，
     `diversity_w` 控制重复结构惩罚力度。未勾选则不写 guided 段（等价关闭）。
   - 高级项：`device`（cuda / cpu，默认 cuda）、`bbox_size`（口袋盒子边长，
     默认 23.0）、`输出目录`（默认 `outputs/`，即 `paths.OUTPUTS`）。
3. **开始 / 停止**
   - `开始采样`：校验参数 → 生成临时 yml → 启动后端 → 实时刷日志。
     运行期间按钮禁用，**重复点击被禁止**；参数非法（数字格式、文件不存在等）
     会弹窗指明具体字段，不会启动。
   - `停止采样`：先向子进程发送 **CTRL_BREAK**（后端捕获 KeyboardInterrupt
     后会照常保存已生成分子到 SMILES.txt / SDF）；若 10 秒后仍存活则强制
     `terminate()`。被强制终止的进程**不会**保存部分结果。
4. **运行日志 / 进度**
   - 后台线程每秒增量读取日志文件，完整回显后端输出；
   - 自动解析 `[Pool] Queue N | Finished M | Failed K` 行，左下角
     `Finished: M` 实时更新已完成分子数。
5. **结果区**
   - 进程退出后自动扫描输出目录中最新的 `sample_for_pdb_*` 会话目录，
     读取 `SMILES.txt` 列表展示；**第 N 行对应会话目录内 `SDF/(N-1).sdf`**。
   - `导出选中 SDF`：选中一行 → 另存为 .sdf（复制 `SDF/(N-1).sdf`）。
   - `打开输出目录`：用资源管理器打开当前会话目录（无会话时打开输出目录）。
   - 双击结果行可复制该条 SMILES 到剪贴板。

## 三、运行机制（ GUI 内部实现约定）

- **临时配置**：点击开始后读取 `configs/sample_for_pdb.yml` 模板，改写
  `sample.seed / num_samples / beam_size / max_steps` 并按需注入 `guided` 段，
  写为 `%TEMP%\Pocket2MolGUI\gui_config_<时间戳>.yml`（**UTF-8 无 BOM、LF 换行**）。
  优先用环境自带的 PyYAML 解析改写；PyYAML 不可用时退化为逐行改写（仍是纯标准库）。
  后端会把该配置拷贝一份进会话目录存档。
- **启动命令**（子进程，`cwd=` 仓库根 `paths.ROOT`）：

  ```bat
  python src\sample_for_pdb.py
      --pdb_path <pdb> --center " x,y,z" --config <临时yml> --outdir <输出目录>
      --bbox_size 23.0 --device cuda
  ```

- **日志重定向**：子进程 stdout/stderr 合并重定向到
  `outputs/gui_run_<时间戳>.log` 的文件句柄。**禁止 `stdout=PIPE`**——本机沙箱
  禁止命名管道，用 PIPE 会 EPERM；读取由后台线程按文件增量完成。
- **健壮性**：所有数字参数转换失败均弹窗提示具体字段；运行中异常（启动失败、
  日志线程异常、导出失败等）会显示在进度框或弹窗；关闭窗口时若仍在采样会先
  确认并终止子进程。
- 路径来源唯一：`app.py` 顶部不再写死盘符，`PYTHON_EXE / BACKEND_SCRIPT /
  CONFIG_TEMPLATE / DEFAULT_OUTDIR / PDB_DIR` 全部由 `src/paths.py` 推导；
  需要换位置时用环境变量（`EONMOL_PYTHON` 等），不要改源码。

## 四、输出结构

```
outputs/
├── gui_run_20250101_120000.log          ← GUI 子进程运行日志（stdout/stderr 合并）
└── sample_for_pdb_2025_01_01__12_00_03/ ← 后端会话目录（每此一次采样一个）
    ├── log.txt                          ← 后端 logger 日志
    ├── gui_config_….yml 的副本           ← 本次使用的采样配置（自动拷入）
    ├── <pdb 文件名>.pdb                  ← 本次使用的口袋（自动拷入）
    ├── SMILES.txt                       ← 每行一个 SMILES（结果区来源）
    ├── SDF/0.sdf, 1.sdf, …              ← 与 SMILES.txt 行序一一对应（0 基）
    └── samples_*.pt                     ← 每步束搜索池快照
```

## 五、常见问题

- **弹窗“参数错误”**：检查对应输入框是否为合法整数/浮点数（num_samples /
  beam_size / max_steps / seed / lam / diversity_w / bbox_size / 中心 xyz）。
- **提示“未找到 PDB”**：`data/targets/`（`paths.TARGETS`）下没有以该 PDB 代码开头的
  `.pdb` 文件（或目录不存在），用“浏览...”手动选择即可。
- **CUDA 不可用**：把高级项 `device` 切为 `cpu`（会显著变慢）。
- **进程刚退出但结果区为空**：多为被强制 terminate（无保存）、未生成任何完整
  分子（`Failed/Duplicate` 太多）、或口袋中心落在蛋白质外导致
  `No atoms found in the bounding box`——查看进度框里的后端报错。
- **日志出现少量乱码**：GUI 已强制子进程 `PYTHONIOENCODING=utf-8`，若手工改动
  请保持；乱码不影响 `[Pool]` 进度解析（该行为 ASCII）。
- **语法自检**（不启动 GUI）：

  ```bat
  python -m py_compile src\gui\app.py
  ```

