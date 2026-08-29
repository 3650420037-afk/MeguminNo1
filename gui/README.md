# Pocket2Mol 采样 GUI 使用说明

基于 **Python 标准库 Tkinter** 的本地图形界面（零额外依赖），用于驱动后端
`sample_for_pdb.py` 完成口袋分子生成：选靶点 → 配参数 → 跑采样 → 看进度 →
浏览 SMILES → 导出 SDF。

- GUI 文件：`gui/app.py`
- 运行解释器：`D:\Miniconda3\envs\Pocket2Mol\python.exe`（GUI 与后端共用该环境）
- 后端工作目录：`D:\MMModel\Pocket2Mol`（GUI 启动子进程时自动设置）

## 一、启动方法

在命令行（cmd / PowerShell）执行：

```bat
D:\Miniconda3\envs\Pocket2Mol\python.exe D:\MMModel\Pocket2Mol\gui\app.py
```

工作目录任意均可；GUI 内部固定以 `D:\MMModel\Pocket2Mol` 作为后端工作目录。

## 二、界面与操作流程

1. **靶点与口袋**
   - `靶点` 下拉框：内置 4 个靶点，选中即自动填入 PDB 路径与口袋中心：

     | 靶点 | PDB 文件（`D:\MMModel\靶点结构\`） | 口袋中心 (x,y,z) | 备注 |
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
     默认 23.0）、`输出目录`（默认 `D:\MMModel\Pocket2Mol\outputs`）。
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
- **启动命令**（子进程，`cwd=D:\MMModel\Pocket2Mol`）：

  ```bat
  D:\Miniconda3\envs\Pocket2Mol\python.exe sample_for_pdb.py
      --pdb_path <pdb> --center " x,y,z" --config <临时yml> --outdir <输出目录>
      --bbox_size 23.0 --device cuda
  ```

- **日志重定向**：子进程 stdout/stderr 合并重定向到
  `outputs/gui_run_<时间戳>.log` 的文件句柄。**禁止 `stdout=PIPE`**——本机沙箱
  禁止命名管道，用 PIPE 会 EPERM；读取由后台线程按文件增量完成。
- **健壮性**：所有数字参数转换失败均弹窗提示具体字段；运行中异常（启动失败、
  日志线程异常、导出失败等）会显示在进度框或弹窗；关闭窗口时若仍在采样会先
  确认并终止子进程。
- 修改 GUI 端路径常量（解释器 / 工作目录 / PDB 目录等）：编辑 `app.py` 顶部
  `WORK_DIR / PYTHON_EXE / PDB_DIR` 等常量。

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
- **提示“未找到 PDB”**：`D:\MMModel\靶点结构\` 下没有以该 PDB 代码开头的
  `.pdb` 文件（或目录不存在），用“浏览...”手动选择即可。
- **CUDA 不可用**：把高级项 `device` 切为 `cpu`（会显著变慢）。
- **进程刚退出但结果区为空**：多为被强制 terminate（无保存）、未生成任何完整
  分子（`Failed/Duplicate` 太多）、或口袋中心落在蛋白质外导致
  `No atoms found in the bounding box`——查看进度框里的后端报错。
- **日志出现少量乱码**：GUI 已强制子进程 `PYTHONIOENCODING=utf-8`，若手工改动
  请保持；乱码不影响 `[Pool]` 进度解析（该行为 ASCII）。
- **语法自检**（不启动 GUI）：

  ```bat
  D:\Miniconda3\envs\Pocket2Mol\python.exe -m py_compile D:\MMModel\Pocket2Mol\gui\app.py
  ```
