#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pocket2Mol 口袋分子生成 GUI（仅 Python 标准库 Tkinter，零额外依赖）。

启动方式（在仓库根目录下）:
    python src\\gui\\app.py
（GUI 本身也在 Pocket2Mol conda 环境的 python 上运行，即 src/paths.py 的
PYTHON；双击运行需保证 .py 关联到该解释器，推荐用上面的命令行方式启动。）

功能:
    * 内置 4 个靶点下拉选择（自动填 PDB 路径与口袋中心）
    * 自定义 PDB 文件选择 + 口袋中心 (x,y,z) 手工输入
    * 采样参数配置: num_samples / beam_size / max_steps / seed /
      guided(enabled + lam / diversity_w)；高级项: device / bbox_size / 输出目录
    * 一键运行后端 sample_for_pdb.py 并实时显示进度日志
    * 结束后浏览 SMILES 列表、导出选中分子 SDF、打开输出目录

运行机制（重要约定）:
    * 后端命令:
        sample_for_pdb.py --pdb_path <pdb> --center " x,y,z" --config <yml>
                          --outdir <dir> [--bbox_size 23.0] [--device cuda]
      其中 --center 的值首字符必须是一个空格（项目记忆约定）。
    * 点击开始后基于 configs/sample_for_pdb.yml 模板动态生成临时 yml
      （UTF-8 无 BOM、LF 换行），写入 %TEMP%\\Pocket2MolGUI\\；
      后端会把该配置拷贝一份进会话目录存档。
    * 子进程 stdout/stderr 重定向到 outputs/gui_run_<时间戳>.log 文件句柄。
      禁止使用 stdout=PIPE：本机沙箱禁止命名管道（会 EPERM）。
    * 后台线程每 1 秒增量读取日志文件新行，解析 "[Pool] ... Finished N"
      更新进度计数。
    * 进程退出后扫描 outdir 下最新的 sample_for_pdb_* 会话目录，
      读取 SMILES.txt 列入结果区。SMILES 第 N 行对应 SDF/(N-1).sdf。
    * 停止按钮先向子进程组发送 CTRL_BREAK（后端捕获 KeyboardInterrupt
      后会保存已生成分子），10 秒后仍存活才强制 terminate()。
"""

import os
import re
import sys
import time
import queue
import shutil
import signal
import tempfile
import threading
import subprocess
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

# ==================== 路径配置（唯一来源: src/paths.py） ====================
# 本文件位于 src/gui/, 先把它上层的 src/ 放进 sys.path, 再取统一路径常量。
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, SRC, CONFIGS, OUTPUTS, TARGETS as TARGETS_DIR, PYTHON  # noqa: E402

PYTHON_EXE = PYTHON                                    # 后端解释器(可用 EONMOL_PYTHON 覆盖)
BACKEND_SCRIPT = os.path.join(SRC, 'sample_for_pdb.py')  # 迁移前在 Pocket2Mol\ 根下
CONFIG_TEMPLATE = os.path.join(CONFIGS, 'sample_for_pdb.yml')
DEFAULT_OUTDIR = OUTPUTS
PDB_DIR = TARGETS_DIR                                  # 迁移前的"靶点结构"目录

# 内置靶点: label 用于下拉框; code 用于在 PDB_DIR 中匹配 "<code>*.pdb" 文件;
# center 为口袋中心 x,y,z。confirmed=False 的中心是按参考配体（链 A）质心
# 估算的（估算方法已用 4EIY/2RH1 两个已知中心交叉验证，误差 < 0.1 A），
# 界面上标注 "[中心待确认]"。
BUILTIN_TARGETS = [
    {'label': 'A2A 腺苷受体 (4EIY)', 'code': '4EIY',
     'center': ('-0.4', '8.5', '17.1'), 'confirmed': True},
    {'label': 'β2-AR 肾上腺素受体 (2RH1)', 'code': '2RH1',
     'center': ('-29.5', '9.2', '6.9'), 'confirmed': True},
    {'label': 'D3 多巴胺受体 (3PBL) [中心待确认]', 'code': '3PBL',
     'center': ('0.1', '-14.8', '10.4'), 'confirmed': False},
    {'label': '5-HT2B 受体 (4IB4) [中心待确认]', 'code': '4IB4',
     'center': ('22.5', '18.3', '11.7'), 'confirmed': False},
]

SESSION_PREFIX = 'sample_for_pdb'   # 后端 get_new_log_dir 的会话目录前缀
POOL_FINISHED_RE = re.compile(r'\[Pool\].*?Finished\s+(\d+)')


# ==================== 纯函数工具（不依赖 Tk，便于自测） ====================

def resolve_builtin_pdb(target):
    """在 PDB_DIR 中按 PDB 代码前缀匹配靶点文件；找不到返回 None。"""
    if not os.path.isdir(PDB_DIR):
        return None
    code = target['code'].lower()
    for name in sorted(os.listdir(PDB_DIR)):
        low = name.lower()
        if low.startswith(code) and low.endswith('.pdb'):
            return os.path.join(PDB_DIR, name)
    return None


def load_template_text(path=CONFIG_TEMPLATE):
    """读取 yml 模板文本（utf-8-sig 以容错可能存在的 BOM）。"""
    with open(path, 'r', encoding='utf-8-sig') as f:
        return f.read()


def build_config_text(template_text, params):
    """基于模板文本生成新的 yml 文本。

    params 键: num_samples, beam_size, max_steps, seed, guided(bool),
               lam, div_w
    优先使用环境自带的 PyYAML（后端依赖，随 conda 环境提供）；
    不可用时退化为逐行改写（纯标准库），结果等价。
    """
    try:
        import yaml  # noqa: PLC0415  conda 环境自带
    except ImportError:
        yaml = None

    if yaml is not None:
        cfg = yaml.safe_load(template_text) if template_text.strip() else {}
        if not isinstance(cfg, dict):
            cfg = {}
        sample = cfg.setdefault('sample', {})
        if not isinstance(sample, dict):
            sample = cfg['sample'] = {}
        sample['seed'] = params['seed']
        sample['num_samples'] = params['num_samples']
        sample['beam_size'] = params['beam_size']
        sample['max_steps'] = params['max_steps']
        if params['guided']:
            # 后端 utils/guidance 引导束搜索: 化学分数(QED/SA) + 多样性惩罚,
            # lam 控制 logp 与化学分数的平衡; qed_w/sa_w 未在 GUI 暴露, 取默认 1.0
            sample['guided'] = {
                'enabled': True,
                'lam': params['lam'],
                'diversity_w': params['div_w'],
                'qed_w': 1.0,
                'sa_w': 1.0,
            }
        else:
            sample.pop('guided', None)
        return yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False,
                              default_flow_style=False)
    return patch_template_lines(template_text, params)


def patch_template_lines(text, params):
    """无 PyYAML 时的兜底实现: 逐行改写 sample 参数并注入 guided 块。"""
    values = {
        'seed': params['seed'],
        'num_samples': params['num_samples'],
        'beam_size': params['beam_size'],
        'max_steps': params['max_steps'],
    }
    lines = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
    out = []
    in_sample = False
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        stripped = line.strip()
        indent = len(line) - len(line.lstrip(' '))
        key = ''
        if stripped and not stripped.startswith('#') and ':' in stripped:
            key = stripped.split(':', 1)[0].strip()
        if key and indent == 0:
            in_sample = (key == 'sample')
        elif key and in_sample and indent > 0:
            if key in values:
                line = line[:indent] + '%s: %s' % (key, values[key])
            elif key == 'guided':
                # 跳过已有 guided 块（缩进比 guided 行更深的所有后续行）
                i += 1
                while i < n:
                    s2 = lines[i].strip()
                    ind2 = len(lines[i]) - len(lines[i].lstrip(' '))
                    if s2 and ind2 <= indent:
                        break
                    i += 1
                continue
        out.append(line)
        i += 1

    if params['guided']:
        block = ['  guided:',
                 '    enabled: true',
                 '    lam: %s' % params['lam'],
                 '    diversity_w: %s' % params['div_w'],
                 '    qed_w: 1.0',
                 '    sa_w: 1.0']
        si = None
        for idx, ln in enumerate(out):
            if ln.strip() == 'sample:':
                si = idx
                break
        if si is None:
            out.extend(['sample:'] + block)
        else:
            j = si + 1
            while j < len(out):
                s = out[j].strip()
                ind = len(out[j]) - len(out[j].lstrip(' '))
                if s and ind == 0:
                    break
                j += 1
            out[j:j] = block
    return '\n'.join(out)


def child_creation_flags():
    """子进程创建标志。

    * CREATE_NEW_PROCESS_GROUP: 允许用 CTRL_BREAK 优雅停止
      （后端捕获 KeyboardInterrupt 后保存已生成分子）。
    * CREATE_NO_WINDOW: 仅当父进程没有控制台（pythonw 启动）时使用，
      避免闪出黑窗；有控制台时不加它，以便 CTRL_BREAK 能送达同一控制台。
    """
    flags = 0
    if hasattr(subprocess, 'CREATE_NEW_PROCESS_GROUP'):
        flags |= subprocess.CREATE_NEW_PROCESS_GROUP
    if hasattr(subprocess, 'CREATE_NO_WINDOW'):
        try:
            import ctypes  # noqa: PLC0415  仅标准库
            if not ctypes.windll.kernel32.GetConsoleWindow():
                flags |= subprocess.CREATE_NO_WINDOW
        except Exception:
            pass
    return flags


# ==================== GUI 主体 ====================

class Pocket2MolGUI(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.title('Pocket2Mol 口袋分子生成 GUI')
        self.geometry('900x760')
        self.minsize(800, 640)

        # ---- 运行状态 ----
        self.proc = None            # subprocess.Popen
        self.running = False
        self.stop_requested = False
        self.kill_timer = None      # threading.Timer: 优雅停止超时后强杀
        self.ui_queue = queue.Queue()   # 后台线程 -> UI 线程 消息通道
        self.line_buf = ''          # 日志半行缓冲
        self.finished_count = 0
        self.run_start_time = 0.0
        self.outdir = None
        self.session_dir = None
        self.smiles_list = []

        self._build_vars()
        self._build_ui()
        self._wire_events()
        self.after(150, self._poll_ui)
        self.protocol('WM_DELETE_WINDOW', self._on_close)

        self._gui_log('就绪。工作目录: %s' % ROOT)
        self._gui_log('后端解释器: %s' % PYTHON_EXE)
        for t in BUILTIN_TARGETS:
            if resolve_builtin_pdb(t) is None:
                self._gui_log('[警告] 未找到靶点 %s 的 PDB 文件 (%s*.pdb @ %s)，'
                              '使用前请用"浏览..."手动选择。' % (t['label'], t['code'], PDB_DIR))
        if BUILTIN_TARGETS:
            self.combo_target.set(BUILTIN_TARGETS[0]['label'])
            self._on_target_selected(initial=True)

    # ---------- 变量 ----------
    def _build_vars(self):
        self.var_pdb = tk.StringVar(value='')
        self.var_cx = tk.StringVar(value='-0.4')
        self.var_cy = tk.StringVar(value='8.5')
        self.var_cz = tk.StringVar(value='17.1')
        self.var_num = tk.StringVar(value='100')      # 默认值取自模板
        self.var_beam = tk.StringVar(value='300')
        self.var_steps = tk.StringVar(value='50')
        self.var_seed = tk.StringVar(value='2020')
        self.var_guided = tk.BooleanVar(value=False)
        self.var_lam = tk.StringVar(value='1.0')
        self.var_div = tk.StringVar(value='0.0')
        self.var_device = tk.StringVar(value='cuda')
        self.var_bbox = tk.StringVar(value='23.0')
        self.var_outdir = tk.StringVar(value=DEFAULT_OUTDIR)
        self.var_status = tk.StringVar(value='Finished: 0')
        self.var_session = tk.StringVar(value='会话目录: (运行结束后显示)')

    # ---------- 界面 ----------
    def _build_ui(self):
        self.columnconfigure(0, weight=1)
        # 行: 0 靶点 / 1 参数 / 2 按钮 / 3 进度(可伸展) / 4 结果(可伸展)
        self.rowconfigure(3, weight=3)
        self.rowconfigure(4, weight=2)

        # == 靶点与口袋 ==
        fr_top = ttk.LabelFrame(self, text=' 靶点与口袋 ')
        fr_top.grid(row=0, column=0, sticky='ew', padx=8, pady=(8, 4))
        fr_top.columnconfigure(1, weight=1)

        ttk.Label(fr_top, text='靶点:').grid(row=0, column=0, sticky='w', padx=6, pady=4)
        self.combo_target = ttk.Combobox(fr_top, state='readonly', width=46,
                                         values=[t['label'] for t in BUILTIN_TARGETS])
        self.combo_target.grid(row=0, column=1, sticky='w', padx=4, pady=4)
        ttk.Label(fr_top, text='选中后自动填入 PDB 路径与中心')\
            .grid(row=0, column=2, sticky='w', padx=6)

        ttk.Label(fr_top, text='PDB 文件:').grid(row=1, column=0, sticky='w', padx=6, pady=4)
        self.ent_pdb = ttk.Entry(fr_top, textvariable=self.var_pdb)
        self.ent_pdb.grid(row=1, column=1, sticky='ew', padx=4)
        self.btn_pdb = ttk.Button(fr_top, text='浏览...', width=10, command=self._browse_pdb)
        self.btn_pdb.grid(row=1, column=2, sticky='w', padx=6, pady=4)

        ttk.Label(fr_top, text='口袋中心 x,y,z:').grid(row=2, column=0, sticky='w', padx=6, pady=4)
        fr_center = ttk.Frame(fr_top)
        fr_center.grid(row=2, column=1, sticky='w', padx=4)
        self.ent_cx = ttk.Entry(fr_center, textvariable=self.var_cx, width=10)
        self.ent_cy = ttk.Entry(fr_center, textvariable=self.var_cy, width=10)
        self.ent_cz = ttk.Entry(fr_center, textvariable=self.var_cz, width=10)
        self.ent_cx.grid(row=0, column=0, padx=(0, 4))
        self.ent_cy.grid(row=0, column=1, padx=(0, 4))
        self.ent_cz.grid(row=0, column=2, padx=(0, 8))
        ttk.Label(fr_center, text='(自定义 PDB 时手工填写; 启动命令的 --center 值首字符空格由程序自动添加)')\
            .grid(row=0, column=3, sticky='w')

        # == 采样参数 ==
        fr_param = ttk.LabelFrame(self, text=' 采样参数 ')
        fr_param.grid(row=1, column=0, sticky='ew', padx=8, pady=4)

        def _plabel(parent, text, col, row):
            ttk.Label(parent, text=text).grid(row=row, column=col, sticky='e', padx=(10, 3), pady=3)

        _plabel(fr_param, 'num_samples:', 0, 0)
        self.ent_num = ttk.Entry(fr_param, textvariable=self.var_num, width=8)
        self.ent_num.grid(row=0, column=1, sticky='w')
        _plabel(fr_param, 'beam_size:', 2, 0)
        self.ent_beam = ttk.Entry(fr_param, textvariable=self.var_beam, width=8)
        self.ent_beam.grid(row=0, column=3, sticky='w')
        _plabel(fr_param, 'max_steps:', 4, 0)
        self.ent_steps = ttk.Entry(fr_param, textvariable=self.var_steps, width=8)
        self.ent_steps.grid(row=0, column=5, sticky='w')
        _plabel(fr_param, 'seed:', 6, 0)
        self.ent_seed = ttk.Entry(fr_param, textvariable=self.var_seed, width=8)
        self.ent_seed.grid(row=0, column=7, sticky='w')

        self.chk_guided = tk.Checkbutton(fr_param, text='启用 guided 引导',
                                         variable=self.var_guided,
                                         command=self._toggle_guided)
        self.chk_guided.grid(row=1, column=0, columnspan=2, sticky='w', padx=(10, 3), pady=3)
        ttk.Label(fr_param, text='lam:').grid(row=1, column=2, sticky='e', padx=(6, 3))
        self.ent_lam = ttk.Entry(fr_param, textvariable=self.var_lam, width=8)
        self.ent_lam.grid(row=1, column=3, sticky='w')
        ttk.Label(fr_param, text='diversity_w:').grid(row=1, column=4, sticky='e', padx=(6, 3))
        self.ent_div = ttk.Entry(fr_param, textvariable=self.var_div, width=8)
        self.ent_div.grid(row=1, column=5, sticky='w')

        ttk.Label(fr_param, text='设备 device:').grid(row=1, column=6, sticky='e', padx=(10, 3))
        self.combo_device = ttk.Combobox(fr_param, textvariable=self.var_device,
                                         state='readonly', width=6, values=['cuda', 'cpu'])
        self.combo_device.grid(row=1, column=7, sticky='w')

        _plabel(fr_param, 'bbox_size:', 0, 2)
        self.ent_bbox = ttk.Entry(fr_param, textvariable=self.var_bbox, width=8)
        self.ent_bbox.grid(row=2, column=1, sticky='w')
        ttk.Label(fr_param, text='输出目录:').grid(row=2, column=2, sticky='e', padx=(10, 3), pady=3)
        self.ent_outdir = ttk.Entry(fr_param, textvariable=self.var_outdir)
        self.ent_outdir.grid(row=2, column=3, columnspan=4, sticky='ew', pady=3)
        self.btn_outdir = ttk.Button(fr_param, text='浏览...', width=10,
                                     command=self._browse_outdir)
        self.btn_outdir.grid(row=2, column=7, sticky='w', padx=(4, 0), pady=3)

        # == 控制(开始/停止) ==
        fr_btn = ttk.Frame(self)
        fr_btn.grid(row=2, column=0, sticky='ew', padx=8, pady=2)
        fr_btn.columnconfigure(0, weight=1)
        self.btn_start = ttk.Button(fr_btn, text='开始采样', command=self.start_sampling)
        self.btn_start.grid(row=0, column=1, sticky='e', padx=4)
        self.btn_stop = ttk.Button(fr_btn, text='停止采样', command=self.stop_sampling,
                                   state='disabled')
        self.btn_stop.grid(row=0, column=2, sticky='e', padx=4)
        ttk.Label(fr_btn, textvariable=self.var_status).grid(row=0, column=0, sticky='w', padx=6)

        # == 进度 ==
        fr_prog = ttk.LabelFrame(self, text=' 运行日志 / 进度 (每秒刷新) ')
        fr_prog.grid(row=3, column=0, sticky='nsew', padx=8, pady=4)
        fr_prog.rowconfigure(0, weight=1)
        fr_prog.columnconfigure(0, weight=1)
        self.txt_log = scrolledtext.ScrolledText(fr_prog, height=12, state='disabled',
                                                 font=('Consolas', 9))
        self.txt_log.grid(row=0, column=0, sticky='nsew')

        # == 结果 ==
        fr_result = ttk.LabelFrame(self, text=' 结果 (SMILES) ')
        fr_result.grid(row=4, column=0, sticky='nsew', padx=8, pady=(4, 8))
        fr_result.rowconfigure(1, weight=1)
        fr_result.columnconfigure(0, weight=1)
        ttk.Label(fr_result, textvariable=self.var_session, wraplength=860,
                  justify='left', foreground='#555555')\
            .grid(row=0, column=0, columnspan=2, sticky='ew', padx=6, pady=(4, 2))

        fr_list = ttk.Frame(fr_result)
        fr_list.grid(row=1, column=0, columnspan=2, sticky='nsew', padx=6)
        fr_list.rowconfigure(0, weight=1)
        fr_list.columnconfigure(0, weight=1)
        self.list_result = tk.Listbox(fr_list, exportselection=False,
                                      font=('Consolas', 9))
        self.list_result.grid(row=0, column=0, sticky='nsew')
        sb = ttk.Scrollbar(fr_list, orient='vertical', command=self.list_result.yview)
        sb.grid(row=0, column=1, sticky='ns')
        self.list_result.configure(yscrollcommand=sb.set)

        fr_res_btn = ttk.Frame(fr_result)
        fr_res_btn.grid(row=2, column=0, columnspan=2, sticky='ew', padx=6, pady=4)
        self.btn_export = ttk.Button(fr_res_btn, text='导出选中 SDF', command=self._export_sdf)
        self.btn_export.grid(row=0, column=0, sticky='w', padx=4)
        self.btn_open = ttk.Button(fr_res_btn, text='打开输出目录', command=self._open_outdir)
        self.btn_open.grid(row=0, column=1, sticky='w', padx=4)
        ttk.Label(fr_res_btn, text='(第 N 行对应会话目录 SDF/(N-1).sdf; 双击行复制 SMILES)')\
            .grid(row=0, column=2, sticky='w', padx=10)

        self._toggle_guided()

    def _wire_events(self):
        self.combo_target.bind('<<ComboboxSelected>>',
                               lambda _e: self._on_target_selected())
        self.list_result.bind('<Double-Button-1>', lambda _e: self._copy_smiles())

    # ---------- 日志（线程安全: 一律经 ui_queue 由 UI 线程写入） ----------
    def _gui_log(self, msg):
        self.ui_queue.put(('gui', msg))

    def _append_text(self, text):
        self.txt_log.configure(state='normal')
        self.txt_log.insert('end', text)
        self.txt_log.see('end')
        self.txt_log.configure(state='disabled')

    def _handle_lines(self, lines):
        for line in lines:
            if '[Pool]' in line:
                m = POOL_FINISHED_RE.search(line)
                if m:
                    self.finished_count = int(m.group(1))
                    self.var_status.set('Finished: %d' % self.finished_count)
        self._append_text('\n'.join(lines) + '\n')

    def _feed_raw(self, chunk):
        self.line_buf += chunk.replace('\r\n', '\n').replace('\r', '\n')
        if '\n' not in self.line_buf:
            return
        lines = self.line_buf.split('\n')
        self.line_buf = lines.pop()   # 最后可能是半行, 留给下一轮
        self._handle_lines(lines)

    def _poll_ui(self):
        try:
            while True:
                try:
                    kind, payload = self.ui_queue.get_nowait()
                except queue.Empty:
                    break
                if kind == 'gui':
                    self._append_text('[%s] %s\n' % (time.strftime('%H:%M:%S'), payload))
                elif kind == 'raw':
                    self._feed_raw(payload)
                elif kind == 'exit':
                    self._finish_run(payload)
        finally:
            try:
                self.after(150, self._poll_ui)
            except tk.TclError:
                pass

    # ---------- 交互 ----------
    def _on_target_selected(self, initial=False):
        label = self.combo_target.get()
        for t in BUILTIN_TARGETS:
            if t['label'] != label:
                continue
            path = resolve_builtin_pdb(t)
            if path is None:
                msg = ('在 %s 下未找到以 %s 开头的 .pdb 文件，\n'
                       '请使用 "浏览..." 手动选择 PDB 文件。' % (PDB_DIR, t['code']))
                self._gui_log('[警告] %s' % msg)
                if not initial:
                    messagebox.showwarning('未找到 PDB', msg)
                path = ''
            self.var_pdb.set(path)
            self.var_cx.set(t['center'][0])
            self.var_cy.set(t['center'][1])
            self.var_cz.set(t['center'][2])
            note = '' if t['confirmed'] else ' (中心待确认, 建议核对)'
            self._gui_log('载入靶点 %s -> PDB: %s | center: %s,%s,%s%s'
                          % (label, path or '<未找到>', *t['center'], note))
            break

    def _browse_pdb(self):
        initial = PDB_DIR if os.path.isdir(PDB_DIR) else ROOT
        path = filedialog.askopenfilename(
            title='选择 PDB 文件', initialdir=initial,
            filetypes=[('PDB 文件', '*.pdb'), ('所有文件', '*.*')])
        if path:
            self.var_pdb.set(os.path.normpath(path))
            self._gui_log('已选择自定义 PDB: %s' % path)

    def _browse_outdir(self):
        path = filedialog.askdirectory(title='选择输出目录',
                                       initialdir=self.var_outdir.get().strip() or ROOT)
        if path:
            self.var_outdir.set(os.path.normpath(path))

    def _toggle_guided(self):
        state = 'normal' if self.var_guided.get() else 'disabled'
        self.ent_lam.configure(state=state)
        self.ent_div.configure(state=state)

    def _copy_smiles(self):
        sel = self.list_result.curselection()
        if sel:
            self.clipboard_clear()
            self.clipboard_append(self.smiles_list[sel[0]])
            self.var_status.set('Finished: %d (SMILES 已复制到剪贴板)' % self.finished_count)

    # ---------- 参数校验 ----------
    def _validate_params(self):
        def as_int(var, name, minimum=None):
            raw = var.get().strip()
            try:
                v = int(raw)
            except ValueError:
                raise ValueError('%s 必须是整数 (当前输入: %r)' % (name, raw))
            if minimum is not None and v < minimum:
                raise ValueError('%s 必须不小于 %d (当前输入: %d)' % (name, minimum, v))
            return v

        def as_float(var, name):
            raw = var.get().strip()
            try:
                return float(raw)
            except ValueError:
                raise ValueError('%s 必须是数字 (当前输入: %r)' % (name, raw))

        try:
            pdb = self.var_pdb.get().strip()
            if not pdb:
                raise ValueError('请先选择靶点或用 "浏览..." 选择 PDB 文件。')
            if not os.path.isfile(pdb):
                raise ValueError('PDB 文件不存在:\n%s' % pdb)
            center = []
            for var, axis in ((self.var_cx, 'x'), (self.var_cy, 'y'), (self.var_cz, 'z')):
                raw = var.get().strip()
                try:
                    float(raw)
                except ValueError:
                    raise ValueError('口袋中心 %s 必须是数字 (当前输入: %r)' % (axis, raw))
                center.append(raw)
            params = {
                'pdb_path': pdb,
                'center': center,
                'num_samples': as_int(self.var_num, 'num_samples', 1),
                'beam_size': as_int(self.var_beam, 'beam_size', 1),
                'max_steps': as_int(self.var_steps, 'max_steps', 1),
                'seed': as_int(self.var_seed, 'seed'),
                'guided': bool(self.var_guided.get()),
                'lam': as_float(self.var_lam, 'lam'),
                'div_w': as_float(self.var_div, 'diversity_w'),
                'bbox_size': as_float(self.var_bbox, 'bbox_size'),
                'device': self.var_device.get().strip() or 'cuda',
                'outdir': self.var_outdir.get().strip() or DEFAULT_OUTDIR,
            }
            if params['bbox_size'] <= 0:
                raise ValueError('bbox_size 必须大于 0 (当前输入: %s)' % params['bbox_size'])
        except ValueError as e:
            messagebox.showerror('参数错误', str(e))
            return None

        checks = [
            (os.path.isdir(ROOT), '后端工作目录不存在: %s' % ROOT),
            (os.path.isfile(PYTHON_EXE), '后端解释器不存在: %s' % PYTHON_EXE),
            (os.path.isfile(BACKEND_SCRIPT), '后端脚本不存在: %s' % BACKEND_SCRIPT),
            (os.path.isfile(CONFIG_TEMPLATE), '配置模板不存在: %s' % CONFIG_TEMPLATE),
        ]
        for ok, msg in checks:
            if not ok:
                messagebox.showerror('环境错误', msg)
                return None
        return params

    # ---------- 运行 ----------
    def _write_config_yml(self, params, ts):
        """基于模板生成临时 yml (UTF-8 无 BOM, LF 换行)。"""
        tmpdir = os.path.join(tempfile.gettempdir(), 'Pocket2MolGUI')
        os.makedirs(tmpdir, exist_ok=True)
        yml_path = os.path.join(tmpdir, 'gui_config_%s.yml' % ts)
        text = build_config_text(load_template_text(), params)
        with open(yml_path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(text)
        return yml_path

    def start_sampling(self):
        if self.running:
            messagebox.showwarning('正在运行', '采样进程仍在运行，请先点击 "停止采样"。')
            return
        params = self._validate_params()
        if params is None:
            return
        try:
            outdir = os.path.abspath(params['outdir'])
            os.makedirs(outdir, exist_ok=True)
            ts = time.strftime('%Y%m%d_%H%M%S')
            yml_path = self._write_config_yml(params, ts)
            log_path = os.path.join(outdir, 'gui_run_%s.log' % ts)

            # 项目记忆约定: --center 值首字符必须是一个空格
            center_str = ' ' + ','.join(params['center'])
            cmd = [
                PYTHON_EXE, BACKEND_SCRIPT,
                '--pdb_path', params['pdb_path'],
                '--center', center_str,
                '--config', yml_path,
                '--outdir', outdir,
                '--bbox_size', str(params['bbox_size']),
                '--device', params['device'],
            ]
            env = dict(os.environ)
            env['PYTHONIOENCODING'] = 'utf-8'   # 子进程日志统一 UTF-8, 便于解析
            env['PYTHONUTF8'] = '1'

            self._gui_log('生成临时配置: %s' % yml_path)
            self._gui_log('启动命令: %s' % subprocess.list2cmdline(cmd))
            self._gui_log('运行日志: %s' % log_path)

            with open(log_path, 'w', encoding='utf-8', errors='replace') as fh:
                # 禁止 stdout=PIPE: 本机沙箱禁止命名管道 (EPERM)，必须重定向到文件
                self.proc = subprocess.Popen(
                    cmd, cwd=ROOT,
                    stdout=fh, stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL,
                    env=env, creationflags=child_creation_flags())
            # 子进程已继承句柄, 父进程立即关闭自身副本

            self.running = True
            self.stop_requested = False
            self.finished_count = 0
            self.var_status.set('Finished: 0')
            self.run_start_time = time.time() - 5.0   # 会话目录 mtime 比较 5s 容差
            self.outdir = outdir
            self.session_dir = None
            self.smiles_list = []
            self.list_result.delete(0, 'end')
            self.var_session.set('会话目录: (运行结束后显示)')
            self.btn_start.configure(state='disabled')
            self.btn_stop.configure(state='normal')
            self._append_text('\n' + '=' * 78 + '\n')
            self._gui_log('采样已启动 (PID %s)。' % self.proc.pid)
            threading.Thread(target=self._reader_thread, args=(log_path,),
                             daemon=True).start()
        except Exception as e:
            self._gui_log('启动失败: %r' % (e,))
            messagebox.showerror('启动失败', '无法启动后端进程:\n%r' % (e,))

    def _reader_thread(self, log_path):
        """后台线程: 每秒增量读取日志文件; 进程退出后发 ('exit', code)。"""
        proc = self.proc
        offset = 0
        try:
            with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
                while True:
                    try:
                        f.seek(offset)
                        chunk = f.read()
                        offset = f.tell()
                    except OSError:
                        chunk = ''
                    if chunk:
                        self.ui_queue.put(('raw', chunk))
                    if proc.poll() is not None:
                        time.sleep(1.0)   # 收尾: 等子进程缓冲区落盘后再读一次
                        try:
                            f.seek(offset)
                            tail = f.read()
                        except OSError:
                            tail = ''
                        if tail:
                            self.ui_queue.put(('raw', tail))
                        break
                    time.sleep(1.0)
        except Exception as e:
            self.ui_queue.put(('gui', '读取日志线程异常: %r' % (e,)))
        finally:
            self.ui_queue.put(('exit', proc.poll() if proc is not None else -1))

    def stop_sampling(self):
        proc = self.proc
        if proc is None or proc.poll() is not None:
            messagebox.showinfo('提示', '当前没有正在运行的采样进程。')
            return
        self.stop_requested = True
        self.btn_stop.configure(state='disabled')
        self._gui_log('请求停止: 先发送 CTRL_BREAK (后端捕获后可保存已生成分子), '
                      '10 秒后仍未退出则强制 terminate...')
        try:
            os.kill(proc.pid, signal.CTRL_BREAK_EVENT)
        except Exception as e:
            self._gui_log('CTRL_BREAK 发送失败 (%r)，直接 terminate。' % (e,))
            try:
                proc.terminate()
            except Exception:
                pass
        if self.kill_timer is not None:
            self.kill_timer.cancel()
        self.kill_timer = threading.Timer(10.0, self._hard_kill_if_alive)
        self.kill_timer.daemon = True
        self.kill_timer.start()

    def _hard_kill_if_alive(self):
        proc = self.proc
        if proc is not None and proc.poll() is None:
            self.ui_queue.put(('gui', '进程 10 秒后仍未退出, 强制 terminate()。'))
            try:
                proc.terminate()
            except Exception as e:
                self.ui_queue.put(('gui', 'terminate 失败: %r' % (e,)))

    def _finish_run(self, returncode):
        if self.kill_timer is not None:
            self.kill_timer.cancel()
            self.kill_timer = None
        self.running = False
        self.btn_start.configure(state='normal')
        self.btn_stop.configure(state='disabled')
        if returncode == 0:
            self._gui_log('后端进程正常退出 (返回码 0)。')
        else:
            self._gui_log('后端进程退出, 返回码 %s (非 0 通常表示被终止或异常)。' % returncode)
        if self.line_buf:
            self._handle_lines([self.line_buf])
            self.line_buf = ''

        self._gui_log('扫描输出目录: %s' % self.outdir)
        session = self._find_latest_session(self.outdir, self.run_start_time) \
            if self.outdir else None
        if session is None:
            self._gui_log('[警告] 未在输出目录找到 sample_for_pdb_* 会话目录, '
                          '请检查上方日志中的报错信息。')
            self.var_session.set('会话目录: 未找到')
            return
        self.session_dir = session
        self._gui_log('会话目录: %s' % session)
        self._load_results(session)

    def _find_latest_session(self, outdir, since_ts):
        """取 outdir 下 sample_for_pdb_* 会话目录中最新者。

        注意: 只接受本次运行之后 (>= since_ts) 产生的目录。
        旧实现 `return best or best_any` 会回退到**任意最新的历史会话**,
        在后端启动失败(配置非法/解释器报错)时把上一次的 SMILES 当成本次结果显示。
        """
        if not outdir or not os.path.isdir(outdir):
            return None
        best, best_mtime = None, -1.0
        for name in os.listdir(outdir):
            path = os.path.join(outdir, name)
            if not name.startswith(SESSION_PREFIX) or not os.path.isdir(path):
                continue
            try:
                mtime = os.path.getmtime(path)
            except OSError:
                continue
            if mtime >= since_ts and mtime > best_mtime:
                best, best_mtime = path, mtime
        return best

    def _load_results(self, session):
        self.list_result.delete(0, 'end')
        self.smiles_list = []
        smiles_path = os.path.join(session, 'SMILES.txt')
        if os.path.isfile(smiles_path):
            with open(smiles_path, 'r', encoding='utf-8', errors='replace') as f:
                self.smiles_list = [s.strip() for s in f if s.strip()]
        if self.smiles_list:
            for i, s in enumerate(self.smiles_list):
                self.list_result.insert('end', '%d. %s' % (i + 1, s))
            self.var_session.set('会话目录: %s | 共 %d 个分子 (第 N 行对应 SDF/(N-1).sdf)'
                                 % (session, len(self.smiles_list)))
            self._gui_log('已载入 %d 个 SMILES。' % len(self.smiles_list))
        else:
            self.var_session.set('会话目录: %s | 未找到 SMILES.txt '
                                 '(可能未生成完整分子, 或进程被强制终止)' % session)
            self._gui_log('[警告] 未找到 SMILES.txt: %s' % smiles_path)

    # ---------- 结果导出 ----------
    def _export_sdf(self):
        if not self.session_dir:
            messagebox.showinfo('提示', '还没有可导出的会话结果，请先运行一次采样。')
            return
        sel = self.list_result.curselection()
        if not sel:
            messagebox.showinfo('提示', '请先在结果列表中选择一个分子。')
            return
        idx = sel[0]
        sdf_src = os.path.join(self.session_dir, 'SDF', '%d.sdf' % idx)
        if not os.path.isfile(sdf_src):
            messagebox.showerror(
                '导出失败',
                '未找到:\n%s\n\n(SMILES 第 %d 行应对应 SDF/%d.sdf，\n'
                '该文件可能未被后端生成。)' % (sdf_src, idx + 1, idx))
            return
        dst = filedialog.asksaveasfilename(
            title='导出 SDF', defaultextension='.sdf',
            initialfile='mol_%d.sdf' % (idx + 1),
            filetypes=[('SDF 文件', '*.sdf'), ('所有文件', '*.*')])
        if not dst:
            return
        try:
            shutil.copyfile(sdf_src, dst)
        except OSError as e:
            messagebox.showerror('导出失败', '复制 SDF 失败:\n%r' % (e,))
            return
        self._gui_log('已导出 SDF: %s' % dst)
        messagebox.showinfo('导出成功', "\n".join(
            ['已导出:', dst, '', 'SMILES:', self.smiles_list[idx]]))

    def _open_outdir(self):
        target = self.session_dir
        if not (target and os.path.isdir(target)):
            target = self.outdir
        if not (target and os.path.isdir(target)):
            target = self.var_outdir.get().strip()
        if not target or not os.path.isdir(target):
            messagebox.showwarning('提示', '输出目录不存在:\n%s' % target)
            return
        startfile = getattr(os, 'startfile', None)
        if startfile is None:
            self._gui_log('[警告] 当前平台不支持 os.startfile，请手动打开: %s' % target)
            return
        try:
            startfile(target)
        except Exception as e:
            self._gui_log('打开目录失败: %r' % (e,))

    # ---------- 关闭 ----------
    def _on_close(self):
        if self.running and self.proc is not None and self.proc.poll() is None:
            if not messagebox.askokcancel(
                    '退出', '采样仍在运行，退出将强制终止后端进程（结果不会保存）。\n确定退出？'):
                return
            try:
                self.proc.terminate()
            except Exception:
                pass
        self.destroy()


def main():
    if sys.platform == 'win32':
        try:   # 高分屏下界面更清晰 (标准库 ctypes, 失败不影响运行)
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
    app = Pocket2MolGUI()
    app.mainloop()


if __name__ == '__main__':
    main()
