# -*- coding: utf-8 -*-
"""7-eonmol · GPCR 靶向分子生成器 (打包用 GUI)

相对旧版 (app_easy.py) 的改进:
  * 现代观感: 圆角控件/卡片/双主题(浅色·深色)切换/自适应 DPI/图标
  * 五个标签页: 生成分子 · 候选库 · 分析结果 · 环境设置 · 使用说明
  * 结果表: 可点表头排序、QED 阈值筛选、关键字搜索、统计摘要、导出 CSV
  * 内置环境自检与路径设置 (写 gui_config.json), 分发到其他机器免改源码

后端计算全部通过子进程调用 Pocket2Mol 环境, 本文件仅依赖 Python 标准库 + Tkinter。
"""
import os, sys, json, csv, re, glob, threading, subprocess, time, base64
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter import font as tkfont

# ---- 路径前置: 本文件在 src/gui/, 先放 gui/ (找 ui_kit), 再放它上层的 src/ ----
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, OUTPUTS, RESULTS, DOCS, TARGETS as TARGETS_DIR, PYTHON  # noqa: E402

import ui_kit
from ui_kit import (THEMES, FONT, FONT_MONO, RoundedButton, TargetCard, SectionTitle, Card,
                     app_icon, measure)

APP_TITLE = "7-eonmol"
APP_SUB = "GPCR 靶向分子生成器"
APP_VER = "v2.0 Pro"
WIN_TITLE = "%s · %s  %s" % (APP_TITLE, APP_SUB, APP_VER)
_HERE = os.path.dirname(os.path.abspath(__file__))
_EXE_DIR = os.path.dirname(sys.executable) if getattr(sys, "frozen", False) else _HERE
# 工程目录(仓库根) —— 默认取自 paths.ROOT, 部署/分发时由下面两个分支覆盖
REPO = ROOT
PY = PYTHON
PDB_DIR = TARGETS_DIR
OUTROOT = OUTPUTS
CFG = os.path.join(_EXE_DIR, "gui_config.json")
LIB_DIR = RESULTS
# 分析结果/报告的来源目录(缺省 docs/, 可用 gui_config.json 的 report_dirs 追加)
REPORT_DIRS = [DOCS]

# 部署模式: exe 同目录存在 repo\ 子文件夹时自动改用包内仓库与环境
if os.path.isdir(os.path.join(_EXE_DIR, "repo")):
    REPO = os.path.join(_EXE_DIR, "repo")
    OUTROOT = os.path.join(REPO, "outputs")
    LIB_DIR = os.path.join(REPO, "library")
    # 靶点结构目录逐个候选探测: 取第一个真的含 .pdb 的目录
    # (旧实现硬编码 structures\, 而仓库里结构其实放在 targets\ —— 部署后四个靶点全部找不到)
    PDB_DIR = os.path.join(REPO, "structures")
    for _d in ("structures", "targets", "pdb", os.path.join("data", "structures")):
        _cand = os.path.join(REPO, _d)
        if glob.glob(os.path.join(_cand, "*.pdb")):
            PDB_DIR = _cand
            break
    for _c in (os.path.join(_EXE_DIR, "env", "python.exe"),
               os.path.join(_EXE_DIR, "env", "Scripts", "python.exe")):
        if os.path.exists(_c):
            PY = _c
            break

try:
    _c = json.load(open(CFG, encoding="utf-8-sig"))
    REPO = _c.get("repo", REPO); PY = _c.get("python", PY)
    PDB_DIR = _c.get("pdb_dir", PDB_DIR); OUTROOT = _c.get("outroot", OUTROOT)
    LIB_DIR = _c.get("library", LIB_DIR)
    REPORT_DIRS = list(_c.get("report_dirs", REPORT_DIRS))
except Exception:
    pass

TARGETS = [
    ("A2A 腺苷受体", "4EIY", " -0.4,8.5,17.1", "帕金森病 · 肿瘤免疫 · 炎症"),
    ("D3 多巴胺受体", "3PBL", " 0.085,-14.828,10.432", "精神分裂症 · 成瘾"),
    ("5-HT2B 血清素受体", "4IB4", " 22.448,18.284,11.726", "偏头痛 · 肺动脉高压"),
]
FLAG = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0


def enable_dpi():
    """Windows 高 DPI 处理.

    源码运行 (python.exe 自带 per-monitor manifest): 声明感知并同步 Tk 缩放 -> 物理像素渲染,
    窗口尺寸按 DPI 放大以保持视觉大小一致。
    打包 exe (PyInstaller manifest 未声明): 交由 Windows 统一缩放 (逻辑像素布局),
    此时不改 Tk scaling, 否则字体会被二次放大。
    返回 DPI 缩放系数, 0.0 表示"不要改 Tk 缩放"。
    """
    if os.name != "nt":
        return 1.0
    if getattr(sys, "frozen", False):
        return 0.0
    try:
        import ctypes
        try:
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            ctypes.windll.user32.SetProcessDPIAware()
        return ctypes.windll.user32.GetDpiForSystem() / 96.0
    except Exception:
        return 1.0


class App:
    def __init__(self, root, theme="light"):
        self.root = root
        self.dpi_scale = getattr(root, "_p2m_scale", 1.0)
        self.theme_name = theme
        self.T = THEMES[theme]
        self.proc = None
        self.log_path = None
        self.log_pos = 0
        self.session_dir = None
        self.lib_dir = None
        self.rows = []
        self.sort_col = None
        self.sort_desc = True
        self._themed = []          # [(widget, {option: palette_key})]
        self._custom = []          # [(callable(T))] 主题切换回调
        self._target_cards = {}
        self.sel_target = TARGETS[0][0]
        self._lib_rows = []
        self._lib_cols = []
        self._build()
        self.set_theme(self.theme_name)

    # ================= 主题 =================
    def reg(self, w, **opts):
        """登记 tk 组件的配色 (option=调色板键), 主题切换时自动重刷"""
        self._themed.append((w, opts))
        return w

    def on_theme(self, fn):
        self._custom.append(fn)

    def set_theme(self, name):
        T = THEMES[name]
        self.theme_name = name
        self.T = T
        ui_kit.apply_theme(self.root, self.style, T)
        for w, opts in self._themed:
            try:
                w.configure(**{k: T[v] for k, v in opts.items()})
            except Exception:
                pass
        for w in self._target_cards.values():
            w.set_theme(T, T["bg"])
        for w in getattr(self, "_all_buttons", []):
            w.set_theme(T)
        if hasattr(self, "tree"):
            self.tree.tag_configure("odd", background=T["tree_alt"])
            self.tree.tag_configure("even", background=T["card"])
        if hasattr(self, "lib_tree"):
            self.lib_tree.tag_configure("odd", background=T["tree_alt"])
            self.lib_tree.tag_configure("even", background=T["card"])
        for fn in self._custom:
            try:
                fn(T)
            except Exception:
                pass
        self.btn_theme.set_text("☾ 深色" if name == "light" else "☀ 浅色")

    def toggle_theme(self):
        self.set_theme("dark" if self.theme_name == "light" else "light")

    # ================= 小部件助手 =================
    def lab(self, parent, text="", size=10, bold=False, bg="bg", fg="fg", **kw):
        w = tk.Label(parent, text=text, font=(FONT, size, "bold" if bold else "normal"),
                     anchor=kw.pop("anchor", "w"), **kw)
        self.reg(w, bg=bg, fg=fg)
        return w

    def btn(self, parent, text, cmd, kind="primary", w=150, h=38, icon="", size=11, bg="bg"):
        b = RoundedButton(parent, text=text, command=cmd, kind=kind, width=w, height=h,
                          icon=icon, font=(FONT, size, "bold"),
                          bg=self.T[bg] if isinstance(bg, str) else self.T["bg"])
        b._pal_key = bg
        self._all_buttons.append(b)
        return b

    def card(self, parent, pad=14):
        c = Card(parent, pad=pad)
        self._custom.append(lambda T, c=c: c.set_theme(T))
        return c

    def section(self, parent, text, sub=None):
        """分区标题 (自动跟随主题切换)"""
        s = SectionTitle(parent, text, theme=self.T, sub=sub)
        self._custom.append(lambda T, s=s: s.set_theme(T))
        return s

    # ================= 骨架 =================
    def _build(self):
        self._all_buttons = []
        self.style = ttk.Style(self.root)
        self.root.title(WIN_TITLE)
        sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        k = max(1.0, self.dpi_scale)
        w = int(min(1280 * k, sw - 60 * k))
        h = int(min(880 * k, sh - 130 * k))
        self.root.geometry("%dx%d+%d+%d" % (w, h, max(0, (sw - w) // 2), max(0, (sh - h) // 3)))
        self.root.minsize(int(1060 * k), int(680 * k))
        app_icon(self.root)

        # ---- 顶栏 ----
        head = tk.Frame(self.root, height=78)
        head.pack(fill="x", side="top")
        head.pack_propagate(False)
        self.reg(head, bg="header")
        self.hbar = tk.Frame(head, width=5)
        self.hbar.pack(side="left", fill="y")
        self.reg(self.hbar, bg="accent")
        left = tk.Frame(head)
        left.pack(side="left", padx=16)
        self.reg(left, bg="header")
        t1 = tk.Label(left, text=APP_TITLE, font=(FONT, 14, "bold"))
        t1.pack(anchor="w", pady=(2, 0))
        self.reg(t1, bg="header", fg="fg")
        t2 = tk.Label(left, text="GPCR 靶向分子生成器 · 引导束搜索 · 类药性筛选",
                      font=(FONT, 9))
        t2.pack(anchor="w")
        self.reg(t2, bg="header", fg="sub")
        right = tk.Frame(head)
        right.pack(side="right", padx=14)
        self.reg(right, bg="header")
        self.lbl_env = tk.Label(right, text="● 环境检测中…", font=(FONT, 9))
        self.lbl_env.pack(side="left", padx=10)
        self.reg(self.lbl_env, bg="header", fg="warn")
        self.btn_theme = self.btn(right, "☾ 深色", self.toggle_theme, kind="ghost", w=100, h=34,
                                  icon="", size=10, bg="header")
        self.btn_theme.pack(side="left", padx=6)
        self.btn_help = self.btn(right, "?", self.show_help, kind="ghost", w=44, h=34,
                                 icon="", size=12, bg="header")
        self.btn_help.pack(side="left")

        # ---- 底部状态栏 (先占位, 避免被 Notebook 挤掉) ----
        self._status = tk.StringVar(value="就绪。选择靶点后点击「开始生成」。")
        bar = tk.Frame(self.root)
        bar.pack(fill="x", side="bottom")
        self.reg(bar, bg="card2")
        sl = tk.Label(bar, textvariable=self._status, font=(FONT, 9), anchor="w")
        sl.pack(fill="x", padx=12, pady=5)
        self.reg(sl, bg="card2", fg="sub")

        # ---- 标签页 ----
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill="both", expand=True, padx=12, pady=(8, 12))
        self.tab_gen = tk.Frame(self.nb)
        self.tab_lib = tk.Frame(self.nb)
        self.tab_rep = tk.Frame(self.nb)
        self.tab_set = tk.Frame(self.nb)
        for f in (self.tab_gen, self.tab_lib, self.tab_rep, self.tab_set):
            self.reg(f, bg="bg")
        self.nb.add(self.tab_gen, text="  生成分子  ")
        self.nb.add(self.tab_lib, text="  候选库浏览  ")
        self.nb.add(self.tab_rep, text="  分析结果  ")
        self.nb.add(self.tab_set, text="  环境设置  ")

        self._build_generate(self.tab_gen)
        self._build_library(self.tab_lib)
        self._build_reports(self.tab_rep)
        self._build_settings(self.tab_set)

    # ================= Tab1 生成分子 =================
    def _build_generate(self, P):
        top = tk.Frame(P)
        top.pack(fill="x", padx=4, pady=(6, 0))
        self.reg(top, bg="bg")

        # 左: 靶点卡片
        lf = tk.Frame(top)
        lf.pack(side="left", fill="both", expand=True)
        self.reg(lf, bg="bg")
        st = self.section(lf, "选择靶点", sub="点击卡片切换")
        st.pack(anchor="w", pady=(4, 8))
        self._custom.append(lambda T, s=st: s.set_theme(T))
        grid = tk.Frame(lf)
        grid.pack(anchor="w")
        self.reg(grid, bg="bg")
        for i, (name, pdb, center, disease) in enumerate(TARGETS):
            c = TargetCard(grid, name, pdb, center, disease, command=self._pick_target,
                           theme=self.T, width=250, height=110)
            c.grid(row=i // 2, column=i % 2, padx=(0, 12), pady=(0, 12))
            self._target_cards[name] = c
        self._target_cards[TARGETS[0][0]].set_selected(True)


        # 右: 任务面板 (固定宽度, 避免挤压左侧卡片)
        rtf = tk.Frame(top, width=452, height=250)
        self.rtf = rtf
        self._rtf_base_h = 250
        rtf.pack(side="left", fill="y", padx=(18, 0))
        rtf.pack_propagate(False)
        self.reg(rtf, bg="bg")
        rt = self.card(rtf, pad=16)
        rt.pack(fill="both", expand=True)
        R = rt.inner
        self.lab(R, "运行生成", size=12, bold=True, bg="card").pack(anchor="w")
        self.summary = self.lab(R, "", size=9, fg="sub", bg="card", justify="left")
        self.summary.pack(anchor="w", pady=(2, 10))
        self._refresh_summary()

        brow = tk.Frame(R)
        brow.pack(fill="x")
        self.reg(brow, bg="card")
        self.btn_start = self.btn(brow, "开始生成分子", self.start, kind="success", w=214, h=54,
                                  icon="▶", size=13, bg="card")
        self.btn_start.pack(side="left")
        col = tk.Frame(brow)
        col.pack(side="left", padx=10)
        self.reg(col, bg="card")
        self.btn_stop = self.btn(col, "停止", self.stop, kind="ghost", w=112, h=25, icon="■",
                                 size=9, bg="card")
        self.btn_stop.pack(pady=(1, 4))
        self.btn_open = self.btn(col, "打开结果目录", self.open_dir, kind="ghost", w=112, h=25,
                                 icon="", size=9, bg="card")
        self.btn_open.pack()
        self.btn_stop.set_state("disabled")
        self.btn_open.set_state("disabled")

        self.progress = ttk.Progressbar(R, mode="determinate", maximum=100)
        self.progress.pack(fill="x", pady=(14, 6))
        self.lbl_run = self.lab(R, "等待开始", size=10, fg="accent", bg="card")
        self.lbl_run.pack(anchor="w")

        # 高级参数
        advh = tk.Frame(R)
        advh.pack(fill="x", pady=(12, 0))
        self.reg(advh, bg="card")
        self.adv_on = tk.BooleanVar(value=False)
        cb = ttk.Checkbutton(advh, text="高级参数", variable=self.adv_on, command=self._toggle_adv)
        cb.pack(side="left")
        self.advfrm = tk.Frame(R)
        self.reg(self.advfrm, bg="card")
        self.adv_h = 0
        self.advvals = {}
        specs = [("num_samples", "50", "目标分子数"), ("beam_size", "100", "搜索宽度"),
                 ("max_steps", "50", "最大步数"), ("lam", "3.0", "引导强度 λ"),
                 ("diversity_w", "0.5", "多样性权重"), ("seed", "9000", "随机种子")]
        for i, (k, v, txt) in enumerate(specs):
            r, c = divmod(i, 2)
            self.lab(self.advfrm, txt, size=9, fg="sub", bg="card").grid(
                row=r * 2, column=c, sticky="w", padx=(0, 10), pady=(6, 0))
            e = ttk.Entry(self.advfrm, width=9)
            e.insert(0, v)
            e.grid(row=r * 2 + 1, column=c, sticky="w", padx=(0, 10))
            self.advvals[k] = e
        self._toggle_adv()
        # 任务面板宽度按内容自适应 (高 DPI 下字体更宽, 固定宽度会裁字)
        try:
            need = max(self.btn_start._cw + 10 + max(self.btn_stop._cw, self.btn_open._cw) + 44,
                       measure(self.advfrm, "多样性权重", FONT, 9) + 12 +
                       measure(self.advfrm, "9000", FONT_MONO, 9) + 76,
                       measure(self.summary, self.summary.cget("text"), FONT, 9) + 44)
            rtf.configure(width=int(need))
        except Exception:
            pass

        # 结果区
        mid = tk.Frame(P)
        mid.pack(fill="both", expand=True, padx=4, pady=(12, 0))
        self.reg(mid, bg="bg")
        mid.columnconfigure(0, weight=1)
        mid.rowconfigure(0, weight=1)
        lg = tk.Frame(mid)
        lg.grid(row=0, column=0, sticky="nsew")
        self.reg(lg, bg="bg")
        lg.columnconfigure(0, weight=1)
        lg.rowconfigure(2, weight=1)
        head = tk.Frame(lg)
        head.grid(row=0, column=0, sticky="ew")
        self.reg(head, bg="bg")
        self.section(head, "候选分子").pack(side="left")
        self.lbl_stat = self.lab(head, "尚无结果", size=9, fg="sub")
        self.lbl_stat.pack(side="left", padx=12)
        self.log_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(head, text="显示运行日志", variable=self.log_on,
                        command=self._toggle_log).pack(side="right")

        filt = tk.Frame(lg)
        filt.grid(row=1, column=0, sticky="ew", pady=(6, 4))
        self.reg(filt, bg="bg")
        self.lab(filt, "类药分 ≥", size=9, fg="sub").pack(side="left")
        self.qed_min = tk.DoubleVar(value=0.0)
        self.scale = ttk.Scale(filt, from_=0.0, to=1.0, orient="horizontal",
                               variable=self.qed_min, length=140, command=self._on_qed)
        self.scale.pack(side="left", padx=6)
        self.lbl_qed = self.lab(filt, "0.00", size=9, fg="accent")
        self.lbl_qed.pack(side="left")
        self.lab(filt, "搜索", size=9, fg="sub").pack(side="left", padx=(14, 4))
        self.q_text = tk.StringVar()
        ent = ttk.Entry(filt, textvariable=self.q_text, width=26)
        ent.pack(side="left")
        ent.bind("<KeyRelease>", lambda e: self._fill_table())

        tw = tk.Frame(lg)
        tw.grid(row=2, column=0, sticky="nsew")
        tw.columnconfigure(0, weight=1)
        tw.rowconfigure(0, weight=1)
        self.reg(tw, bg="bg")
        cols = ("rank", "qed", "mw", "logp", "sa", "scaf", "smiles")
        self.tree = ttk.Treeview(tw, columns=cols, show="headings", height=4)
        for c, w, txt in (("rank", 52, "序号"), ("qed", 118, "类药分 QED"), ("mw", 76, "分子量"),
                          ("logp", 62, "LogP"), ("sa", 54, "SA"), ("scaf", 110, "骨架"),
                          ("smiles", 330, "分子结构编码 (SMILES)")):
            self.tree.heading(c, text=txt, command=lambda cc=c: self._sort(cc))
            self.tree.column(c, width=max(w, measure(self.tree, txt, FONT, 9, True) + 36),
                             anchor="w")
        ys = ttk.Scrollbar(tw, orient="vertical", command=self.tree.yview)
        xs = ttk.Scrollbar(tw, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        ys.grid(row=0, column=1, sticky="ns")
        xs.grid(row=1, column=0, sticky="ew")
        self.tree.bind("<<TreeviewSelect>>", self._show_mol)

        bar = tk.Frame(lg)
        bar.grid(row=3, column=0, sticky="ew", pady=(6, 0))
        self.reg(bar, bg="bg")
        self.btn_sdf = self.btn(bar, "导出选中分子 3D 结构", self.export_sdf, kind="primary",
                                w=250, h=40, size=10)
        self.btn_sdf.pack(side="left")
        self.btn_copy = self.btn(bar, "复制分子编码", self.copy_smiles, kind="ghost", w=190,
                                 h=40, size=10)
        self.btn_copy.pack(side="left", padx=10)
        self.btn_exp = self.btn(bar, "导出当前列表为 CSV", self.export_csv, kind="ghost", w=210,
                                h=40, size=10)
        self.btn_exp.pack(side="left")

        # 右: 2D 预览 + 选中分子的操作
        rp = tk.Frame(mid, width=390)
        rp.grid(row=0, column=1, sticky="ns", padx=(14, 0))
        rp.grid_propagate(False)
        self.reg(rp, bg="bg")
        self.section(rp, "结构预览").pack(anchor="w", pady=(0, 6))
        self.mol_canvas = tk.Label(rp, text="在左侧点击一个分子\n这里显示 2D 结构式",
                                   width=34, height=11, relief="flat")
        self.mol_canvas.pack()
        self.reg(self.mol_canvas, bg="card", fg="sub")
        self.mol_img = None
        self.mol_info = self.lab(rp, "", size=9, fg="sub", justify="left", wraplength=360)
        self.mol_info.pack(anchor="w", pady=8)

        # 日志 (默认收起, 勾选后展开并挤压结果区)
        self.logfrm = tk.Frame(P)
        self.reg(self.logfrm, bg="bg")
        self.logtxt = tk.Text(self.logfrm, height=7, font=(FONT_MONO, 9), relief="flat",
                              wrap="none")
        self.logtxt.pack(fill="both", expand=True)
        self.reg(self.logtxt, bg="log_bg", fg="log_fg", insertbackground="log_fg")
        self.logtxt.configure(state="disabled")

    def _on_qed(self, *_):
        try:
            self.lbl_qed.configure(text="%.2f" % self.qed_min.get())
        except Exception:
            pass
        self._fill_table()

    def _toggle_adv(self):
        if self.adv_on.get():
            self.advfrm.pack(fill="x")
            self.root.update_idletasks()
            self.adv_h = self.advfrm.winfo_reqheight()
        else:
            self.advfrm.pack_forget()
            self.adv_h = 0
        try:
            self.rtf.configure(height=self._rtf_base_h + self.adv_h)
        except Exception:
            pass
        self._refresh_summary()

    def _toggle_log(self):
        if self.log_on.get():
            self.logfrm.pack(fill="x", padx=4, pady=(4, 8))
        else:
            self.logfrm.pack_forget()

    def _refresh_summary(self):
        try:
            v = {k: e.get() for k, e in self.advvals.items()}
            self.summary.configure(
                text="目标 %s 个 · 束宽 %s · 最大步数 %s\n引导 λ=%s · 多样性 %s · 预计 8-15 分钟"
                     % (v["num_samples"], v["beam_size"], v["max_steps"], v["lam"],
                        v["diversity_w"]))
        except Exception:
            pass

    def _pick_target(self, card):
        for name, c in self._target_cards.items():
            c.set_selected(c is card)
        self.sel_target = card.title

        self._status.set("已选靶点: %s（%s）  —  点「开始生成分子」即可" % (card.title, card.disease))

    # ================= Tab2 候选库 =================
    def _build_library(self, P):
        top = tk.Frame(P)
        top.pack(fill="x", padx=4, pady=(8, 6))
        self.reg(top, bg="bg")
        self.section(top, "已有化合物库",
                     sub="读取 %s 下的分析结果" % LIB_DIR).pack(side="left")
        row = tk.Frame(P)
        row.pack(fill="x", padx=4)
        self.reg(row, bg="bg")
        self.lab(row, "数据源", size=9, fg="sub").pack(side="left")
        self.lib_src = tk.StringVar()
        self.cmb = ttk.Combobox(row, textvariable=self.lib_src, width=34, state="readonly")
        self.cmb.pack(side="left", padx=6)
        self.cmb.bind("<<ComboboxSelected>>", lambda e: self.load_source())
        self.btn_refresh = self.btn(row, "刷新列表", self.refresh_sources, kind="ghost", w=104,
                                    h=30, size=9)
        self.btn_refresh.pack(side="left", padx=4)
        self.lib_info = self.lab(row, "", size=9, fg="sub")
        self.lib_info.pack(side="left", padx=10)
        self.btn_lib_exp = self.btn(row, "导出当前视图", self.lib_export, kind="ghost", w=124,
                                    h=30, size=9)
        self.btn_lib_exp.pack(side="right")

        mid = tk.Frame(P)
        mid.pack(fill="both", expand=True, padx=4, pady=8)
        self.reg(mid, bg="bg")
        mid.columnconfigure(0, weight=1)
        mid.rowconfigure(0, weight=1)
        tw = tk.Frame(mid)
        tw.grid(row=0, column=0, sticky="nsew")
        self.reg(tw, bg="bg")
        self.lib_tree = ttk.Treeview(tw, show="headings", height=16)
        ys = ttk.Scrollbar(tw, orient="vertical", command=self.lib_tree.yview)
        xs = ttk.Scrollbar(tw, orient="horizontal", command=self.lib_tree.xview)
        self.lib_tree.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        self.lib_tree.pack(side="top", fill="both", expand=True)
        xs.pack(side="bottom", fill="x")
        ys.pack(side="right", fill="y")
        self.lib_tree.bind("<<TreeviewSelect>>", self._lib_show_mol)
        rp = tk.Frame(mid, width=390)
        rp.grid(row=0, column=1, sticky="ns", padx=(14, 0))
        rp.grid_propagate(False)
        self.reg(rp, bg="bg")
        self.section(rp, "结构预览").pack(anchor="w", pady=(0, 6))
        self.lib_canvas = tk.Label(rp, text="点击左侧任意一行\n这里显示 2D 结构式", width=40,
                                   height=17)
        self.lib_canvas.pack()
        self.reg(self.lib_canvas, bg="card", fg="sub")
        self.lib_img = None
        self.lib_detail = self.lab(rp, "", size=9, fg="sub", justify="left", wraplength=360)
        self.lib_detail.pack(anchor="w", pady=8)
        self.refresh_sources()

    # ================= Tab3 报告 =================
    def _build_reports(self, P):
        top = tk.Frame(P)
        top.pack(fill="x", padx=4, pady=(8, 6))
        self.reg(top, bg="bg")
        self.section(top, "分析结果与文档",
                     sub="双击在系统默认程序中打开").pack(side="left")
        body = tk.Frame(P)
        body.pack(fill="both", expand=True, padx=4, pady=4)
        self.reg(body, bg="bg")
        lf = tk.Frame(body)
        lf.pack(side="left", fill="both", expand=True)
        self.reg(lf, bg="bg")
        self.rep_tree = ttk.Treeview(lf, columns=("kind", "size", "time"), show="headings",
                                     height=18)
        for c, w, t in (("kind", 300, "文件"), ("size", 90, "大小"), ("time", 140, "修改时间")):
            self.rep_tree.heading(c, text=t)
            self.rep_tree.column(c, width=max(w, measure(self.rep_tree, t, FONT, 9, True) + 36),
                                 anchor="w")
        ys = ttk.Scrollbar(lf, orient="vertical", command=self.rep_tree.yview)
        self.rep_tree.configure(yscrollcommand=ys.set)
        self.rep_tree.pack(side="left", fill="both", expand=True)
        ys.pack(side="left", fill="y")
        self.rep_tree.bind("<Double-1>", lambda e: self.open_report())
        rf = tk.Frame(body, width=440)
        rf.pack(side="left", fill="both", padx=(14, 0))
        rf.pack_propagate(False)
        self.reg(rf, bg="bg")
        self.rep_path = self.lab(rf, "未选择文件", size=9, fg="sub", wraplength=420)
        self.rep_path.pack(anchor="w")
        btns = tk.Frame(rf)
        btns.pack(fill="x", pady=6)
        self.reg(btns, bg="bg")
        self.btn_openrep = self.btn(btns, "打开文件", self.open_report, kind="primary", w=118,
                                    h=32, size=10)
        self.btn_openrep.pack(side="left")
        self.btn_openfolder = self.btn(btns, "打开所在文件夹", self.open_report_folder,
                                       kind="ghost", w=140, h=32, size=10)
        self.btn_openfolder.pack(side="left", padx=8)
        self.rep_preview = tk.Text(rf, height=22, font=(FONT_MONO, 8), relief="flat", wrap="none")
        self.rep_preview.pack(fill="both", expand=True)
        self.reg(self.rep_preview, bg="card", fg="fg", insertbackground="fg")
        self.rep_preview.configure(state="disabled")
        self.rep_tree.bind("<<TreeviewSelect>>", lambda e: self.preview_report())
        self.refresh_reports()
        self._custom.append(lambda T: self._tag_report_rows())

    def _tag_report_rows(self):
        for i, iid in enumerate(self.rep_tree.get_children()):
            self.rep_tree.item(iid, tags=("odd" if i % 2 else "even",))
        self.rep_tree.tag_configure("odd", background=self.T["tree_alt"])
        self.rep_tree.tag_configure("even", background=self.T["card"])

    # ================= Tab4 设置 =================
    def _build_settings(self, P):
        top = tk.Frame(P)
        top.pack(fill="x", padx=4, pady=(8, 6))
        self.reg(top, bg="bg")
        self.section(top, "路径与环境",
                     sub="分发到其他机器时在此修改并保存").pack(side="left")
        body = tk.Frame(P)
        body.pack(fill="both", expand=True, padx=4)
        self.reg(body, bg="bg")
        lf = tk.Frame(body)
        lf.pack(side="left", fill="both", expand=True)
        self.reg(lf, bg="bg")
        self.set_entries = {}
        for k, txt, dflt in (("repo", "Pocket2Mol 仓库目录", REPO),
                             ("python", "Python 解释器 (含 torch)", PY),
                             ("pdb_dir", "靶点 PDB 目录", PDB_DIR),
                             ("outroot", "输出目录", OUTROOT),
                             ("library", "化合物库目录", LIB_DIR)):
            self.lab(lf, txt, size=10, bold=True).pack(anchor="w", pady=(8, 2))
            e = ttk.Entry(lf)
            e.insert(0, dflt)
            e.pack(fill="x")
            self.set_entries[k] = e
        brow = tk.Frame(lf)
        brow.pack(fill="x", pady=14)
        self.reg(brow, bg="bg")
        self.btn_save = self.btn(brow, "保存设置", self.save_settings, kind="success", w=130,
                                 h=36, size=10)
        self.btn_save.pack(side="left")
        self.btn_check = self.btn(brow, "运行环境自检", self.run_env_check, kind="primary",
                                  w=140, h=36, size=10)
        self.btn_check.pack(side="left", padx=10)
        self.btn_cfg = self.btn(brow, "打开配置文件", self.open_cfg, kind="ghost", w=130, h=36,
                                size=10)
        self.btn_cfg.pack(side="left")

        rf = tk.Frame(body)
        rf.pack(side="left", fill="both", expand=True, padx=(16, 0))
        self.reg(rf, bg="bg")
        self.lab(rf, "自检结果", size=11, bold=True).pack(anchor="w", pady=(6, 4))
        self.chk_tree = ttk.Treeview(rf, columns=("v",), show="headings", height=18)
        self.chk_tree.heading("v", text="项目")
        self.chk_tree.column("v", width=760, anchor="w")
        ysc = ttk.Scrollbar(rf, orient="vertical", command=self.chk_tree.yview)
        xsc = ttk.Scrollbar(rf, orient="horizontal", command=self.chk_tree.xview)
        self.chk_tree.configure(yscrollcommand=ysc.set, xscrollcommand=xsc.set)
        self.chk_tree.pack(side="top", fill="both", expand=True)
        xsc.pack(side="bottom", fill="x")
        ysc.pack(side="right", fill="y")
        self.run_env_check()

    # ================= 帮助窗口 (顶栏 ? 按钮) =================
    def show_help(self):
        win = getattr(self, "_help_win", None)
        if win is not None and win.winfo_exists():
            win.lift()
            win.focus_force()
            return
        T = self.T
        win = tk.Toplevel(self.root)
        self._help_win = win
        win.title("使用说明 — %s · %s" % (APP_TITLE, APP_SUB))
        k = max(1.0, self.dpi_scale)
        w, h = int(860 * k), int(700 * k)
        try:
            self.root.update_idletasks()
            px = self.root.winfo_rootx() + max(0, (self.root.winfo_width() - w) // 2)
            py = self.root.winfo_rooty() + max(0, (self.root.winfo_height() - h) // 3)
        except Exception:
            px = py = 80
        win.geometry("%dx%d+%d+%d" % (w, h, px, py))
        win.transient(self.root)
        win.configure(bg=T["bg"])
        try:
            app_icon(win)
        except Exception:
            pass
        head = tk.Frame(win)
        head.pack(fill="x", padx=16, pady=(14, 4))
        self.reg(head, bg="bg")
        self.section(head, "使用说明").pack(side="left")
        txt = tk.Text(win, font=(FONT, 10), relief="flat", wrap="word", padx=18, pady=14)
        txt.pack(fill="both", expand=True, padx=16, pady=(4, 8))
        self.reg(txt, bg="card", fg="fg", insertbackground="fg")
        txt.insert("1.0", HELP_TEXT)
        txt.configure(state="disabled")
        btn = self.btn(win, "关闭", win.destroy, kind="primary", w=130, h=38, size=10, bg="bg")
        btn.pack(pady=(0, 16))
        win.protocol("WM_DELETE_WINDOW", win.destroy)

    # ================= 生成流程 =================
    def validate(self):
        tgt = next(t for t in TARGETS if t[0] == self.sel_target)
        pdb = glob.glob(os.path.join(PDB_DIR, tgt[1] + "*.pdb"))
        if not pdb:
            messagebox.showerror("找不到靶点文件", "在 %s 下未找到 %s*.pdb" % (PDB_DIR, tgt[1]))
            return None, None
        if not os.path.exists(PY):
            messagebox.showerror("环境错误", "Python 解释器不存在:\n" + PY)
            return None, None
        try:
            p = {k: float(e.get()) for k, e in self.advvals.items()}
        except ValueError:
            messagebox.showerror("参数错误", "高级参数必须是数字 (整数或小数)")
            return None, None
        if p["num_samples"] < 1 or p["beam_size"] < 1 or p["max_steps"] < 1:
            messagebox.showerror("参数错误", "目标分子数 / 搜索宽度 / 最大步数必须 ≥ 1")
            return None, None
        return tgt, p

    def start(self):
        if self.proc:
            return
        tgt, p = self.validate()
        if tgt is None:
            return
        pdb = glob.glob(os.path.join(PDB_DIR, tgt[1] + "*.pdb"))
        ts = time.strftime("%Y%m%d_%H%M%S")
        self.session_dir = os.path.join(OUTROOT, "easy_" + ts)
        os.makedirs(self.session_dir, exist_ok=True)
        cfgp = os.path.join(self.session_dir, "config.yml")
        gen = os.path.join(REPO, "src", "scripts", "gen_sample_config.py")
        tmpl = os.path.join(REPO, "configs", "sample_for_pdb_guided_l3.yml")
        if not os.path.exists(gen):
            messagebox.showerror("环境错误", "缺少配置生成器:\n" + gen); return
        gen_cmd = [PY, gen, "--template", tmpl, "--out", cfgp,
                   "--seed", str(int(p["seed"])), "--num-samples", str(int(p["num_samples"])),
                   "--beam", str(int(p["beam_size"])), "--max-steps", str(int(p["max_steps"])),
                   "--lam", str(p["lam"]), "--diversity-w", str(p["diversity_w"]),
                   "--guided", "1", "--relax-output", "1"]
        try:
            g = subprocess.run(gen_cmd, cwd=REPO, timeout=120, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, creationflags=FLAG)
            if g.returncode != 0 or not os.path.exists(cfgp):
                messagebox.showerror("配置生成失败",
                                     (g.stdout or b"").decode("utf-8", "ignore")[-400:] or "返回非零")
                return
        except Exception as e:
            messagebox.showerror("配置生成失败", "%s: %s" % (type(e).__name__, e)); return

        self.log_path = os.path.join(self.session_dir, "run.log")
        lf = open(self.log_path, "w", encoding="utf-8")
        cmd = [PY, os.path.join(REPO, "src", "sample_for_pdb.py"),
               "--pdb_path", pdb[0], "--center", tgt[2],
               "--config", cfgp.replace("\\", "/"), "--outdir", self.session_dir.replace("\\", "/")]
        try:
            self.proc = subprocess.Popen(cmd, cwd=REPO, stdout=lf, stderr=subprocess.STDOUT,
                                         creationflags=FLAG)
        except Exception as e:
            messagebox.showerror("启动失败", "%s: %s" % (type(e).__name__, e)); return
        self.log_pos = 0
        self.rows = []
        self.btn_start.set_state("disabled")
        self.btn_stop.set_state("normal")
        self.btn_open.set_state("disabled")
        self.tree.delete(*self.tree.get_children())
        self.mol_canvas.configure(image="", text="生成中…")
        self.progress["value"] = 0
        self.lbl_run.configure(text="正在生成分子…")
        self.logtxt.configure(state="normal"); self.logtxt.delete("1.0", "end")
        self.logtxt.configure(state="disabled")
        self._log_line("运行目录: %s" % self.session_dir)
        self._poll()

    def _log_line(self, s):
        self.logtxt.configure(state="normal")
        self.logtxt.insert("end", s.rstrip() + "\n")
        self.logtxt.see("end")
        self.logtxt.configure(state="disabled")

    def _poll(self):
        if not self.proc:
            return
        if self.log_path and os.path.exists(self.log_path):
            try:
                with open(self.log_path, encoding="utf-8", errors="ignore") as f:
                    f.seek(self.log_pos)
                    new = f.read()
                    self.log_pos = f.tell()
                for line in new.splitlines():
                    if "[Pool]" in line:
                        msg = line[line.index("[Pool]"):].strip()
                        self.lbl_run.configure(text="生成中 · " + msg)
                        self._status.set("%s · %s" % (self.sel_target, msg))
                        try:
                            fin = int(msg.split("Finished")[1].split("|")[0])
                            tgt = float(self.advvals["num_samples"].get())
                            self.progress["value"] = min(95, 100.0 * fin / max(1, tgt))
                        except Exception:
                            pass
                    elif "Success:" in line:
                        self._log_line(line[line.index("Success:"):].strip())
                    elif "Error" in line or "Traceback" in line or "WARNING" in line:
                        self._log_line(line.strip())
            except Exception:
                pass
        rc = self.proc.poll()
        if rc is None:
            self.root.after(1000, self._poll)
        else:
            self.proc = None
            self.btn_stop.set_state("disabled")
            if rc == 0:
                self.lbl_run.configure(text="采样完成，正在过滤与计算类药指标…")
                self._status.set("采样完成，正在过滤入库…")
                threading.Thread(target=self._postprocess, daemon=True).start()
            else:
                self.lbl_run.configure(text="采样进程异常退出 (代码 %d)" % rc)
                self._status.set("采样失败，详见日志")
                self.btn_start.set_state("normal")
                self.btn_open.set_state("normal")

    def _postprocess(self):
        lib = os.path.join(self.session_dir, "library")
        try:
            r = subprocess.run([PY, os.path.join(REPO, "src", "scripts", "build_library.py"),
                                "--runs", self.session_dir, "--library", lib],
                               cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               timeout=1800, creationflags=FLAG)
            csvp = os.path.join(lib, "compounds.csv")
            if r.returncode != 0 or not os.path.exists(csvp):
                self.root.after(0, lambda: (self.lbl_run.configure(text="入库异常，显示原始分子"),
                                            self._load_raw(),
                                            self.btn_start.set_state("normal"),
                                            self.btn_open.set_state("normal")))
                return
            rows = []
            with open(csvp, encoding="utf-8-sig") as f:
                for row in csv.DictReader(f):
                    try:
                        row["_qed"] = float(row["qed"])
                    except (TypeError, ValueError):
                        row["_qed"] = 0.0
                    rows.append(row)
            rows.sort(key=lambda x: -x["_qed"])
            self.root.after(0, lambda: self._fill(rows, lib))
        except Exception as e:
            msg = "%s: %s" % (type(e).__name__, e)
            self.root.after(0, lambda: (self.lbl_run.configure(text="过滤阶段异常"),
                                        self._load_raw(),
                                        messagebox.showwarning("过滤异常", msg),
                                        self.btn_start.set_state("normal"),
                                        self.btn_open.set_state("normal")))

    def _load_raw(self):
        f = glob.glob(os.path.join(self.session_dir, "**", "SMILES.txt"), recursive=True)
        if not f:
            return
        with open(f[0], encoding="utf-8-sig") as fh:
            lines = [x for x in fh.read().splitlines() if x.strip()]
        self.rows = [{"smiles": s, "qed": "", "mw": "", "logp": "", "sa": "",
                      "murcko_scaffold": "", "_qed": 0.0} for s in lines]
        self._fill_table()
        self.lbl_run.configure(text="已显示 %d 个原始分子（未过滤）" % len(lines))

    def _fill(self, rows, lib):
        self.lib_dir = lib
        self.rows = rows
        self._fill_table()
        q = [r["_qed"] for r in rows if r["_qed"] > 0]
        med = sorted(q)[len(q) // 2] if q else 0.0
        self.lbl_stat.configure(text="共 %d 个候选 · QED 中位 %.3f" % (len(rows), med))
        self.lbl_run.configure(text="完成！%d 个分子通过类药性过滤" % len(rows))
        self._status.set("生成完成：%d 个候选分子" % len(rows))
        self.progress["value"] = 100
        self.btn_start.set_state("normal")
        self.btn_open.set_state("normal")

    def _fill_table(self):
        if not hasattr(self, "tree"):
            return
        self.tree.delete(*self.tree.get_children())
        try:
            qmin = self.qed_min.get()
        except Exception:
            qmin = 0.0
        kw = self.q_text.get().strip().lower()
        data = []
        for r in self.rows:
            if r.get("qed", "") == "":          # 原始未过滤分子: 不参与 QED 阈值筛选
                data.append(r)
            elif r.get("_qed", 0.0) >= qmin:
                data.append(r)
        if kw:
            data = [r for r in data if kw in ",".join(str(v) for v in r.values()).lower()]
        if self.sort_col:
            key = self.sort_col

            def kf(r):
                v = r.get(key, "")
                try:
                    return (0, float(v))
                except (TypeError, ValueError):
                    return (1, str(v))
            data.sort(key=kf, reverse=self.sort_desc)
        else:
            data.sort(key=lambda r: -r.get("_qed", 0.0))
        for i, r in enumerate(data, 1):
            def num(v, f="%.2f"):
                try:
                    return f % float(v)
                except (TypeError, ValueError):
                    return "-"
            scaf = str(r.get("murcko_scaffold", "") or "")
            self.tree.insert("", "end", iid=str(i), tags=("odd" if i % 2 else "even",),
                             values=(i, num(r.get("qed"), "%.3f"), num(r.get("mw"), "%.0f"),
                                     num(r.get("logp")), num(r.get("sa")),
                                     scaf if scaf else "-", r.get("smiles", "")))
        if self.rows and len(data) != len(self.rows):
            self.lbl_stat.configure(text="显示 %d / %d 个候选" % (len(data), len(self.rows)))

    def _sort(self, col):
        if self.sort_col == col:
            self.sort_desc = not self.sort_desc
        else:
            self.sort_col, self.sort_desc = col, True
        self._fill_table()

    def _selected_smiles(self, tree=None):
        tree = tree or self.tree
        sel = tree.selection()
        if not sel:
            return None
        return tree.set(sel[0], "smiles") if tree is self.tree else None

    def _show_mol(self, event=None):
        sel = self.tree.selection()
        if not sel:
            return
        vals = self.tree.item(sel[0], "values")
        smi = vals[6]
        self.mol_info.configure(text="QED %s · MW %s · LogP %s · SA %s" % (vals[1], vals[2],
                                                                          vals[3], vals[4]))
        self.mol_canvas.configure(text="渲染中…", image="")
        self.mol_canvas._img = None
        # 按预览区实际宽度渲染, 否则 440x360 的图会被 Label 裁掉大半
        try:
            self.root.update_idletasks()
            avail = self.mol_canvas.master.winfo_width() - 10
            w = max(260, min(380, avail))
            h = int(w * 0.84)
        except Exception:
            w, h = 360, 300
        threading.Thread(target=self._render, args=(smi, self.mol_canvas, (w, h)),
                         daemon=True).start()

    def _render(self, smi, canvas, size=(440, 360)):
        png = os.path.join(os.environ.get("TEMP", "."), "_p2m_gui_%d_%d.png" % (os.getpid(),
                                                                               abs(hash(smi)) % 99999))
        try:
            code = ("import sys;from rdkit import Chem,RDLogger;RDLogger.DisableLog('rdApp.*');"
                    "from rdkit.Chem.Draw import rdMolDraw2D;"
                    "m=Chem.MolFromSmiles(sys.argv[1]);"
                    "d=rdMolDraw2D.MolDraw2DCairo(%d,%d);"
                    "o=d.drawOptions();o.bondLineWidth=2;"
                    "d.DrawMolecule(m);d.FinishDrawing();"
                    "open(sys.argv[2],'wb').write(d.GetDrawingText())" % size)
            r = subprocess.run([PY, "-c", code, smi, png], cwd=REPO, timeout=120,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               creationflags=FLAG)
            if os.path.exists(png):
                self.root.after(0, lambda: self._show_png(png, canvas))
            else:                       # 失败原因直接显示, 便于排查
                msg = (r.stdout or b"").decode("utf-8", "ignore").strip()[-160:] or "渲染进程无输出"
                self.root.after(0, lambda m=msg: canvas.configure(text="2D 渲染失败\n" + m, image=""))
        except Exception:
            pass

    def _show_png(self, path, canvas):
        try:
            img = tk.PhotoImage(file=path)
            # Label 若保留字符单位 width/height 会把图片裁掉, 须按图片实际像素重设
            canvas.configure(image=img, text="", width=img.width(), height=img.height())
            canvas._img = img
        except Exception:
            canvas.configure(text="2D 渲染失败", image="")

    def export_sdf(self):
        smi = self.tree.set(self.tree.selection()[0], "smiles") if self.tree.selection() else None
        if not smi:
            messagebox.showinfo("提示", "请先在列表中选择一个分子"); return
        dst = filedialog.asksaveasfilename(defaultextension=".sdf", initialfile="candidate.sdf",
                                           filetypes=[("3D 结构", "*.sdf")], title="导出 3D 结构")
        if not dst:
            return
        code = ("import sys;from rdkit import Chem,RDLogger;RDLogger.DisableLog('rdApp.*');"
                "from rdkit.Chem import AllChem;m=Chem.MolFromSmiles(sys.argv[1]);"
                "m=Chem.AddHs(m);AllChem.EmbedMolecule(m,randomSeed=42);"
                "AllChem.MMFFOptimizeMolecule(m);m=Chem.RemoveHs(m);"
                "open(sys.argv[2],'w',encoding='utf-8').write(Chem.MolToMolBlock(m))")
        try:
            subprocess.run([PY, "-c", code, smi, dst], cwd=REPO, timeout=120,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
                           creationflags=FLAG)
            messagebox.showinfo("导出完成", "已导出 3D 结构:\n" + dst)
        except Exception as e:
            messagebox.showerror("导出失败", "%s: %s" % (type(e).__name__, e))

    def export_csv(self):
        if not self.rows:
            messagebox.showinfo("提示", "当前没有可导出的结果"); return
        dst = filedialog.asksaveasfilename(defaultextension=".csv", initialfile="candidates.csv",
                                           filetypes=[("表格", "*.csv")], title="导出候选分子表")
        if not dst:
            return
        cols = ["smiles", "qed", "mw", "logp", "sa", "murcko_scaffold"]
        try:
            with open(dst, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(cols)
                for iid in self.tree.get_children():
                    v = self.tree.item(iid, "values")
                    w.writerow([v[6], v[1], v[2], v[3], v[4], v[5]])
            messagebox.showinfo("导出完成", "已导出 %d 行:\n%s" % (len(self.tree.get_children()), dst))
        except Exception as e:
            messagebox.showerror("导出失败", "%s: %s" % (type(e).__name__, e))

    def copy_smiles(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("提示", "请先选择一个分子"); return
        smi = self.tree.set(sel[0], "smiles")
        self.root.clipboard_clear(); self.root.clipboard_append(smi)
        self._status.set("已复制 SMILES 到剪贴板")

    def open_dir(self):
        if self.session_dir and os.path.isdir(self.session_dir):
            os.startfile(self.session_dir)

    def stop(self):
        if self.proc:
            self.proc.terminate()
            self.proc = None
            self.lbl_run.configure(text="已停止")
            self._status.set("已手动停止")
            self.btn_start.set_state("normal")
            self.btn_stop.set_state("disabled")
            self.btn_open.set_state("normal")

    # ================= Tab2 逻辑 =================
    def refresh_sources(self):
        srcs = []
        for base in (LIB_DIR, OUTROOT):
            if os.path.isdir(base):
                for p in sorted(glob.glob(os.path.join(base, "*.csv"))):
                    srcs.append(p)
        if os.path.isdir(LIB_DIR):
            for p in sorted(glob.glob(os.path.join(LIB_DIR, "**", "*.csv"), recursive=True)):
                if p not in srcs:
                    srcs.append(p)
        self._sources = srcs
        names = [os.path.relpath(p, LIB_DIR) if p.startswith(LIB_DIR) else p for p in srcs]
        self.cmb["values"] = names
        if names and not self.lib_src.get():
            self.cmb.current(0)
            self.load_source()
        self.lib_info.configure(text="共 %d 个数据文件" % len(srcs))

    def load_source(self):
        i = self.cmb.current()
        if i < 0 or i >= len(self._sources):
            return
        path = self._sources[i]
        try:
            with open(path, encoding="utf-8-sig", errors="ignore") as f:
                rd = csv.reader(f)
                header = next(rd)
                data = [r for r in rd]
        except Exception as e:
            messagebox.showerror("读取失败", "%s: %s" % (type(e).__name__, e)); return
        self._lib_rows, self._lib_cols = data, header
        self.lib_tree.delete(*self.lib_tree.get_children())
        self.lib_tree["columns"] = header
        for c in header:
            self.lib_tree.heading(c, text=c)
            self.lib_tree.column(c, width=max(110, measure(self.lib_tree, c, FONT, 9, True) + 36),
                                 anchor="w")
        for i2, r in enumerate(data[:3000], 1):
            self.lib_tree.insert("", "end", iid=str(i2), tags=("odd" if i2 % 2 else "even",),
                                 values=r)
        self.lib_info.configure(text="%s · %d 行 × %d 列" % (os.path.basename(path), len(data),
                                                            len(header)))
        self._lib_path = path
        if "smiles" in [c.lower() for c in header]:
            self._lib_smi_col = [c for c in header if c.lower() == "smiles"][0]
        else:
            self._lib_smi_col = next((c for c in header if "smiles" in c.lower()), None)

    def _lib_show_mol(self, _=None):
        if not getattr(self, "_lib_smi_col", None):
            return
        sel = self.lib_tree.selection()
        if not sel:
            return
        vals = self.lib_tree.item(sel[0], "values")
        idx = list(self.lib_tree["columns"]).index(self._lib_smi_col)
        smi = vals[idx] if idx < len(vals) else ""
        info = " · ".join("%s=%s" % (c, v) for c, v in zip(self.lib_tree["columns"], vals)
                          if c != self._lib_smi_col and v)[:400]
        self.lib_detail.configure(text=info)
        self.lib_canvas.configure(text="渲染中…", image="")
        self.lib_canvas._img = None
        if smi:
            threading.Thread(target=self._render, args=(smi, self.lib_canvas), daemon=True).start()

    def lib_export(self):
        if not self._lib_rows:
            messagebox.showinfo("提示", "没有可导出的数据"); return
        dst = filedialog.asksaveasfilename(defaultextension=".csv", initialfile="export.csv",
                                           filetypes=[("表格", "*.csv")])
        if not dst:
            return
        with open(dst, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f)
            w.writerow(self._lib_cols)
            w.writerows(self._lib_rows)
        messagebox.showinfo("导出完成", "已导出 %d 行:\n%s" % (len(self._lib_rows), dst))

    # ================= Tab3 逻辑 =================
    def refresh_reports(self):
        self.rep_tree.delete(*self.rep_tree.get_children())
        self._reps = []
        items = []
        for d in [LIB_DIR] + REPORT_DIRS:
            if not os.path.isdir(d):
                continue
            for ext in ("*.csv", "*.md", "*.txt", "*.pptx", "*.pdf"):
                for p in sorted(glob.glob(os.path.join(d, ext))):
                    items.append(p)
        for d in glob.glob(os.path.join(OUTROOT, "docking_*")):
            for p in sorted(glob.glob(os.path.join(d, "summary.csv"))):
                items.append(p)
        for i, p in enumerate(items, 1):
            try:
                st = os.stat(p)
                sz = "%.1f KB" % (st.st_size / 1024.0) if st.st_size < 1048576 else \
                     "%.1f MB" % (st.st_size / 1048576.0)
                tm = time.strftime("%Y-%m-%d %H:%M", time.localtime(st.st_mtime))
            except Exception:
                sz, tm = "?", "?"
            self._reps.append(p)
            parent = os.path.basename(os.path.dirname(p))
            label = ("%s/%s" % (parent, os.path.basename(p))) if parent else os.path.basename(p)
            self.rep_tree.insert("", "end", iid=str(i), tags=("odd" if i % 2 else "even",),
                                 values=(label, sz, tm))
        self._tag_report_rows()

    def _rep_selected(self):
        sel = self.rep_tree.selection()
        if not sel:
            return None
        return self._reps[int(sel[0]) - 1]

    def preview_report(self):
        p = self._rep_selected()
        if not p:
            return
        self.rep_path.configure(text=p)
        try:
            if p.lower().endswith((".csv", ".txt", ".md")):
                with open(p, encoding="utf-8-sig", errors="ignore") as f:
                    head = "".join(f.readlines()[:60])
                if p.lower().endswith(".csv"):
                    try:
                        import io
                        rd = list(csv.reader(io.StringIO(head)))
                        w = [max(len(str(r[i])) for r in rd if i < len(r)) for i in range(len(rd[0]))]
                        w = [min(x, 28) for x in w]
                        head = "\n".join("  ".join(str(r[i])[:w[i]].ljust(w[i])
                                                   for i in range(len(r))) for r in rd)
                    except Exception:
                        pass
            else:
                head = "(二进制或非文本文件, 请用「打开文件」查看)"
        except Exception as e:
            head = "读取失败: %s" % e
        self.rep_preview.configure(state="normal")
        self.rep_preview.delete("1.0", "end")
        self.rep_preview.insert("1.0", head)
        self.rep_preview.configure(state="disabled")

    def open_report(self):
        p = self._rep_selected()
        if p and os.path.exists(p):
            try:
                os.startfile(p)
            except Exception as e:
                messagebox.showerror("打开失败", str(e))

    def open_report_folder(self):
        p = self._rep_selected()
        if p:
            try:
                subprocess.Popen(["explorer", "/select,", os.path.normpath(p)], creationflags=FLAG)
            except Exception as e:
                messagebox.showerror("打开失败", str(e))

    # ================= Tab4 逻辑 =================
    def save_settings(self):
        global REPO, PY, PDB_DIR, OUTROOT, LIB_DIR
        vals = {k: e.get().strip() for k, e in self.set_entries.items()}
        try:
            with open(CFG, "w", encoding="utf-8") as f:
                json.dump(vals, f, ensure_ascii=False, indent=2)
        except Exception as e:
            messagebox.showerror("保存失败", str(e)); return
        REPO, PY = vals["repo"], vals["python"]
        PDB_DIR, OUTROOT, LIB_DIR = vals["pdb_dir"], vals["outroot"], vals["library"]
        self.refresh_sources()
        self.run_env_check()
        messagebox.showinfo("已保存", "设置已写入:\n%s\n（立即生效，无需重启）" % CFG)

    def open_cfg(self):
        if not os.path.exists(CFG):
            self.save_settings()
        try:
            os.startfile(CFG)
        except Exception as e:
            messagebox.showerror("打开失败", str(e))

    def run_env_check(self):
        """环境自检: 全程在后台线程执行 (torch 子进程要 5-15 秒),
        窗口与界面先显示出来, 结果稍后回填 —— 否则启动会被阻塞。"""
        self.chk_tree.delete(*self.chk_tree.get_children())
        self.chk_tree.insert("", "end", values=("检测中…",))
        threading.Thread(target=self._env_worker, daemon=True).start()

    def _env_worker(self):
        items = []
        status = None
        try:
            items.append(("Python 解释器", os.path.exists(PY), PY))
            ck = os.path.join(REPO, "models", "pretrained_Pocket2Mol.pt")
            items.append(("预训练权重", os.path.exists(ck), ck))
            items.append(("采样脚本 sample_for_pdb.py",
                          os.path.exists(os.path.join(REPO, "src", "sample_for_pdb.py")), REPO))
            items.append(("配置生成器 gen_sample_config.py",
                          os.path.exists(os.path.join(REPO, "src", "scripts", "gen_sample_config.py")), ""))
            items.append(("过滤入库 build_library.py",
                          os.path.exists(os.path.join(REPO, "src", "scripts", "build_library.py")), ""))
            for name, pdb, c, d in TARGETS:
                hit = glob.glob(os.path.join(PDB_DIR, pdb + "*.pdb"))
                items.append(("靶点 %s (%s)" % (name, pdb), bool(hit), hit[0] if hit else "未找到"))
            items.append(("输出目录可写", os.access(OUTROOT, os.W_OK) if os.path.isdir(OUTROOT)
                          else os.path.isdir(os.path.dirname(OUTROOT)), OUTROOT))
            items.append(("化合物库目录", os.path.isdir(LIB_DIR), LIB_DIR))
            try:
                r = subprocess.run([PY, "-c", "import torch;print('torch',torch.__version__,"
                                    "'cuda',torch.cuda.is_available())"], timeout=180,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   creationflags=FLAG)
                out = (r.stdout or b"").decode("utf-8", "ignore").strip().splitlines()
                last = out[-1] if out else ""
                items.append(("PyTorch / CUDA", r.returncode == 0, last))
                if r.returncode == 0 and "cuda True" in last:
                    status = ("● 环境正常 (GPU 可用)", self.T["ok"])
                elif r.returncode == 0:
                    status = ("● 环境正常 (CPU 模式)", self.T["warn"])
                else:
                    status = ("● 环境异常", self.T["err"])
            except Exception as e:
                items.append(("PyTorch / CUDA", False, str(e)))
                status = ("● 环境异常", self.T["err"])
        except Exception as e:
            items.append(("自检异常", False, str(e)))
        self.root.after(0, lambda: self._render_env(items, status))

    def _render_env(self, items, status):
        try:
            self.chk_tree.delete(*self.chk_tree.get_children())
            for i, (name, ok, detail) in enumerate(items):
                self.chk_tree.insert("", "end", tags=("odd" if i % 2 else "even",),
                                     values=("%s  %s%s" % ("✔" if ok else "✘", name,
                                                           ("   —   " + detail) if detail else ""),))
            self.chk_tree.tag_configure("odd", background=self.T["tree_alt"])
            self.chk_tree.tag_configure("even", background=self.T["card"])
            if status:
                self.lbl_env.configure(text=status[0], fg=status[1])
        except Exception:
            pass

    # ================= 关闭 =================
    def on_close(self):
        if self.proc is not None:
            if not messagebox.askyesno("确认退出", "正在生成分子，退出将终止本次采样。\n"
                                       "已产出的中间快照会保留。是否退出？"):
                return
            try:
                self.proc.terminate()
                try:
                    self.proc.wait(timeout=10)
                except Exception:
                    self.proc.kill()
            except Exception:
                pass
            self.proc = None
        try:
            self.root.destroy()
        except Exception:
            pass


HELP_TEXT = """7-eonmol · GPCR 靶向分子生成器 · 使用说明

【它能做什么】
输入一个 GPCR 靶点的口袋结构，用 Pocket2Mol 模型自回归地逐个原子生成小分子，
并通过引导束搜索（QED 类药性 + 合成可及性 + 多样性惩罚）优先生成类药分子，
再自动完成药物相似性过滤、构象精修与三维结构输出。

【操作】
1. 「生成分子」页点击靶点卡片选择靶点（A2A / D3 / 5-HT2B；B2AR 已于 2026-09-30 按项目决策停用）
2. 点击绿色「开始生成分子」，等待 8–15 分钟（GPU）
   · 进度条与状态栏实时显示束搜索池状态 Finished/Queue
   · 中途可「停止」，已产出的分子会保留
3. 完成后在下方表格浏览候选分子：点击任意行查看 2D 结构式
   · 类药分 QED 越高越好（>0.7 为优良）
   · 可用「类药分 ≥」滑块与搜索框筛选；点击表头可按该列排序
   · 选中后可「导出 3D 结构 (SDF)」「复制编码」「导出 CSV」

【其它功能】
· 候选库浏览：打开已有化合物库/分析结果 CSV，排序、筛选、看结构、导出
· 分析结果：集中列出化合物库 CSV、对接结果、桌面上的分析报告，双击即可打开
· 环境设置：修改仓库/解释器/靶点目录等路径并保存；一键运行环境自检
· 右上角「?」随时打开本说明，「☾/☀」切换浅色 / 深色主题

【高级参数】（一般无需修改）
· 目标分子数：本次希望生成并通过过滤的分子数量（默认 50）
· 搜索宽度 beam：束搜索保留的候选数，越大越慢但探索更充分（默认 100）
· 最大步数：单分子最多生长步数，也是运行时长的主要上限（默认 50）
· 引导强度 λ：化学分数对采样概率的影响权重（默认 3.0，越大越偏向高 QED）
· 多样性权重：对与已生成分子相似的结构施加的惩罚（默认 0.5）
· 随机种子：相同种子 + 相同参数可复现同一批分子

【常见问题】
· 提示找不到靶点文件：在「环境设置」页检查「靶点 PDB 目录」是否正确
· 提示环境错误：确认 Python 解释器路径存在且装有 torch/rdkit
· 结果为空：可能采样未产生完整分子，勾选「显示运行日志」查看后端输出
· GPU 显存不足：关闭其他占用显卡的程序；模型较小（8GB 显存足够）

【输出位置】
结果保存在「输出目录/easy_<时间戳>/」下：
  config.yml（本次参数）· SMILES.txt（分子编码）· SDF/（三维结构）
  library/compounds.csv（过滤后的候选库与各项类药指标）
"""


def main():
    import argparse
    ap = argparse.ArgumentParser(description=APP_TITLE)
    ap.add_argument("--tab", type=int, default=0, help="启动时打开的标签页序号 0-3")
    ap.add_argument("--theme", default="light", choices=["light", "dark"], help="初始主题")
    args = ap.parse_args()
    scale = enable_dpi()
    root = tk.Tk()
    root._p2m_scale = scale if scale and scale > 0 else 1.0
    if scale and scale > 0:
        try:
            root.tk.call("tk", "scaling", scale * 96.0 / 72.0)
        except Exception:
            pass
    app = App(root, theme=args.theme)
    try:
        if 0 <= args.tab < 4:
            app.nb.select(args.tab)
    except Exception:
        pass
    root.protocol("WM_DELETE_WINDOW", app.on_close)
    root.mainloop()


if __name__ == "__main__":
    main()
