# -*- coding: utf-8 -*-
"""UI 组件库 (纯 Tkinter 实现现代观感)

- 双主题 (浅色 / 深色) 调色板与 ttk 样式配置
- 圆角按钮 (Canvas 自绘, 支持 hover/press/disabled/图标)
- 靶点选择卡片
- 分区标题、徽章、卡片容器
"""
import tkinter as tk
from tkinter import ttk
from tkinter import font as tkfont
import base64
import os, sys

FONT = "Microsoft YaHei UI"
FONT_MONO = "Consolas"

THEMES = {
    "light": dict(
        bg="#EEF2F7", card="#FFFFFF", card2="#F8FAFC", border="#DCE3EC", fg="#0F172A",
        sub="#5B6B82", accent="#2563EB", accent_dk="#1D4ED8", accent_lt="#E3EDFF",
        ok="#16A34A", ok_dk="#12813C", ok_lt="#E6F6EC", warn="#D97706", err="#DC2626",
        header="#1E293B", tree_alt="#F7FAFD", sel="#2563EB", entry="#FFFFFF",
        log_bg="#0E1726", log_fg="#CBD5E1", disabled="#A9B4C4", chip="#EAF0F8",
    ),
    "dark": dict(
        bg="#0B1220", card="#141E30", card2="#0F1829", border="#26344A", fg="#E8EEF8",
        sub="#93A4BF", accent="#3B82F6", accent_dk="#2563EB", accent_lt="#1B2C4A",
        ok="#22C55E", ok_dk="#16803C", ok_lt="#12301F", warn="#F59E0B", err="#EF4444",
        header="#0A1120", tree_alt="#111B2C", sel="#2563EB", entry="#0E1727",
        log_bg="#05090F", log_fg="#C7D3E4", disabled="#4A5A73", chip="#1B2739",
    ),
}


def round_rect(cv, x1, y1, x2, y2, r, **kw):
    """在 Canvas 上画圆角矩形 (smooth spline 近似), 返回 item id"""
    r = max(1, min(r, (x2 - x1) / 2, (y2 - y1) / 2))
    pts = [x1 + r, y1, x2 - r, y1, x2, y1, x2, y1 + r, x2, y2 - r, x2, y2,
           x2 - r, y2, x1 + r, y2, x1, y2, x1, y2 - r, x1, y1 + r, x1, y1]
    return cv.create_polygon(pts, smooth=True, **kw)


def measure(widget, text, family=FONT, size=10, bold=False):
    """按当前 DPI/缩放实测文本像素宽度 (高 DPI 下字体变宽, 固定宽度会被裁)"""
    try:
        f = tkfont.Font(root=widget, family=family, size=size,
                        weight="bold" if bold else "normal")
        return f.measure(text)
    except Exception:
        return int(len(text) * size * 0.75)


class RoundedButton(tk.Canvas):
    """圆角按钮. kind: primary / success / ghost / danger / warn"""

    def __init__(self, master, text="", command=None, kind="primary", width=170, height=40,
                 radius=12, font=None, theme=None, icon="", bg=None):
        self.T = theme or THEMES["light"]
        self._parent_bg = bg or master.cget("bg") if isinstance(master, (tk.Frame, tk.Canvas, tk.Tk)) else self.T["bg"]
        self.font = font or (FONT, 11, "bold")
        self.text = text
        self.icon = icon
        label = (icon + "  " + text).strip() if icon else text
        # 宽度自适应: 保证高 DPI 下文字不被裁
        need = measure(master, label, self.font[0], self.font[1],
                       len(self.font) > 2 and "bold" in str(self.font[2])) + 40
        height = max(height, self.font[1] * 3 + 4)
        width = max(width, need)
        super().__init__(master, width=width, height=height, highlightthickness=0,
                         bd=0, bg=self._parent_bg)
        self.command = command
        self.kind = kind
        self.radius = radius
        self._enabled = True
        self._hover = False
        self._press = False
        self._cw, self._ch = width, height
        self.bind("<Configure>", self._on_conf)
        self.bind("<Enter>", self._on_enter)
        self.bind("<Leave>", self._on_leave)
        self.bind("<ButtonPress-1>", self._on_press)
        self.bind("<ButtonRelease-1>", self._on_release)
        self.configure(cursor="hand2")
        self._redraw()

    # -- 颜色 --
    def _palette(self):
        T = self.T
        if not self._enabled:
            return T["disabled"], T["card"] if T["card"] == "#FFFFFF" else "#FFFFFF"
        if self.kind == "primary":
            bg, fg = T["accent"], "#FFFFFF"
        elif self.kind == "success":
            bg, fg = T["ok"], "#FFFFFF"
        elif self.kind == "danger":
            bg, fg = T["err"], "#FFFFFF"
        elif self.kind == "warn":
            bg, fg = T["warn"], "#FFFFFF"
        else:  # ghost: 白底描边, hover 时浅蓝
            bg, fg = (T["accent_lt"] if self._hover else T["card"]), T["fg"]
        if self._press:
            bg = self._shade(bg, 0.86)
        elif self._hover:
            bg = self._shade(bg, 1.10)
        return bg, fg

    @staticmethod
    def _shade(hexcol, k):
        hexcol = hexcol.lstrip("#")
        r, g, b = (int(hexcol[i:i + 2], 16) for i in (0, 2, 4))
        f = lambda v: max(0, min(255, int(v * k)))
        return "#%02X%02X%02X" % (f(r), f(g), f(b))

    def _on_conf(self, ev):
        self._cw, self._ch = ev.width, ev.height
        self._redraw()

    def _redraw(self):
        T = self.T
        self.delete("all")
        bg, fg = self._palette()
        outline = T["border"] if (self.kind == "ghost" and self._enabled) else bg
        round_rect(self, 1, 1, self._cw - 1, self._ch - 1, self.radius, fill=bg,
                   outline=outline, width=1)
        label = (self.icon + "  " + self.text).strip() if self.icon else self.text
        self.create_text(self._cw / 2, self._ch / 2 + 1, text=label, fill=fg, font=self.font)

    # -- 事件 --
    def _on_enter(self, _):
        if self._enabled:
            self._hover = True
            self._redraw()

    def _on_leave(self, _):
        self._hover = self._press = False
        self._redraw()

    def _on_press(self, _):
        if self._enabled:
            self._press = True
            self._redraw()

    def _on_release(self, ev):
        if not self._enabled:
            return
        was = self._press
        self._press = False
        self._redraw()
        if was and 0 <= ev.x <= self._cw and 0 <= ev.y <= self._ch and self.command:
            self.command()

    # -- 公开 API --
    def set_state(self, state):
        self._enabled = (state != "disabled")
        self.configure(cursor="hand2" if self._enabled else "arrow")
        self._redraw()

    def set_text(self, text):
        self.text = text
        self._redraw()

    def set_kind(self, kind):
        self.kind = kind
        self._redraw()

    def set_theme(self, T, parent_bg=None):
        self.T = T
        if parent_bg:
            self._parent_bg = parent_bg
            self.configure(bg=parent_bg)
        self._redraw()


class TargetCard(tk.Canvas):
    """靶点选择卡片 (标题 / PDB / 适应症), 选中时高亮"""

    def __init__(self, master, title, pdb, center, disease, command=None, theme=None,
                 width=250, height=110):
        self.T = theme or THEMES["light"]
        # 宽度自适应: 标题 + PDB 徽章 / 口袋中心 / 适应症三行都要放得下
        title_w = measure(master, title, FONT, 12, True)
        badge_w = measure(master, pdb, FONT_MONO, 9, True) + 22
        need = max(18 + title_w + 16 + badge_w + 18,
                   18 + measure(master, "口袋中心 " + center.strip(), FONT_MONO, 8) + 18,
                   18 + measure(master, disease, FONT, 9) + 100,
                   18 + measure(master, disease, FONT, 9) + 18)
        width = max(width, need)
        height = max(height, 112)
        super().__init__(master, width=width, height=height, highlightthickness=0, bd=0,
                         bg=self.T["bg"])
        self.title, self.pdb, self.center, self.disease = title, pdb, center, disease
        self.command = command
        self.selected = False
        self._hover = False
        self._cw, self._ch = width, height
        self.bind("<Configure>", self._on_conf)
        self.bind("<Enter>", lambda e: self._set_hover(True))
        self.bind("<Leave>", lambda e: self._set_hover(False))
        self.bind("<Button-1>", lambda e: self.command and self.command(self))
        self.configure(cursor="hand2")
        self._redraw()

    def _on_conf(self, ev):
        self._cw, self._ch = ev.width, ev.height
        self._redraw()

    def _set_hover(self, v):
        self._hover = v
        self._redraw()

    def set_selected(self, v):
        self.selected = v
        self._redraw()

    def set_theme(self, T, parent_bg):
        self.T = T
        self.configure(bg=parent_bg)
        self._redraw()

    def _redraw(self):
        T = self.T
        self.delete("all")
        w, h = self._cw, self._ch
        border = T["accent"] if self.selected else (T["accent_lt"] if self._hover else T["border"])
        fill = T["accent_lt"] if self.selected else (T["card2"] if self._hover else T["card"])
        round_rect(self, 1, 1, w - 1, h - 1, 14, fill=fill, outline=border, width=2 if self.selected else 1)
        if self.selected:
            round_rect(self, 0, 8, 5, h - 8, 2, fill=T["accent"], outline=T["accent"])
        x = 18
        # 第一行: 靶点名称 (左) + PDB 徽章 (右)
        self.create_text(x, 26, text=self.title, anchor="w", fill=T["fg"],
                         font=(FONT, 12, "bold"))
        bw = measure(self, self.pdb, FONT_MONO, 9, True) + 20
        round_rect(self, w - 18 - bw, 14, w - 18, 38, 8, fill=T["accent_lt"], outline=T["accent_lt"])
        self.create_text(w - 18 - bw / 2, 26, text=self.pdb, fill=T["accent"],
                         font=(FONT_MONO, 9, "bold"))
        # 第二行: 口袋中心
        self.create_text(x, 58, text="口袋中心 " + self.center.strip(), anchor="w",
                         fill=T["sub"], font=(FONT_MONO, 8))
        # 第三行: 适应症
        self.create_text(x, 86, text=self.disease, anchor="w", fill=T["sub"], font=(FONT, 9))
        if self.selected:
            self.create_text(w - 18, 86, text="✓ 已选", anchor="e", fill=T["ok"],
                             font=(FONT, 9, "bold"))


class SectionTitle(tk.Frame):
    """左侧色条 + 标题 + 可选副标题"""

    def __init__(self, master, text, theme=None, sub=None):
        T = theme or THEMES["light"]
        super().__init__(master, bg=T["bg"])
        self.T = T
        bar = tk.Frame(self, bg=T["accent"], width=4, height=17)
        bar.pack(side="left", padx=(0, 8))
        self.lbl = tk.Label(self, text=text, bg=T["bg"], fg=T["fg"], font=(FONT, 12, "bold"))
        self.lbl.pack(side="left")
        self.sub = None
        if sub:
            self.sub = tk.Label(self, text=sub, bg=T["bg"], fg=T["sub"], font=(FONT, 9))
            self.sub.pack(side="left", padx=10)

    def set_theme(self, T):
        self.T = T
        self.configure(bg=T["bg"])
        self.lbl.configure(bg=T["bg"], fg=T["fg"])
        for c in self.winfo_children():
            if isinstance(c, tk.Frame):
                c.configure(bg=T["accent"])
        if self.sub:
            self.sub.configure(bg=T["bg"], fg=T["sub"])


class Card(tk.Frame):
    """带边框的卡片容器 (内部再放内容)"""

    def __init__(self, master, theme=None, pad=14):
        T = theme or THEMES["light"]
        super().__init__(master, bg=T["card"], highlightthickness=1,
                         highlightbackground=T["border"], highlightcolor=T["border"], bd=0)
        self.T = T
        self.inner = tk.Frame(self, bg=T["card"])
        self.inner.pack(fill="both", expand=True, padx=pad, pady=pad)

    def set_theme(self, T):
        self.T = T
        self.configure(bg=T["card"], highlightbackground=T["border"], highlightcolor=T["border"])
        self.inner.configure(bg=T["card"])


def apply_theme(root, style, T):
    """把调色板套到 ttk 样式上 (可反复调用以实现主题切换)"""
    style.theme_use("clam")
    root.configure(bg=T["bg"])

    style.configure(".", background=T["bg"], foreground=T["fg"], font=(FONT, 10),
                    fieldbackground=T["entry"], bordercolor=T["border"], focuscolor=T["accent"])
    for name in ("TFrame", "TLabel", "TLabelframe", "TLabelframe.Label"):
        style.configure(name, background=T["bg"], foreground=T["fg"])
    style.configure("TLabel", font=(FONT, 10))
    style.configure("Sub.TLabel", foreground=T["sub"])
    style.configure("Card.TFrame", background=T["card"])
    style.configure("Card.TLabel", background=T["card"], foreground=T["fg"])
    style.configure("CardSub.TLabel", background=T["card"], foreground=T["sub"])

    style.configure("TButton", padding=(12, 7), relief="flat", background=T["chip"],
                    foreground=T["fg"], borderwidth=0, focusthickness=0)
    style.map("TButton", background=[("active", T["accent_lt"]), ("disabled", T["card2"])],
              foreground=[("disabled", T["disabled"])])

    style.configure("TCheckbutton", background=T["bg"], foreground=T["fg"], focuscolor=T["bg"])
    style.map("TCheckbutton", background=[("active", T["bg"])])
    style.configure("TRadiobutton", background=T["bg"], foreground=T["fg"], focuscolor=T["bg"])
    style.map("TRadiobutton", background=[("active", T["bg"])])

    style.configure("TEntry", fieldbackground=T["entry"], foreground=T["fg"],
                    insertcolor=T["fg"], bordercolor=T["border"], lightcolor=T["border"],
                    darkcolor=T["border"], padding=5)
    style.configure("TCombobox", fieldbackground=T["entry"], background=T["chip"],
                    foreground=T["fg"], arrowcolor=T["fg"], bordercolor=T["border"], padding=4)
    style.map("TCombobox", fieldbackground=[("readonly", T["entry"])],
              foreground=[("readonly", T["fg"])])
    root.option_add("*TCombobox*Listbox.background", T["card"])
    root.option_add("*TCombobox*Listbox.foreground", T["fg"])
    root.option_add("*TCombobox*Listbox.selectBackground", T["accent"])

    style.configure("Treeview", background=T["card"], fieldbackground=T["card"],
                    foreground=T["fg"], rowheight=30, borderwidth=0, font=(FONT, 9))
    style.map("Treeview", background=[("selected", T["sel"])],
              foreground=[("selected", "#FFFFFF")])
    style.configure("Treeview.Heading", background=T["header"], foreground="#FFFFFF",
                    relief="flat", font=(FONT, 9, "bold"), padding=(6, 7))
    style.map("Treeview.Heading", background=[("active", T["accent"])])

    style.configure("Horizontal.TProgressbar", troughcolor=T["chip"], background=T["accent"],
                    bordercolor=T["chip"], lightcolor=T["accent"], darkcolor=T["accent"],
                    thickness=16, borderwidth=0)
    style.configure("Vertical.TScrollbar", background=T["chip"], troughcolor=T["bg"],
                    bordercolor=T["bg"], arrowcolor=T["sub"])
    style.configure("Horizontal.TScrollbar", background=T["chip"], troughcolor=T["bg"],
                    bordercolor=T["bg"], arrowcolor=T["sub"])

    style.configure("TNotebook", background=T["bg"], borderwidth=0, tabmargins=(6, 6, 6, 0))
    style.configure("TNotebook.Tab", background=T["chip"], foreground=T["sub"],
                    padding=(20, 10), font=(FONT, 10, "bold"), borderwidth=0)
    style.map("TNotebook.Tab",
              background=[("selected", T["accent"]), ("active", T["accent_lt"])],
              foreground=[("selected", "#FFFFFF"), ("active", T["fg"])])

    style.configure("TSeparator", background=T["border"])
    return style


def app_icon(root):
    """设置窗口图标 (运行时从内嵌 base64 PNG 解码, 无需外部文件)"""
    try:
        from icon_b64 import ICON_B64
    except Exception:
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from icon_b64 import ICON_B64
        except Exception:
            return None
    try:
        img = tk.PhotoImage(data=base64.b64decode(ICON_B64))
        root.iconphoto(True, img)
        root._icon_ref = img      # 防 GC
        return img
    except Exception:
        return None
