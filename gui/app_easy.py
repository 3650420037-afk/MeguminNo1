# -*- coding: utf-8 -*-
"""GPCR 小分子药物生成器 — 外行友好版 GUI (打包用)

依赖: 仅 Python 标准库(Tkinter)。计算全部通过子进程调用 Pocket2Mol 环境。
默认后端: D:\\Miniconda3\\envs\\Pocket2Mol\\python.exe (可在本文件常量或 config.json 修改)。
"""
import os, sys, json, csv, re, threading, subprocess, time, glob
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

APP_TITLE = "GPCR 小分子药物生成器"
_HERE = os.path.dirname(os.path.abspath(__file__))
_EXE_DIR = os.path.dirname(sys.executable) if getattr(sys, "frozen", False) else _HERE
REPO = r"D:\MMModel\Pocket2Mol"
PY = r"D:\Miniconda3\envs\Pocket2Mol\python.exe"
PDB_DIR = r"D:\MMModel\靶点结构"
OUTROOT = os.path.join(REPO, "outputs")
CFG = os.path.join(_HERE, "gui_config.json")

# 部署模式: exe 同目录存在 repo\ 子文件夹时, 自动改用包内仓库与环境 (分发免配置)
if os.path.isdir(os.path.join(_EXE_DIR, "repo")):
    REPO = os.path.join(_EXE_DIR, "repo")
    PDB_DIR = os.path.join(REPO, "structures")
    OUTROOT = os.path.join(REPO, "outputs")
    for _cand in (os.path.join(_EXE_DIR, "env", "python.exe"),
                  os.path.join(_EXE_DIR, "env", "Scripts", "python.exe")):
        if os.path.exists(_cand):
            PY = _cand
            break

# 允许 exe 旁的 config.json 覆盖路径 (分发到其他机器时改这里)
try:
    _c = json.load(open(os.path.join(_EXE_DIR, "gui_config.json"), encoding="utf-8"))
    REPO = _c.get("repo", REPO); PY = _c.get("python", PY); OUTROOT = _c.get("outroot", OUTROOT)
    PDB_DIR = _c.get("pdb_dir", PDB_DIR)
except Exception:
    pass

TARGETS = [
    ("A2A 腺苷受体", "4EIY", " -0.4,8.5,17.1", "帕金森病 / 肿瘤免疫 / 炎症"),
    ("β2 肾上腺素受体", "2RH1", " -29.5,9.2,6.9", "哮喘 / 心血管疾病"),
    ("D3 多巴胺受体", "3PBL", " 0.085,-14.828,10.432", "精神分裂症 / 成瘾"),
    ("5-HT2B 血清素受体", "4IB4", " 22.448,18.284,11.726", "偏头痛 / 肺动脉高压"),
]

class App:
    def __init__(self, root):
        self.root = root
        root.title(APP_TITLE)
        root.geometry("980x720")
        self.proc = None
        self.log_path = None
        self.log_pos = 0
        self.session_dir = None
        self._build()

    # ---------- UI ----------
    def _build(self):
        pad = dict(padx=8, pady=6)
        frm = ttk.Frame(self.root); frm.pack(fill="both", expand=True, **pad)

        ttk.Label(frm, text="第一步：选择药物靶点", font=("Microsoft YaHei", 12, "bold")).pack(anchor="w")
        self.tgt = tk.StringVar(value=TARGETS[0][0])
        row = ttk.Frame(frm); row.pack(fill="x", **pad)
        for i, (name, pdb, center, disease) in enumerate(TARGETS):
            rb = ttk.Radiobutton(row, text=name, value=name, variable=self.tgt, command=self._on_target)
            rb.grid(row=0, column=i, sticky="w", padx=6)
        self.disease_var = tk.StringVar(value="相关疾病: " + TARGETS[0][3])
        ttk.Label(frm, textvariable=self.disease_var, foreground="#555").pack(anchor="w", padx=12)

        ttk.Separator(frm).pack(fill="x", pady=6)
        ttk.Label(frm, text="第二步：点击开始（默认使用已验证的最佳参数，约 10 分钟）",
                  font=("Microsoft YaHei", 12, "bold")).pack(anchor="w")
        btnrow = ttk.Frame(frm); btnrow.pack(fill="x", **pad)
        self.btn_start = tk.Button(btnrow, text="▶  开始生成药物分子", font=("Microsoft YaHei", 13, "bold"),
                                   bg="#2f8f4f", fg="white", height=2, command=self.start)
        self.btn_start.pack(side="left", padx=(0, 8), ipadx=20)
        self.btn_stop = tk.Button(btnrow, text="■ 停止", font=("Microsoft YaHei", 11),
                                  state="disabled", command=self.stop)
        self.btn_stop.pack(side="left", padx=4)
        self.btn_open = tk.Button(btnrow, text="📂 打开结果文件夹", state="disabled", command=self.open_dir)
        self.btn_open.pack(side="left", padx=4)

        self.progress = ttk.Progressbar(frm, mode="determinate", maximum=100)
        self.progress.pack(fill="x", **pad)
        self.status = tk.StringVar(value="就绪。选择靶点后点击开始。")
        ttk.Label(frm, textvariable=self.status, foreground="#0a4d8c").pack(anchor="w")

        ttk.Separator(frm).pack(fill="x", pady=6)
        ttk.Label(frm, text="第三步：查看生成的分子（点击行查看 2D 结构）",
                  font=("Microsoft YaHei", 12, "bold")).pack(anchor="w")
        mid = ttk.Frame(frm); mid.pack(fill="both", expand=True, **pad)
        cols = ("rank", "qed", "mw", "smiles")
        self.tree = ttk.Treeview(mid, columns=cols, show="headings", height=10)
        for c, w, txt in (("rank", 50, "序号"), ("qed", 70, "类药分"), ("mw", 70, "分子量"), ("smiles", 420, "分子结构编码(SMILES)")):
            self.tree.heading(c, text=txt); self.tree.column(c, width=w, anchor="w")
        self.tree.pack(side="left", fill="both", expand=True)
        self.tree.bind("<<TreeviewSelect>>", self._show_mol)
        right = ttk.Frame(mid); right.pack(side="left", fill="y", padx=6)
        self.mol_img = None
        self.mol_canvas = tk.Label(right, text="点击左侧分子\n此处显示 2D 结构", width=34, height=16, relief="groove")
        self.mol_canvas.pack()
        btns = ttk.Frame(right); btns.pack(fill="x", pady=4)
        ttk.Button(btns, text="导出该分子 3D 结构(SDF)", command=self.export_sdf).pack(fill="x", pady=2)
        ttk.Button(btns, text="复制分子编码", command=self.copy_smiles).pack(fill="x")

        adv = ttk.LabelFrame(frm, text="高级设置（一般无需修改）")
        adv.pack(fill="x", **pad)
        self.adv_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(adv, text="显示高级参数", variable=self.adv_on, command=self._toggle_adv).pack(anchor="w")
        self.advfrm = ttk.Frame(adv)
        self.advvals = {}
        for i, (k, v, txt) in enumerate([("num_samples", "50", "目标分子数"), ("beam_size", "100", "搜索宽度"),
                                         ("max_steps", "50", "最大步数"), ("lam", "3.0", "引导强度λ"),
                                         ("diversity_w", "0.5", "多样性权重"), ("seed", "9000", "随机种子")]):
            ttk.Label(self.advfrm, text=txt).grid(row=0, column=i*2, padx=4, sticky="e")
            e = ttk.Entry(self.advfrm, width=7); e.insert(0, v)
            e.grid(row=0, column=i*2+1, padx=4)
            self.advvals[k] = e
        self._toggle_adv()

        self.log = tk.Text(frm, height=6, state="disabled", font=("Consolas", 9))
        self.log.pack(fill="x", **pad)

    def _toggle_adv(self):
        if self.adv_on.get(): self.advfrm.pack(fill="x", padx=4, pady=2)
        else: self.advfrm.pack_forget()

    def _on_target(self):
        for name, pdb, center, disease in TARGETS:
            if name == self.tgt.get():
                self.disease_var.set("相关疾病: " + disease)

    def _log(self, s):
        self.log.configure(state="normal")
        self.log.insert("end", s.rstrip() + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    # ---------- 流程 ----------
    def start(self):
        if self.proc: return
        tgt = next(t for t in TARGETS if t[0] == self.tgt.get())
        pdb = glob.glob(os.path.join(PDB_DIR, tgt[1] + "*.pdb"))
        if not pdb:
            messagebox.showerror("错误", "找不到靶点蛋白文件: " + tgt[1]); return
        try:
            p = {k: float(e.get()) for k, e in self.advvals.items()}
        except ValueError:
            messagebox.showerror("错误", "高级参数必须是数字"); return
        ts = time.strftime("%Y%m%d_%H%M%S")
        self.session_dir = os.path.join(OUTROOT, "easy_" + ts)
        os.makedirs(self.session_dir, exist_ok=True)
        cfgp = os.path.join(self.session_dir, "config.yml")
        src = open(os.path.join(REPO, "configs", "sample_for_pdb_guided_l3.yml"), encoding="utf-8-sig").read()
        import re
        cfg = re.sub(r"seed: \d+", "seed: %d" % int(p["seed"]), src)
        cfg = re.sub(r"num_samples: \d+", "num_samples: %d" % int(p["num_samples"]), cfg)
        cfg = re.sub(r"beam_size: \d+", "beam_size: %d" % int(p["beam_size"]), cfg)
        cfg = re.sub(r"max_steps: \d+", "max_steps: %d" % int(p["max_steps"]), cfg)
        if "diversity_w" not in cfg:
            cfg = cfg.rstrip() + "\n    diversity_w: %s\n" % p["diversity_w"]
        else:
            cfg = re.sub(r"diversity_w: [\d.]+", "diversity_w: %s" % p["diversity_w"], cfg)
        open(cfgp, "w", encoding="utf-8").write(cfg)

        self.log_path = os.path.join(self.session_dir, "run.log")
        lf = open(self.log_path, "w", encoding="utf-8")
        cmd = [PY, "sample_for_pdb.py", "--pdb_path", pdb[0], "--center", tgt[2],
               "--config", cfgp.replace("\\", "/"), "--outdir", self.session_dir.replace("\\", "/")]
        try:
            self.proc = subprocess.Popen(cmd, cwd=REPO, stdout=lf, stderr=subprocess.STDOUT,
                                         creationflags=subprocess.CREATE_NO_WINDOW)
        except Exception as e:
            messagebox.showerror("启动失败", str(e)); return
        self.btn_start.configure(state="disabled"); self.btn_stop.configure(state="normal")
        self.btn_open.configure(state="disabled")
        self.tree.delete(*self.tree.get_children())
        self.mol_canvas.configure(text="点击左侧分子\n此处显示 2D 结构")
        self.progress["value"] = 0
        self._poll()

    def _poll(self):
        if not self.proc: return
        if self.log_path and os.path.exists(self.log_path):
            try:
                with open(self.log_path, encoding="utf-8", errors="ignore") as f:
                    f.seek(self.log_pos); new = f.read(); self.log_pos = f.tell()
                for line in new.splitlines():
                    if "[Pool]" in line:
                        self.status.set("正在生成分子… " + line[line.index("[Pool]"):].strip())
                    elif "Success:" in line:
                        self._log(line[line.index("Success:"):].strip())
                if "num_samples" in new and new.count("[Pool]"):
                    m = [l for l in new.splitlines() if "[Pool]" in l][-1]
                    try:
                        fin = int(m.split("Finished")[1].split("|")[0])
                        target = float(self.advvals["num_samples"].get())
                        self.progress["value"] = min(95, 100 * fin / max(1, target))
                    except Exception: pass
            except Exception:
                pass
        rc = self.proc.poll()
        if rc is None:
            self.root.after(1000, self._poll)
        else:
            self.proc = None
            self.btn_stop.configure(state="disabled")
            if rc == 0:
                self.status.set("生成完成，正在自动过滤与计算类药指标…")
                threading.Thread(target=self._postprocess, daemon=True).start()
            else:
                self.status.set("采样进程退出(代码 %d)。详见日志文件。" % rc)
                self.btn_start.configure(state="normal"); self.btn_open.configure(state="normal")

    def _postprocess(self):
        lib = os.path.join(self.session_dir, "library")
        r = subprocess.run([PY, os.path.join(REPO, "scripts", "build_library.py"),
                            "--runs", self.session_dir, "--library", lib],
                           cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600,
                           creationflags=subprocess.CREATE_NO_WINDOW)
        csvp = os.path.join(lib, "compounds.csv")
        if r.returncode != 0 or not os.path.exists(csvp):
            self.root.after(0, lambda: (self.status.set("过滤完成但入库异常，显示原始分子。"),
                                        self._load_raw(), self.btn_start.configure(state="normal"),
                                        self.btn_open.configure(state="normal")))
            return
        rows = list(csv.DictReader(open(csvp, encoding="utf-8")))
        rows.sort(key=lambda x: -float(x["qed"]))
        self.root.after(0, lambda: self._fill(rows, lib))

    def _load_raw(self):
        f = glob.glob(os.path.join(self.session_dir, "**", "SMILES.txt"), recursive=True)
        if not f: return
        for i, s in enumerate([x for x in open(f[0], encoding="utf-8").read().splitlines() if x.strip()], 1):
            self.tree.insert("", "end", iid=str(i), values=(i, "-", "-", s))

    def _fill(self, rows, lib):
        self.lib_dir = lib
        self.tree.delete(*self.tree.get_children())
        for i, r in enumerate(rows, 1):
            self.tree.insert("", "end", iid=str(i),
                             values=(i, "%.3f" % float(r["qed"]), "%.0f" % float(r["mw"]), r["smiles"]))
        self.status.set("完成！共 %d 个通过药物过滤的候选分子。点击行查看结构。" % len(rows))
        self.progress["value"] = 100
        self.btn_start.configure(state="normal"); self.btn_open.configure(state="normal")

    def _selected_smiles(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("提示", "请先在列表中点击一个分子"); return None
        return self.tree.set(sel[0], "smiles")

    def _show_mol(self):
        smi = self._selected_smiles()
        if not smi: return
        png = os.path.join(self.session_dir, "_mol.png")
        threading.Thread(target=self._render, args=(smi, png), daemon=True).start()

    def _render(self, smi, png):
        try:
            code = ("import sys; from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*');"
                    "from rdkit.Chem.Draw import rdMolDraw2D;"
                    "m=Chem.MolFromSmiles(sys.argv[1]);"
                    "d=rdMolDraw2D.MolDraw2DCairo(420,340);"
                    "rdMolDraw2D.PrepareAndDrawMolecule(d,m); d.FinishDrawing();"
                    "open(sys.argv[2],'wb').write(d.GetDrawingText())")
            subprocess.run([PY, "-c", code, smi, png], cwd=REPO, timeout=60,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           creationflags=subprocess.CREATE_NO_WINDOW)
            if os.path.exists(png):
                self.png_path = png
                self.root.after(0, self._show_png)
        except Exception:
            pass

    def _show_png(self):
        try:
            self.mol_img = tk.PhotoImage(file=self.png_path)
            self.mol_canvas.configure(image=self.mol_img, text="")
        except Exception:
            self.mol_canvas.configure(text="2D 渲染失败")

    def export_sdf(self):
        smi = self._selected_smiles()
        if not smi: return
        dst = filedialog.asksaveasfilename(defaultextension=".sdf", initialfile="candidate.sdf",
                                           filetypes=[("3D 结构", "*.sdf")])
        if not dst: return
        code = ("import sys; from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*');"
                "from rdkit.Chem import AllChem; m=Chem.MolFromSmiles(sys.argv[1]);"
                "m=Chem.AddHs(m); AllChem.EmbedMolecule(m, randomSeed=42);"
                "AllChem.MMFFOptimizeMolecule(m); m=Chem.RemoveHs(m);"
                # 用 Python open() + MolToMolBlock 写文件: RDKit 的 C++ 文件 API
                # 在 Windows 上打不开含中文的路径 (用户导出的目标目录可能含中文)
                "open(sys.argv[2], 'w', encoding='utf-8').write(Chem.MolToMolBlock(m))")
        try:
            subprocess.run([PY, "-c", code, smi, dst], cwd=REPO, timeout=60,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
                           creationflags=subprocess.CREATE_NO_WINDOW)
            messagebox.showinfo("完成", "已导出 3D 结构:\n" + dst)
        except Exception as e:
            messagebox.showerror("导出失败", str(e))

    def copy_smiles(self):
        smi = self._selected_smiles()
        if smi:
            self.root.clipboard_clear(); self.root.clipboard_append(smi)
            self.status.set("已复制到剪贴板")

    def open_dir(self):
        if self.session_dir: os.startfile(self.session_dir)

    def stop(self):
        if self.proc:
            self.proc.terminate(); self.proc = None
            self.status.set("已停止。")
            self.btn_start.configure(state="normal"); self.btn_stop.configure(state="disabled")
            self.btn_open.configure(state="normal")

if __name__ == "__main__":
    root = tk.Tk()
    try:
        from tkinter import font
        default = font.nametofont("TkDefaultFont"); default.configure(family="Microsoft YaHei", size=10)
    except Exception:
        pass
    app = App(root)
    root.mainloop()
