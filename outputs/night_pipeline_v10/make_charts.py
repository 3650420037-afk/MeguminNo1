# -*- coding: utf-8 -*-
"""为 PPT 重新生成 4 张数据图（原图为旧权重数据的光栅图，无法编辑）。"""
import csv, os, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
RDLogger.DisableLog("rdApp.*")

ROOT = r"D:\MMModel"
OUT = os.path.join(ROOT, "outputs", "night_pipeline_v10", "charts")
os.makedirs(OUT, exist_ok=True)

CJK = "Microsoft YaHei"
plt.rcParams["font.sans-serif"] = [CJK, "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

BG = "#0C0F1A"
FG = "#E5E7EB"
MUTED = "#9AA6B2"
TEAL = "#1E7A70"
TEAL_L = "#6BC5BD"
GOLD = "#D9A441"


def style(ax):
    ax.set_facecolor(BG)
    for s in ax.spines.values():
        s.set_color("#2A3444")
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.yaxis.label.set_color(MUTED)
    ax.xaxis.label.set_color(MUTED)


# ---------------- 1) P15 生成批次统计 ----------------
tgt = ["A2A", "D3", "5-HT2B"]
out_n = [583, 254, 314]
in_n = [395, 106, 206]
batch = [11, 5, 6]
fig, ax = plt.subplots(figsize=(5.38, 3.20), dpi=200)
fig.patch.set_facecolor(BG)
style(ax)
x = np.arange(len(tgt))
w = 0.36
b1 = ax.bar(x - w / 2, out_n, w, color=TEAL_L, label="产出量")
b2 = ax.bar(x + w / 2, in_n, w, color=TEAL, label="入库数量")
for bars in (b1, b2):
    for r in bars:
        ax.annotate("%d" % r.get_height(), (r.get_x() + r.get_width() / 2, r.get_height()),
                    ha="center", va="bottom", color=FG, fontsize=8)
ax.set_ylabel("数量(个)", fontsize=9)
ax.set_ylim(0, 700)
ax.set_xticks(x); ax.set_xticklabels(tgt, color=FG, fontsize=9)
ax2 = ax.twinx()
ax2.plot(x, batch, color=GOLD, marker="o", ms=4, lw=1.4, label="生成批次")
for xi, v in zip(x, batch):
    ax2.annotate("%d" % v, (xi, v), textcoords="offset points", xytext=(0, 6), ha="center", color=GOLD, fontsize=8)
ax2.set_ylabel("批次", fontsize=9, color=MUTED)
ax2.set_ylim(0, 14)
ax2.tick_params(colors=MUTED, labelsize=8)
for s in ax2.spines.values():
    s.set_color("#2A3444")
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
leg = ax.legend(h1 + h2, l1 + l2, loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=3,
                frameon=False, fontsize=8, labelcolor=FG)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "p15_batch.png"), facecolor=BG, bbox_inches="tight")
plt.close(fig)

# ---------------- 2) P17 QED 分布 ----------------
labels = ["0.4-0.6", "0.6-0.8", "0.8-1.0"]
vals = [39, 241, 115]
fig, ax = plt.subplots(figsize=(5.50, 2.91), dpi=200)
fig.patch.set_facecolor(BG)
style(ax)
bars = ax.bar(labels, vals, 0.55, color=[TEAL, TEAL_L, TEAL_L])
for r in bars:
    ax.annotate("%d" % r.get_height(), (r.get_x() + r.get_width() / 2, r.get_height()),
                ha="center", va="bottom", color=FG, fontsize=9)
ax.set_ylim(0, 280)
ax.set_yticks([0, 70, 140, 210, 280])
ax.tick_params(axis="x", colors=FG, labelsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "p17_qed.png"), facecolor=BG, bbox_inches="tight")
plt.close(fig)

# ---------------- 读取新颖性分布数据（A2A 395 vs 737 已知） ----------------
known = [s.strip() for s in open(os.path.join(ROOT, "data", "known_drugs", "A2A_known_set.txt"),
                                encoding="utf-8") if s.strip()]
kfps, kmols = [], []
for s in known:
    m = Chem.MolFromSmiles(s)
    if not m:
        continue
    kmols.append(m)
    kfps.append(AllChem.GetMorganFingerprint(m, 2))
lib = []
for p in [os.path.join(ROOT, "results", "compounds_A2A.csv")]:
    for r in csv.DictReader(open(p, encoding="utf-8-sig", newline="")):
        lib.append(r["smiles"])
lfps, lmols = [], []
sims = []
for s in lib:
    m = Chem.MolFromSmiles(s)
    if not m:
        continue
    fp = AllChem.GetMorganFingerprint(m, 2)
    lmols.append(m); lfps.append(fp)
    sims.append(max(DataStructs.BulkTanimotoSimilarity(fp, kfps)))
sims = np.array(sorted(sims))
print("A2A 库 %d 分子, 相似度中位 %.3f 最大 %.3f, <0.4 %.1f%%" % (
    len(sims), np.median(sims), sims[-1], 100 * (sims < 0.4).mean()))

# ---------------- 3) P20 Tanimoto 分布 ----------------
fig, ax = plt.subplots(figsize=(5.79, 2.06), dpi=200)
fig.patch.set_facecolor(BG)
style(ax)
bins = np.linspace(0, 1, 41)
ax.hist(sims, bins=bins, color=TEAL_L, alpha=0.85, edgecolor="none")
ax.set_xlim(0, 1); ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticks([])
ax.axvline(float(np.median(sims)), color=GOLD, ls="--", lw=1.2)
ax.annotate("中位 %.3f" % float(np.median(sims)), (float(np.median(sims)), ax.get_ylim()[1] * 0.82),
            color=GOLD, fontsize=8, ha="left", xytext=(4, 0), textcoords="offset points")
for s in ["left"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "p20_tanimoto.png"), facecolor=BG, bbox_inches="tight")
plt.close(fig)

# ---------------- 4) P20 化学空间散点（PCA） ----------------
def fp_array(mols):
    arr = np.zeros((len(mols), 2048), dtype=np.float32)
    for i, m in enumerate(mols):
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr[i])
    return arr
Xk, Xl = fp_array(kmols), fp_array(lmols)
X = np.vstack([Xk, Xl])
X = X - X.mean(axis=0)
u, s, vt = np.linalg.svd(X, full_matrices=False)
P = X @ vt[:2].T
fig, ax = plt.subplots(figsize=(5.79, 2.06), dpi=200)
fig.patch.set_facecolor(BG)
style(ax)
ax.scatter(P[:len(Xk), 0], P[:len(Xk), 1], s=3, c="#5A6B7A", alpha=0.7, label="对照集(已知分子)")
ax.scatter(P[len(Xk):, 0], P[len(Xk):, 1], s=4, c=TEAL_L, alpha=0.85, label="生成分子(A2A 库)")
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.legend(loc="lower right", frameon=False, fontsize=7, labelcolor=MUTED, markerscale=2)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "p20_cluster.png"), facecolor=BG, bbox_inches="tight")
plt.close(fig)

# ---------------- 5) P16 漏斗 ----------------
stages = [("采样产出", 583), ("规范化去重", 576), ("PAINS 警示", 560), ("BRENK 警示", 478),
          ("MW [250,500]", 429), ("LogP [1,5]", 395), ("QED/SA/元素/环", 395)]
drops = [None, 7, 16, 82, 49, 34, 0]
fig, ax = plt.subplots(figsize=(6.72, 8.80), dpi=200)
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)
ax.axis("off")
n = len(stages)
h = 1.0
top_w, bot_w = 1.0, 0.42
for i, ((name, val), d) in enumerate(zip(stages, drops)):
    y0 = n - i
    w0 = top_w - (top_w - bot_w) * (i / n)
    w1 = top_w - (top_w - bot_w) * ((i + 1) / n)
    xs = [-w0 / 2, w0 / 2, w1 / 2, -w1 / 2]
    ys = [y0, y0, y0 - h * 0.86, y0 - h * 0.86]
    ax.fill(xs, ys, color=TEAL, alpha=0.35 + 0.08 * i, edgecolor=TEAL_L, lw=1.2)
    ax.text(0, y0 - h * 0.43, "%s | %d" % (name, val), ha="center", va="center",
            color=FG, fontsize=15, weight="bold")
    if d:
        ax.annotate("↓ -%d" % d, (bot_w / 2 + 0.06, y0 - h * 0.93), color=GOLD, fontsize=12,
                    ha="left", va="center")
ax.set_xlim(-0.62, 0.78)
ax.set_ylim(0.1, n + 0.35)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "p16_funnel.png"), facecolor=BG, bbox_inches="tight")
plt.close(fig)
print("charts written to", OUT)
for f in sorted(os.listdir(OUT)):
    print("  ", f, os.path.getsize(os.path.join(OUT, f)))
