# -*- coding: utf-8 -*-
"""生成两张文字小图（原为设计工具烘焙的文字图片）+ 紧凑版漏斗图。"""
import os, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"D:\MMModel\outputs\night_pipeline_v10\charts"
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

def chip(text, px_w, px_h, fontsize, out, color="#FFFFFF", weight="bold"):
    fig = plt.figure(figsize=(px_w / 100, px_h / 100), dpi=300)
    fig.patch.set_alpha(0.0)
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.text(0.0, 0.5, text, color=color, fontsize=fontsize, weight=weight,
            va="center", ha="left", transform=ax.transAxes)
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print("chip:", out, text)

chip("四口袋配对", 2.07, 0.46, 11, os.path.join(OUT, "p10_chip.png"))
chip("同一管线 · 三靶点 · 零失败", 5.09, 0.45, 12, os.path.join(OUT, "p22_chip.png"))

# 紧凑漏斗（7 级，含每级剔除数），画布 7.01x9.18cm
BG = "#0C0F1A"; FG = "#E5E7EB"; TEAL = "#1E7A70"; TEAL_L = "#6BC5BD"; GOLD = "#D9A441"
stages = [("采样产出", 583), ("规范化去重", 576), ("PAINS 警示", 560), ("BRENK 警示", 478),
          ("MW [250,500]", 429), ("LogP [1,5]", 395), ("QED/SA/元素/环", 395)]
drops = [None, 7, 16, 82, 49, 34, 0]
fig, ax = plt.subplots(figsize=(7.01 / 2.54, 9.18 / 2.54), dpi=200)
fig.patch.set_facecolor(BG); ax.set_facecolor(BG); ax.axis("off")
n = len(stages)
slot = 1.0                      # 每级占 1 个单位高度
bar_h = 0.66                    # 梯形高度
top_w, bot_w = 1.0, 0.60
for i, ((name, val), d) in enumerate(zip(stages, drops)):
    y0 = (n - i) * slot
    w0 = top_w - (top_w - bot_w) * (i / n)
    w1 = top_w - (top_w - bot_w) * ((i + 1) / n)
    ax.fill([-w0/2, w0/2, w1/2, -w1/2], [y0, y0, y0-bar_h, y0-bar_h],
            color=TEAL, alpha=0.30 + 0.07*i, edgecolor=TEAL_L, lw=1.0)
    ax.text(0, y0 - bar_h/2, "%s | %d" % (name, val), ha="center", va="center",
            color=FG, fontsize=10.5, weight="bold")
    if d:
        ax.text(0, y0 - bar_h - (slot - bar_h)/2, "↓ 剔除 %d" % d, ha="center", va="center",
                color=GOLD, fontsize=9)
ax.set_xlim(-0.78, 0.78)
ax.set_ylim(0.25, n + 0.55)
fig.subplots_adjust(left=0.02, right=0.98, top=0.99, bottom=0.01)
fig.savefig(os.path.join(OUT, "p16_funnel.png"), facecolor=BG)
plt.close(fig)
print("funnel regenerated")
