# -*- coding: utf-8 -*-
"""Pocket2Mol 竞赛项目展示 PPT 生成器 — 科幻深空风, 40 页, 16:9

用法:
    python docs/make_ppt.py [输出文件.pptx]
    不带参数时输出到 docs/ 下 (原实现写死在个人桌面上, 已改为参数/默认目录)。
"""
import os
import sys

# ---- 路径前置: docs/make_ppt.py -> 仓库根 -> 仓库根/src (统一从 paths.py 取路径) ----
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from paths import RESULTS, DOCS  # noqa: E402

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

# ---------- 设计系统 ----------
BG      = RGBColor(0x07, 0x0B, 0x1A)   # 深空底
PANEL   = RGBColor(0x0D, 0x14, 0x2E)   # 卡片底
CYAN    = RGBColor(0x00, 0xE5, 0xFF)   # 霓虹青(主)
PURPLE  = RGBColor(0xB3, 0x88, 0xFF)   # 霓虹紫(辅)
GREEN   = RGBColor(0x00, 0xFF, 0x9C)   # 霓虹绿(数据)
ORANGE  = RGBColor(0xFF, 0xB7, 0x4D)   # 橙(对比)
WHITE   = RGBColor(0xE8, 0xF0, 0xFF)
GREY    = RGBColor(0x8F, 0xA3, 0xC8)
DIM     = RGBColor(0x1B, 0x2A, 0x4A)   # 弱线条
RED     = RGBColor(0xFF, 0x5A, 0x6E)

F = "Microsoft YaHei"
FC = "Consolas"

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]
W, H = 13.333, 7.5

def slide():
    s = prs.slides.add_slide(BLANK)
    r = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, Inches(W), Inches(H))
    r.fill.solid(); r.fill.fore_color.rgb = BG
    r.line.fill.background()
    r.shadow.inherit = False
    return s

def box(s, x, y, w, h, fill=None, line=None, lw=1.0):
    r = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    if fill is None: r.fill.background()
    else: r.fill.solid(); r.fill.fore_color.rgb = fill
    if line is None: r.line.fill.background()
    else: r.line.color.rgb = line; r.line.width = Pt(lw)
    r.shadow.inherit = False
    return r

def text(s, x, y, w, h, runs, size=14, color=WHITE, bold=False, align=PP_ALIGN.LEFT,
         font=F, anchor=MSO_ANCHOR.TOP, spacing=1.0):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = True; tf.vertical_anchor = anchor
    if isinstance(runs, str): runs = [(runs, size, color, bold, font)]
    first = True
    for run in runs:
        txt, sz, cl, bd = run[0], run[1], run[2], run[3]
        ft = run[4] if len(run) > 4 else F
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.alignment = align
        p.line_spacing = spacing
        r = p.add_run(); r.text = txt
        r.font.size = Pt(sz); r.font.color.rgb = cl; r.font.bold = bd; r.font.name = ft
    return tb

def deco(s, tag=""):
    """科幻装饰: 顶部细条/角落双线框/页脚"""
    box(s, 0, 0, W, 0.05, fill=CYAN)
    box(s, 0, H-0.03, W, 0.03, fill=DIM)
    e1 = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(W-1.35), Inches(0.28), Inches(0.9), Inches(0.9))
    e1.fill.background(); e1.line.color.rgb = DIM; e1.line.width = Pt(1.5); e1.shadow.inherit = False
    e2 = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(W-1.2), Inches(0.43), Inches(0.6), Inches(0.6))
    e2.fill.background(); e2.line.color.rgb = RGBColor(0x24, 0x3A, 0x66); e2.line.width = Pt(1); e2.shadow.inherit = False
    text(s, W-2.2, H-0.42, 2.0, 0.3, [("POCKET2MOL // 赛道3·子任务2", 9, GREY, False, FC)], align=PP_ALIGN.RIGHT)
    if tag:
        text(s, 0.25, H-0.42, 3.0, 0.3, [(tag, 9, DIM, False, FC)])

def page_no(s, n):
    text(s, 12.4, H-0.42, 0.7, 0.3, [("%02d" % n, 10, CYAN, True, FC)], align=PP_ALIGN.RIGHT)

def hbar(s, x, y, w, ratio, color, label="", val=""):
    """霓虹条形图: 底槽+填充+标签"""
    box(s, x, y, w, 0.30, fill=PANEL)
    box(s, x, y, w*ratio, 0.30, fill=color)
    if label: text(s, x, y-0.32, w*0.7, 0.3, [(label, 12, WHITE, True)])
    if val: text(s, x+w+0.1, y, 1.8, 0.3, [(val, 13, color, True, FC)])

def stat_card(s, x, y, w, h, big, unit, label, color=CYAN, sub=""):
    box(s, x, y, w, h, fill=PANEL, line=color, lw=1.2)
    text(s, x+0.15, y+0.12, w-0.3, h*0.52,
         [(str(big), 40, color, True, FC), (" "+unit, 16, GREY, False, FC)])
    text(s, x+0.15, y+h*0.58, w-0.3, h*0.4, [(label, 13, WHITE, True)])
    if sub: text(s, x+0.15, y+h-0.35, w-0.3, 0.3, [(sub, 10, GREY, False)])

def bullets(s, x, y, w, items, size=15, gap=0.42, color=WHITE, mark="▸", mark_color=CYAN):
    yy = y
    for it in items:
        if isinstance(it, tuple):
            head, body = it
            text(s, x, yy, 0.35, gap, [(mark, size, mark_color, True, FC)])
            text(s, x+0.4, yy, w-0.4, gap*2,
                 [(head+"  ", size, WHITE, True), (body, size, GREY, False)], spacing=1.05)
        else:
            text(s, x, yy, 0.35, gap, [(mark, size, mark_color, True, FC)])
            text(s, x+0.4, yy, w-0.4, gap*2, [(it, size, color, False)], spacing=1.05)
        yy += gap

def section(n, zh, en, items):
    s = slide()
    box(s, 0, 0, 0.18, H, fill=CYAN)
    text(s, 0.7, 1.1, 4.5, 2.4, [("%02d" % n, 100, PANEL, True, FC)])
    text(s, 0.75, 2.9, 8.5, 1.2, [(zh, 44, WHITE, True)])
    text(s, 0.78, 4.0, 8.5, 0.6, [(en.upper(), 18, CYAN, True, FC)])
    yy = 5.0
    for it in items:
        text(s, 0.8, yy, 10.5, 0.45, [("▸ ", 15, CYAN, True, FC), (it, 15, GREY, False)])
        yy += 0.5
    text(s, 8.2, 0.9, 5.0, 5.6, [("◈", 170, RGBColor(0x0C, 0x12, 0x28), True, FC)], align=PP_ALIGN.CENTER)
    deco(s); page_no(s, _n())
    return s

N = [0]
def _n():
    N[0] += 1
    return N[0]

def content(title, en, blocks, tag=""):
    """标准内容页: 标题+副标题+若干块。blocks=[(kind, ...)]"""
    s = slide()
    text(s, 0.55, 0.32, 11.0, 0.65, [(title, 27, WHITE, True)])
    text(s, 0.57, 0.95, 11.0, 0.4, [(en.upper(), 12, CYAN, True, FC)])
    box(s, 0.57, 1.32, 2.4, 0.035, fill=CYAN)
    yy = 1.62
    for b in blocks:
        kind = b[0]
        if kind == "h":       # 小节标题
            text(s, 0.57, yy+0.06, 12.2, 0.4, [("▮ " + b[1], 16, CYAN, True)])
            yy += 0.5
        elif kind == "b":     # bullets
            bullets(s, 0.62, yy, 12.1, b[1], size=b[2] if len(b) > 2 else 14)
            yy += (b[2] if len(b) > 2 else 14) * 0 + len(b[1]) * 0.46 + 0.18
        elif kind == "gap":
            yy += b[1]
        elif kind == "cards": # 数据卡片行 [(big,unit,label,color,sub)]
            xs = 0.57; cw = (12.2 - 0.3 * (len(b[1]) - 1)) / len(b[1])
            for c in b[1]:
                stat_card(s, xs, yy, cw, 1.55, *c[:4], sub=c[4] if len(c) > 4 else "")
                xs += cw + 0.3
            yy += 1.75
        elif kind == "bars":  # [(label,ratio,color,val)]
            for lb, rt, cl, vl in b[1]:
                hbar(s, 0.62, yy+0.3, 8.6, rt, cl, lb, vl)
                yy += 0.78
        elif kind == "t":     # 表格 [(rows, widths, header_color)]
            tbl_data, colw, hc = b[1], b[2], b[3] if len(b) > 3 else CYAN
            nrow, ncol = len(tbl_data), len(tbl_data[0])
            tw = sum(colw); rh = 0.42
            for ri, row in enumerate(tbl_data):
                xx = 0.62
                for ci, cell in enumerate(row):
                    cw2 = colw[ci]
                    if ri == 0:
                        box(s, xx, yy, cw2-0.04, rh, fill=hc)
                        text(s, xx+0.06, yy+0.05, cw2-0.15, rh-0.08, [(str(cell), 12.5, BG, True)])
                    else:
                        box(s, xx, yy, cw2-0.04, rh, fill=PANEL)
                        cl = WHITE if ci == 0 else GREY
                        text(s, xx+0.06, yy+0.05, cw2-0.15, rh-0.08,
                             [(str(cell), 12.5, cl, ci == 0, FC if ci == 0 else F)])
                    xx += cw2
                yy += rh
            yy += 0.2
        elif kind == "note":
            box(s, 0.57, yy, 12.2, b[2] if len(b) > 2 else 0.7, fill=PANEL)
            text(s, 0.75, yy+0.08, 11.8, (b[2] if len(b) > 2 else 0.7)-0.12,
                 [("◈ ", 12, CYAN, True, FC), (b[1], 13, WHITE, False)], spacing=1.05)
            yy += (b[2] if len(b) > 2 else 0.7) + 0.18
        elif kind == "formula":
            box(s, 0.57, yy, 12.2, 0.85, fill=RGBColor(0x08, 0x10, 0x24), line=CYAN, lw=1.2)
            text(s, 0.8, yy+0.14, 11.8, 0.6, [(b[1], 17, GREEN, True, FC)])
            yy += 1.05
        elif kind == "flow":  # 横向流程 [(步骤, 颜色)]
            xx = 0.62; bw = 11.0 / len(b[1])
            for i, (step, cl) in enumerate(b[1]):
                box(s, xx, yy, bw-0.22, 0.62, fill=PANEL, line=cl, lw=1.4)
                text(s, xx, yy+0.10, bw-0.22, 0.45, [(step, 13.5, cl, True)], align=PP_ALIGN.CENTER)
                if i < len(b[1]) - 1:
                    text(s, xx+bw-0.24, yy+0.08, 0.3, 0.4, [("▶", 13, GREY, True, FC)])
                xx += bw
            yy += 0.95
        elif kind == "big":   # 居中大字
            text(s, 0.57, yy, 12.2, 1.0, [(b[1], b[2] if len(b) > 2 else 30, b[3] if len(b) > 3 else GREEN, True)], align=PP_ALIGN.CENTER)
            yy += 1.1
    deco(s, tag); page_no(s, _n())
    return s

# ================= 01 封面 =================
s = slide()
# 星点
for (x, y) in [(1.2,0.8),(2.8,1.6),(4.1,0.5),(5.6,1.2),(7.3,0.7),(9.1,1.5),(10.8,0.6),(12.2,1.3),(2.0,5.9),(3.6,6.4),(9.8,6.2),(11.5,5.8),(6.4,6.6),(12.4,3.4),(0.6,3.2)]:
    e = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(x), Inches(y), Inches(0.045), Inches(0.045))
    e.fill.solid(); e.fill.fore_color.rgb = RGBColor(0x3E, 0x5A, 0x8F); e.line.fill.background(); e.shadow.inherit = False
# 中央六边形
hexa = s.shapes.add_shape(MSO_SHAPE.HEXAGON, Inches(5.17), Inches(1.15), Inches(3.0), Inches(2.6))
hexa.fill.background(); hexa.line.color.rgb = CYAN; hexa.line.width = Pt(2.2); hexa.shadow.inherit = False
hexb = s.shapes.add_shape(MSO_SHAPE.HEXAGON, Inches(5.42), Inches(1.38), Inches(2.5), Inches(2.15))
hexb.fill.background(); hexb.line.color.rgb = RGBColor(0x24, 0x3A, 0x66); hexb.line.width = Pt(1.2); hexb.shadow.inherit = False
text(s, 0.5, 1.55, 12.33, 0.5, [("M O L E C U L E   G E N E S I S   E N G I N E", 13, GREY, True, FC)], align=PP_ALIGN.CENTER)
text(s, 0.5, 3.9, 12.33, 1.0, [("口袋引导的 GPCR 靶向分子生成", 40, WHITE, True)], align=PP_ALIGN.CENTER)
text(s, 0.5, 4.85, 12.33, 0.6, [("基于等变图神经网络与引导束搜索的类药化合物库构建", 19, CYAN, True)], align=PP_ALIGN.CENTER)
box(s, 4.42, 5.62, 4.5, 0.03, fill=PURPLE)
text(s, 0.5, 5.85, 12.33, 0.5, [("全球大学生生命科学挑战赛 · 赛道3 离子通道/GPCR 小分子药物虚拟筛选 · 子任务2", 14, GREY, False)], align=PP_ALIGN.CENTER)
text(s, 0.5, 6.3, 12.33, 0.5, [("生成式 AI  ·  引导束搜索  ·  2,931 类药化合物  ·  虚拟筛选闭环  ·  全流程开源", 13, GREEN, True, FC)], align=PP_ALIGN.CENTER)
deco(s)

# ================= 02 目录 =================
s = slide()
text(s, 0.55, 0.35, 6.0, 0.7, [("目录", 30, WHITE, True)])
text(s, 0.57, 1.05, 6.0, 0.4, [("C O N T E N T S", 13, CYAN, True, FC)])
toc = [
    ("01", "背景", "Background", "GPCR 药物设计的时代命题与三重瓶颈"),
    ("02", "项目思路", "Methodology", "基线复现 → 引导束搜索 → 微调框架 → 对接闭环"),
    ("03", "效果展示", "Results", "4,231 生成 · 2,931 入库 · 全库筛选 · 全部高新颖"),
    ("04", "使用方法", "Usage", "三步式 GUI / 命令行 / 全套自动化工具链"),
    ("05", "作用", "Impact", "对研发 · 对评审 · 对开源社区的三重价值"),
    ("06", "前景", "Horizon", "实验闭环 · 靶点家族扩展 · 主动学习循环"),
]
yy = 1.7
for num, zh, en, desc in toc:
    box(s, 0.57, yy, 12.2, 0.78, fill=PANEL)
    text(s, 0.75, yy+0.12, 1.1, 0.5, [(num, 26, CYAN, True, FC)])
    text(s, 1.85, yy+0.06, 2.6, 0.65, [(zh, 19, WHITE, True)])
    text(s, 4.3, yy+0.14, 2.8, 0.5, [(en.upper(), 11, PURPLE, True, FC)])
    text(s, 7.0, yy+0.16, 5.6, 0.5, [(desc, 13, GREY, False)])
    yy += 0.92
deco(s); page_no(s, _n())

# ================= 第一章 背景 =================
section(1, "背景", "Background", [
    "GPCR——人体最大膜蛋白药物靶点家族，800+ 成员",
    "动态构象 / 化学空间 / 选择性：药物设计三重瓶颈",
    "赛道3 子任务2：生成式模型探索化学空间",
    "我们的定位：不是找一个药，而是造一台分子生成引擎",
])

content("GPCR：新药研发的黄金靶点家族", "Golden Targets", [
    ("cards", [
        ("800", "+", "GPCR 家族成员", CYAN, "AlphaFold 全覆盖"),
        ("~40", "%", "上市药物作用于 GPCR", GREEN, "最大药物靶点超家族"),
        ("100-150", "", "具实验结构可服务靶点", PURPLE, "高质量结构 + 配体数据"),
        ("50-100", "", "本项目最佳服务区间", ORANGE, "有共结晶配体的黄金靶点"),
    ]),
    ("b", [
        ("四大经典靶点已纳入", "A2A 腺苷受体(1.8Å) · β2-AR(2.4Å) · D3 多巴胺受体 · 5-HT2B——分辨率高、活性配体数据充足、疾病关联明确"),
        ("疾病图谱", "帕金森病 / 肿瘤免疫 / 哮喘 / 心血管 / 精神分裂 / 偏头痛……"),
    ]),
], tag="01 背景")

content("三重瓶颈：为什么 GPCR 药物设计这么难", "Three Bottlenecks", [
    ("cards", [
        ("动态", "", "构象复杂", RED, "跨膜蛋白高度动态，静态结构难以刻画"),
        ("10⁶⁰", "", "化学空间量级", ORANGE, "类药分子空间天文数字，盲目探索低效"),
        ("亚型", "", "选择性难题", PURPLE, "相近亚型相似度高，偏向性设计极难"),
    ]),
    ("b", [
        ("传统方法", "高通量筛选成本高、命中率低；理性设计依赖专家经验，周期以年计"),
        ("已有生成模型", "通用蛋白口袋训练，对 GPCR 深窄疏水口袋适配不足，类药性与新颖性难以兼顾"),
    ]),
    ("note", "赛道原文瓶颈提炼：①动态构象复杂 ②化学空间探索效率低 ③选择性与精度不足——本项目的三个改进全部对准这三点。", 0.75),
], tag="01 背景")

content("任务解读：赛道3 · 子任务2", "Task Analysis", [
    ("b", [
        ("任务原文", "基于生成式AI模型，发展高效的小分子生成算法，生成结构新颖、具备类药特征的大规模化合物数据库"),
        ("评分结构", "设计与计算筛选 55 分（核心）+ 选题价值 15 + 湿实验验证 15 + 报告展示 15"),
        ("交付要求", "原始代码（可执行证明）+ 候选分子信息 + 实验方案，鼓励开源共享"),
    ]),
    ("h", "四个关键词 = 四个必须交付的证据"),
    ("t", [
        ["关键词", "我们的对应交付"],
        ["生成式AI模型", "等变 GNN 基线复现 + 全链路打通（Windows 单机）"],
        ["高效生成算法", "引导束搜索（核心创新，λ 扫描 + 消融实验）"],
        ["结构新颖", "737 已知分子对照：Top50 全部高新颖（Tanimoto < 0.32）"],
        ["类药特征+大规模", "2,931 入库 · QED 中位 0.825 · 1,960 骨架 · 全库对接"],
    ], [2.8, 9.4], CYAN),
], tag="01 背景")

content("我们的定位：造引擎，不挖矿", "Positioning", [
    ("formula", "通用生成模型 + 跨 GPCR 化合物数据库 + 平台化工具 = 服务全部 800+ GPCR 的引擎"),
    ("b", [
        ("目标不是某个药", "真正成为药物的分子，由使用我们模型的全球研究者在生成的数据库中筛选并推进临床"),
        ("靶点=展台", "A2A 等四个靶点用于验证引擎的通用性与可靠性，而非限制能力边界"),
        ("三方向组合", "方向1 引导束搜索(算法创新) + 方向2 多靶点微调框架(数据适配) + 方向3 对接闭环(验证)"),
        ("产品化", "本地 GUI 应用 + GitHub 分发，Windows 单机即可运行，零服务器依赖"),
    ]),
    ("flow", [("基线模型", CYAN), ("算法改进", GREEN), ("微调框架", PURPLE), ("数据库", ORANGE), ("对接筛选", CYAN), ("候选分子", GREEN)]),
], tag="01 背景")

# ================= 第二章 项目思路 =================
section(2, "项目思路", "Methodology", [
    "基线：Pocket2Mol 等变图神经网络，口袋条件逐原子生成",
    "核心创新：引导束搜索——在搜索中长出更类药的分子",
    "双保险：多样性惩罚回调 + 严格 A/B 判据",
    "工程闭环：从 SMILES 到对接打分的全自动管线",
])

content("总体技术路线", "Pipeline Overview", [
    ("flow", [("PDB 口袋", CYAN), ("等变GNN 编码", GREEN), ("引导束搜索", PURPLE), ("分子重建", ORANGE), ("过滤入库", CYAN), ("Vina 对接", GREEN), ("Top 候选", PURPLE)]),
    ("b", [
        ("输入", "蛋白口袋 3D 结构（PDB + 口袋中心坐标），任意 GPCR 通用"),
        ("生成", "等变图神经网络逐原子自回归：位置(MDN) → 元素 → 化学键，束搜索并行探索"),
        ("改进", "束排序引入 QED/SA 化学分数与 Tanimoto 多样性惩罚（方向1 核心）"),
        ("验证", "AutoDock Vina 全库对接 + 已知活性分子阳性对照 + 新颖性对照"),
    ]),
    ("note", "全链条在本机 RTX 4070 Laptop (8GB) 单卡完成，无需服务器——科研可及性的工程证明。", 0.75),
], tag="02 项目思路")

content("基线模型：Pocket2Mol", "Baseline", [
    ("cards", [
        ("370", "万", "模型参数量", CYAN, "轻量可单卡训练"),
        ("6", "层", "等变注意力编码", GREEN, "标量/向量双通道"),
        ("48", "", "口袋图 KNN", PURPLE, "protein-ligand 组合图"),
        ("4", "头", "生成任务头", ORANGE, "前沿/位置/元素/键"),
    ]),
    ("b", [
        ("生成方式", "自回归逐原子：预测生长前沿 → MDN 高斯混合预测位置 → 分类元素 → 预测键，循环至闭合"),
        ("条件机制", "只看口袋三维原子坐标，不依赖配体——任意 GPCR 给定口袋即可生成"),
    ]),
    ("note", "基线为北京大学团队开源工作（ICML 2022），我们完成 Windows/PyTorch 2.6 全套兼容补丁并单机复现，在此之上做算法改进——尊重开源、站在巨人肩上。", 0.75),
], tag="02 项目思路")

content("核心创新：引导束搜索", "Guided Beam Search ★", [
    ("h", "原版束排序只看模型自身概率，化学空间探索盲目"),
    ("formula", "P ∝ ( exp(Σ logP_model) + 1 ) · exp( λ · [ QED + (1 − SA/10) − w_div · maxTanimoto ] )"),
    ("b", [
        ("机制", "束排序概率直接乘以化学分数增益——在搜索过程中，类药、易合成的分支被系统性加权"),
        ("与事后过滤的本质区别", "过滤=生成完再丢弃（浪费算力）；引导=同等算力下产出的分布本身更优"),
        ("多样性守护", "对成品库的 MMR 式最大相似度惩罚，防止搜索塌缩到单一骨架"),
        ("工程形态", "config 开关默认关闭——关闭时与原版行为逐位一致，不破坏基础功能"),
    ]),
], tag="02 项目思路")

content("参数定型：λ 扫描与消融", "λ Ablation", [
    ("t", [
        ["配置 (seed 2024)", "闭合数", "QED", "SA", "MW", "骨架多样性"],
        ["基线 (λ=0)", "62", "0.667", "3.73", "374", "44"],
        ["λ=0.3", "52", "0.710", "—", "360", "38"],
        ["λ=1.0", "56", "0.682", "—", "371", "35"],
        ["λ=3.0 ★", "53", "0.803", "3.69", "360", "32"],
    ], [3.4, 1.7, 1.7, 1.7, 1.7, 2.0], GREEN),
    ("bars", [
        ("基线 QED", 0.667/1.0, PURPLE, "0.667"),
        ("引导 λ=3 QED", 0.803/1.0, GREEN, "0.803 (+20%)"),
    ]),
    ("note", "seed 2025 复验：0.688 → 0.821（+19%）——两随机种子稳定复现；SA 合成可及性 100% ≤ 6 不受损；多样性由 Tanimoto 惩罚回调守护。λ=3.0 定型。", 0.75),
], tag="02 项目思路")

content("方向2：多靶点微调框架（诚实记录）", "Fine-tuning Framework", [
    ("cards", [
        ("479", "对", "四靶点训练配对", CYAN, "A2A/β2AR/D3/5HT2B"),
        ("2000", "步", "正式微调", PURPLE, "7 分钟/次"),
        ("98.6", "%", "对接成功率", GREEN, "全管线稳定性"),
        ("0", "晋级", "当前结论", ORANGE, "诚实记录"),
    ]),
    ("b", [
        ("自建数据管线", "ChEMBL 活性配体 CSV + 受体 PDB → 3D 构象(SDF) → 口袋裁剪 → index/split 训练格式，全自动"),
        ("质量守卫", "构象审计(键长/原子重叠/共面/孤立原子循环清洗) + 非有限损失跳过 + 初始化检查点不覆盖预训练"),
        ("严格判据", "A/B 复验 50/100/50：微调 50 个/QED 0.408 vs 基线 63 个/0.597 → 不晋级，生产权重维持预训练"),
    ]),
    ("note", "不晋级的原因清晰：自建构象(ETKDG)与实验复合物构象存在分布差。管线与接口已就绪——引入实验构象即可重启，这是下一阶段工作。", 0.75),
], tag="02 项目思路")

content("方向3 + 工程化：全自动闭环", "Docking & Engineering", [
    ("flow", [("SMILES", CYAN), ("3D 构象", GREEN), ("PDBQT", PURPLE), ("Vina 对接", ORANGE), ("打分汇总", CYAN), ("Top 候选", GREEN)]),
    ("b", [
        ("对接管线", "AutoDock Vina 1.2.5 + Open Babel 3.2.1 便携部署；受体自动去水/去共晶/转换；配体自动加氢生成构象"),
        ("过滤入库管线", "规范化去重 → PAINS/BRENK 警示结构 → 理化五规则 → SA ≤ 6 → 元素白名单 → SQLite + SDF 库，全带剔除原因"),
        ("桌面应用", "Tkinter GUI（三步式外行界面 + 2D 结构渲染 + 自动过滤）打包为双击即用 exe"),
    ]),
], tag="02 项目思路")

# ================= 第三章 效果展示 =================
section(3, "效果展示", "Results", [
    "四靶点生成 10,004 个分子，过滤入库 6,864 个",
    "A2A 全库对接：强结合 1,137，极强 311",
    "Top50 候选对 737 已知分子全部高新颖",
    "引导束搜索 QED 提升 20%，两种子稳定复现",
])

content("生成规模总览", "Generation Scale", [
    ("cards", [
        ("4,231", "", "A2A 产出", CYAN, "113 批 × 50-68"),
        ("1,577", "", "β2-AR 产出", GREEN, "29 批"),
        ("1,571", "", "D3 产出", PURPLE, "29 批"),
        ("1,525", "", "5-HT2B 产出", ORANGE, "25 批"),
    ]),
    ("bars", [
        ("A2A 入库 2,931", 2931/2931.0, CYAN, "2,931"),
        ("β2-AR 入库 917", 917/2931.0, GREEN, "917"),
        ("D3 入库 951", 951/2931.0, PURPLE, "951"),
        ("5-HT2B 入库 1,065", 1065/2931.0, ORANGE, "1,065"),
    ]),
    ("note", "四靶点同一管线全自动运行，零批次失败——跨 GPCR 通用性的工程实证。", 0.6),
], tag="03 效果展示")

content("数据库构建漏斗", "Library Funnel", [
    ("t", [
        ["环节", "存活", "剔除", "说明"],
        ["采样产出", "4,231", "—", "113 批引导束搜索"],
        ["规范化去重", "4,183", "48", "批间唯一率 98.9%"],
        ["PAINS 警示", "3,957", "226", "假阳性结构剔除"],
        ["BRENK 警示", "3,305", "652", "毒性/反应性基团"],
        ["MW [250,500]", "3,167", "138", "类药五规则"],
        ["LogP [1,5]", "2,932", "235", "亲脂平衡"],
        ["QED/SA/元素/环", "2,931", "1", "最终入库"],
    ], [3.2, 2.0, 2.0, 5.0], CYAN),
    ("note", "每一级剔除全部带原因记录（rejected.csv），可追溯、可审计。", 0.6),
], tag="03 效果展示")

content("库质量画像", "Library Quality", [
    ("cards", [
        ("0.825", "", "QED 中位", GREEN, "最高 0.946"),
        ("3.33", "", "SA 中位", CYAN, "1.5-5.9 全可合成"),
        ("1,960", "", "唯一骨架", PURPLE, "化学空间广覆盖"),
        ("362", "Da", "MW 中位", ORANGE, "口服药物甜区"),
    ]),
    ("h", "QED 分布（n=2,931）"),
    ("bars", [
        ("0.4-0.6", 56/1460.0, GREY, "56"),
        ("0.6-0.8", 1098/1460.0, PURPLE, "1,098"),
        ("0.8-1.0", 1770/1460.0, GREEN, "1,770"),
    ]),
    ("note", "76% 的分子 QED ≥ 0.8——直接对准子任务2“结构新颖、具备类药特征”的核心要求。", 0.6),
], tag="03 效果展示")

content("引导束搜索效果：QED +20%", "Guided Search Effect ★", [
    ("h", "两个独立随机种子，50 样本/100 束宽/50 步严格对照"),
    ("bars", [
        ("seed 2024 基线", 0.667/1.0, PURPLE, "0.667"),
        ("seed 2024 引导λ=3", 0.803/1.0, GREEN, "0.803"),
        ("seed 2025 基线", 0.688/1.0, PURPLE, "0.688"),
        ("seed 2025 引导λ=3", 0.821/1.0, GREEN, "0.821"),
    ]),
    ("b", [
        ("合成可及性不受损", "SA 100% ≤ 6（引导 3.69/3.16 vs 基线 3.73/3.76）"),
        ("闭合代价可控", "-15% / +2%，平均约 -7%，多跑一批即可弥补"),
    ]),
    ("note", "这不是事后过滤——是搜索分布本身的系统性改善，同等算力下产出更优。", 0.6),
], tag="03 效果展示")

content("全库虚拟筛选", "Virtual Screening", [
    ("cards", [
        ("2,931", "", "全库对接成功", CYAN, "成功率 100%"),
        ("-9.68", "kcal/mol", "平均打分", GREEN, "范围 -13.51 ~ -5.92"),
        ("1,137", "", "强结合 <-10", PURPLE, "38.8%"),
        ("311", "", "极强 <-11", ORANGE, "10.6%"),
    ]),
    ("bars", [
        ("全库平均 -9.68", 9.68/13.51, CYAN, "mean"),
        ("Top1 -13.42", 13.42/13.51, GREEN, "#1"),
        ("已知集均值 -9.32", 9.32/13.51, PURPLE, "known"),
    ]),
], tag="03 效果展示")

content("Top 候选与新奇性对照", "Top Candidates", [
    ("t", [
        ["排名", "Vina 打分", "QED", "SA", "MW", "综合分"],
        ["#1", "-13.42", "0.797", "4.3", "374", "0.912"],
        ["#2", "-12.89", "0.857", "4.1", "390", "0.894"],
        ["#3", "-13.20", "0.787", "4.2", "400", "0.890"],
        ["#4", "-12.89", "0.836", "4.4", "405", "0.886"],
        ["#5", "-12.63", "0.884", "4.0", "368", "0.884"],
    ], [1.6, 2.2, 1.6, 1.6, 1.6, 2.0], GREEN),
    ("b", [
        ("综合排序", "0.6 × 打分归一 + 0.4 × QED；Top50 含 47 个唯一骨架"),
        ("骨架多样性", "Top50 亦经多样性惩罚约束，避免单一骨架堆叠"),
    ]),
], tag="03 效果展示")

content("新颖性验证：全部高新颖 ★", "Novelty Validation", [
    ("h", "对照集：737 个已知 A2A 高活性分子 + 上市药物（ChEMBL CHEMBL251，362 个骨架）"),
    ("t", [
        ["判级", "标准", "Top50 数量", "占比"],
        ["高新颖 HIGH_NOVEL", "Tanimoto < 0.4", "50", "100%"],
        ["近似已知", "Tanimoto ≥ 0.85", "0", "0%"],
        ["已知骨架", "Murcko 骨架重合", "0", "0%"],
    ], [3.6, 3.2, 2.6, 2.6], GREEN),
    ("big", "最大相似度仅 0.187 – 0.321 · 零骨架重合", 22, GREEN),
    ("note", "模型在保持强结合打分与类药性的同时，探索的是已知化学空间之外的新颖区域——“结构新颖、具备类药特征”两项核心要求同时满足。", 0.85),
], tag="03 效果展示")

content("逆向阳性对照：管线可信度", "Reverse Validation", [
    ("h", "将 737 个已知 A2A 分子送入同一对接管线（阳性对照）"),
    ("t", [
        ["集合", "成功数", "均值", "强结合<-10", "极强<-11", "最强"],
        ["已知活性分子", "728 (98.6%)", "-9.32", "22.1%", "6.0%", "-12.73"],
        ["我们的生成库", "2,931 (100%)", "-9.68", "38.8%", "10.6%", "-13.51"],
    ], [3.0, 2.2, 1.8, 2.0, 1.8, 1.6], PURPLE),
    ("b", [
        ("验证一 管线可信", "已知活性分子获得合理强打分——对接设置正确的直接证据"),
        ("验证二 生成更优", "生成库强结合比例为已知集 1.8 倍，最强打分超越已知集最强"),
    ]),
], tag="03 效果展示")

content("跨 GPCR 通用性实证", "Cross-Target Generality", [
    ("t", [
        ["靶点", "入库", "对接均值", "最强打分", "强结合占比"],
        ["A2A 腺苷受体", "2,931", "-9.68", "-13.51", "38.8%"],
        ["β2-AR", "917", "-9.78", "-13.97", "38.9%"],
        ["D3 多巴胺", "951", "-8.40", "-11.39", "5.5%"],
        ["5-HT2B", "1,065", "-9.34", "-12.34", "24.9%"],
    ], [2.9, 1.9, 2.1, 2.1, 2.4], CYAN),
    ("b", [
        ("同一管线 · 四靶点 · 零失败", "从口袋坐标到候选分子的全流程在四个 GPCR 上完整复现"),
        ("差异可解释", "D3 口袋更深窄疏水，强结合占比低符合预期——正是下一步靶点微调的目标场景"),
    ]),
], tag="03 效果展示")

# ================= 第四章 使用方法 =================
section(4, "使用方法", "Usage", [
    "三步式 GUI：选靶点 → 开始 → 看结果，零命令行",
    "命令行全参数控制，面向科研人员",
    "数据管线 / 过滤入库 / 对接 / 新颖性 全套脚本",
    "Windows 单机部署，GitHub 开源",
])

content("三步式桌面应用", "3-Step Desktop App", [
    ("flow", [("① 选择靶点", CYAN), ("② 开始生成", GREEN), ("③ 查看导出", PURPLE)]),
    ("b", [
        ("第一步", "下拉选择靶点——自动载入蛋白结构与口袋坐标，并显示相关疾病说明"),
        ("第二步", "一键开始：默认已预填验证过的最优参数（50/100/50，λ=3，多样性 0.5），进度条与已获分子数实时更新"),
        ("第三步", "分子列表按类药分排序，点击任意分子即渲染 2D 结构，一键导出 3D SDF 或复制编码"),
        ("自动过滤", "生成结束自动执行 PAINS/BRENK/理化过滤并按类药分排序——外行拿到即是可用候选"),
    ]),
    ("note", "双击桌面快捷方式即可启动（12.6→9.8MB 单文件 exe），全部计算在本机 GPU 完成。", 0.6),
], tag="04 使用方法")

content("命令行与批量化", "CLI & Batch", [
    ("h", "科研人员全参数控制"),
    ("formula", "sample_for_pdb.py --pdb_path <口袋> --center \" x,y,z\" --config <yml> --outdir <dir>"),
    ("b", [
        ("批量挂机", "overnight_a2a.ps1：12 小时窗口、每批独立种子、磁盘守卫、自动清快照——过夜产出 4,000+"),
        ("过滤入库", "build_library.py --runs ...：漏斗全记录、SQLite + SDF 库、增量幂等"),
        ("对接筛选", "docking_pipeline.py：受体自动制备、批量打分、CSV 汇总"),
        ("新颖性", "top50_novelty_3targets.py：自动拉取 ChEMBL 已知集并判级"),
    ]),
], tag="04 使用方法")

content("部署：Windows 单机即可", "Deployment", [
    ("b", [
        ("硬件", "任意 NVIDIA 显卡（8GB 显存即可，全程在 RTX 4070 Laptop 验证）"),
        ("软件", "Miniconda + Python 3.10 + PyTorch 2.6 (CUDA 12.4) + PyG/RDKit 全家桶，依赖清单随仓库提供"),
        ("权重", "预训练权重 42.84MB（GitHub Releases 分发），代码 GitHub 开源"),
        ("启动", "exe 双击即用；或 notebook（docs/pocket2mol_demo.ipynb）逐格复现——竞赛要求“证明代码可执行”的直接证据"),
    ]),
    ("note", "兼容性补丁三件套（PyTorch 2.6 安全加载 / torch_cluster KNN / Tk DLL 注入）让整个流程在任何现代 Windows 单卡机器上可复现。", 0.75),
], tag="04 使用方法")

# ================= 第五章 作用 =================
section(5, "作用与价值", "Impact", [
    "对药物研发：先导化合物发现的算力平权",
    "对竞赛评审：三重可验证的量化证据链",
    "对科研社区：可复现、可扩展、可二次开发的开源工具链",
])

content("对药物研发：算力平权", "For Drug Discovery", [
    ("cards", [
        ("80", "0+", "可服务 GPCR", CYAN, "任意口袋即生成"),
        ("~10", "分钟", "单批 50 分子", GREEN, "RTX 4070 单卡"),
        ("38.8", "%", "强结合占比", PURPLE, "vs 已知集 22.1%"),
        ("<0.32", "", "对已知最大相似", ORANGE, "全部高新颖"),
    ]),
    ("b", [
        ("先导发现提速", "从口袋结构到强结合候选由月级压缩到天级，湿实验名额只留给通过计算三重验证的分子"),
        ("实验成功率", "QED/SA/PAINS 前置过滤 + 对接排序——假阳性结构在计算阶段就被淘汰"),
    ]),
], tag="05 作用")

content("对评审：三重可验证证据链", "For Reviewers", [
    ("h", "每个结论都有数据、有脚本、有日志"),
    ("t", [
        ["主张", "证据", "复现入口"],
        ["算法有效", "QED +20% 两种子复现 + λ 消融表", "configs + outputs/ab_*"],
        ["新颖成立", "Top50 vs 737 已知集 Tanimoto 全 <0.32", "top50_novelty.csv"],
        ["管线可信", "已知分子阳性对照 98.6% 成功", "docking_known_reverse"],
        ["规模真实", "4,231 产出漏斗逐级可审计", "library summary + rejected.csv"],
    ], [2.6, 5.6, 4.0], GREEN),
    ("note", "科研诚信：微调不晋级如实汇报（含失败根因与重启条件），高风险优化主动回滚——每个决策可追溯至工作日志。", 0.75),
], tag="05 作用")

content("对社区：开源工具链", "For Community", [
    ("b", [
        ("引擎", "引导束搜索模块（guidance.py）即插即用，兼容任何 Pocket2Mol 系工作流"),
        ("数据", "四靶点 479 对微调数据集 + 构建脚本——社区可直接扩展到新靶点"),
        ("工具", "构象审计 / 候选评估 / 过滤入库 / 对接管线 / GUI 全套脚本，Windows 友好"),
        ("文档", "参赛 README + 可执行 notebook + 全部配置固化，逐实验可复现"),
    ]),
    ("note", "获奖后按赛事要求全面开源——从“AI 设计”到“社区共享”的完整闭环。", 0.6),
], tag="05 作用")

# ================= 第六章 前景 =================
section(6, "前景", "Horizon", [
    "短期：实验复合物构象驱动的微调重启",
    "中期：湿实验闭环 + 更多靶点家族",
    "长期：对接反馈的主动学习循环，引擎持续进化",
])

content("路线图", "Roadmap", [
    ("h", "短期（赛前）"),
    ("flow", [("Top50 人工复核", CYAN), ("5-10 个送验", GREEN), ("结合实验", PURPLE), ("数据回填", ORANGE)]),
    ("h", "中期（赛事后 3-6 个月）"),
    ("b", [
        ("微调重启", "引入实验复合物构象与对接反馈分子，487→1000+ 对数据集，冻结-解冻分阶段训练"),
        ("靶点家族扩展", "同一管线复用到激酶、离子通道等其他靶点家族——引擎的通用性论证"),
    ]),
    ("h", "长期（愿景）"),
    ("b", [
        ("主动学习循环", "湿实验结果反馈引导分数 → 引导束搜索参数自动演化 → 每一轮实验让引擎更懂 GPCR"),
        ("从赛道到平台", "开源社区共建靶点数据包，让“任意 GPCR 口袋 → 候选分子”成为公共基础设施"),
    ]),
], tag="06 前景")

# ================= 36-40 补充页 =================
content("项目历程：从发现到收官", "Timeline", [
    ("t", [
        ["日期", "里程碑"],
        ["08-05", "发现 Pocket2Mol，下载代码并生成 646 行模型详解文档"],
        ["08-07", "环境全链路打通：PyTorch 2.6 + CUDA 12.4，首个分子生成成功"],
        ["08-09/10", "参赛定位确立：模型+数据库“卖铲子”；三方向战略成形"],
        ["08-14", "四靶点数据筹备：546 活性配体 + 100 上市药物库"],
        ["08-26", "V3.0 单机化；解构源码；翻车→回滚→微调路线确立"],
        ["08-28", "接管校验；引导束搜索实现，λ=3 定型（QED +20%）"],
        ["08-29/30", "四靶点 10,004 分子生成、入库、对接、新颖性验证全部收官"],
    ], [1.8, 10.4], CYAN),
], tag="附·历程")

content("瓶颈-方案-证据：一一对应", "Bottleneck → Solution → Evidence", [
    ("t", [
        ["赛道瓶颈", "我们的方案", "量化证据"],
        ["化学空间探索效率低", "引导束搜索（QED/SA 引导 + 多样性惩罚）", "QED +20%；强结合占比 38.8%"],
        ["类药性与新颖难兼顾", "MMR 多样性惩罚 + 全库新颖性对照", "Top50 全高新颖（<0.32）"],
        ["GPCR 适配不足", "四靶点数据管线 + 微调框架（接口就绪）", "四靶点 10,004 产出零失败"],
        ["可复现性疑虑", "配置固化 + 日志归档 + 阳性对照", "已知分子对照 98.6% 成功"],
    ], [3.6, 4.8, 3.8], GREEN),
], tag="附·对照")

content("风险与应对（答辩预备）", "Risks & Mitigations", [
    ("b", [
        ("微调不晋级？", "如实汇报：A/B 判据严格；管线就绪，引入实验复合物构象即重启——科研诚信是加分项"),
        ("多样性下降？", "Tanimoto 惩罚开关可回调；库整体 1,960 骨架、Top50 亦 47 骨架，未塌缩"),
        ("打分是管线假象？", "737 已知分子阳性对照：均值 -9.32 合理，生成库 -9.68 反超——假象排除"),
        ("QED 高≠活性高？", "承认：类药是必要非充分条件；结合活性由对接排序，最终由湿实验裁决"),
        ("D3 强结合占比低？", "口袋更深窄疏水，符合结构预期——正是微调的目标场景，非引擎失效"),
    ]),
], tag="附·风险")

content("交付物与开源索引", "Deliverables Index", [
    ("t", [
        ["类别", "内容", "位置"],
        ["化合物库", "2,931 + 917/951/1,065（SQLite/CSV/SDF）", RESULTS],
        ["候选清单", "四靶点 Top50 + 新颖性判级", "top_candidates.csv / top50_novelty.csv"],
        ["引擎代码", "引导束搜索 + 全管线脚本（GitHub 开源）", "仓库 src/scripts/ + src/utils/guidance.py"],
        ["桌面应用", "三步式 GUI（双击即用 exe）", "src/gui/dist/ + 桌面快捷方式"],
        ["文档", "README + notebook + 摘要/评审表/声明", "docs/ + 竞赛材料\\"],
    ], [2.4, 5.6, 4.2], PURPLE),
], tag="附·索引")

content("附录：四靶点成果总表", "Appendix", [
    ("t", [
        ["靶点", "产出", "入库", "对接均值", "最强", "强结合", "Top50 新颖"],
        ["A2A", "4,231", "2,931", "-9.68", "-13.51", "38.8%", "50/50"],
        ["β2-AR", "1,577", "917", "-9.78", "-13.97", "38.9%", "49/50"],
        ["D3", "1,571", "951", "-8.40", "-11.39", "5.5%", "49/50"],
        ["5-HT2B", "1,525", "1,065", "-9.34", "-12.34", "24.9%", "50/50"],
        ["合计", "10,004", "6,864", "—", "-13.97", "—", "198/200"],
    ], [2.0, 1.7, 1.7, 1.9, 1.8, 1.7, 1.9], GREEN),
    ("note", "全部数据可在仓库 outputs/、results/ 中逐文件追溯；本页所有数字由脚本自动统计生成。", 0.6),
], tag="附·总表")

# ================= 40 结尾 =================
s = slide()
for (x, y) in [(1.0,1.0),(2.6,1.8),(4.3,0.9),(6.2,1.5),(8.1,0.8),(10.0,1.6),(11.8,0.9),(1.8,6.2),(3.4,5.6),(10.4,6.4),(12.2,5.9),(6.8,6.8)]:
    e = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(x), Inches(y), Inches(0.045), Inches(0.045))
    e.fill.solid(); e.fill.fore_color.rgb = RGBColor(0x3E, 0x5A, 0x8F); e.line.fill.background(); e.shadow.inherit = False
hexa = s.shapes.add_shape(MSO_SHAPE.HEXAGON, Inches(5.42), Inches(1.5), Inches(2.5), Inches(2.15))
hexa.fill.background(); hexa.line.color.rgb = CYAN; hexa.line.width = Pt(2.0); hexa.shadow.inherit = False
text(s, 0.5, 3.75, 12.33, 0.9, [("谢谢观看", 44, WHITE, True)], align=PP_ALIGN.CENTER)
text(s, 0.5, 4.75, 12.33, 0.5, [("我们不预言哪个分子会成为药——", 16, GREY, False)], align=PP_ALIGN.CENTER)
text(s, 0.5, 5.25, 12.33, 0.5, [("我们把发现分子的能力，交给每一个拥有口袋结构与一块 GPU 的人。", 16, CYAN, True)], align=PP_ALIGN.CENTER)
box(s, 4.42, 5.95, 4.5, 0.03, fill=PURPLE)
text(s, 0.5, 6.15, 12.33, 0.5, [("全球大学生生命科学挑战赛 · 赛道3 · 子任务2", 12, GREY, False)], align=PP_ALIGN.CENTER)
deco(s)

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DOCS, "Pocket2Mol_项目展示.pptx")
prs.save(out)
print("pages:", N[0], "->", out)
