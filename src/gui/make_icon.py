# -*- coding: utf-8 -*-
"""生成 GUI 图标 (icon.ico 多尺寸 + app_icon.png + 内嵌 base64 供运行时窗口图标)

设计: 蓝色渐变圆角方块 + 白色苯环分子图案 (六边形 + 节点 + 键)。
运行: python src/gui/make_icon.py
"""
import os, base64
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
ICO = os.path.join(HERE, "icon.ico")
PNG = os.path.join(HERE, "app_icon.png")
B64PY = os.path.join(HERE, "icon_b64.py")

TOP = (59, 130, 246)      # #3B82F6
BOT = (29, 78, 216)       # #1D4ED8
S = 512


def gradient_bg(size):
    img = Image.new("RGB", (size, size))
    d = ImageDraw.Draw(img)
    for y in range(size):
        t = y / max(1, size - 1)
        # 轻微非线性, 上部更亮
        t2 = t ** 0.85
        c = tuple(int(TOP[i] + (BOT[i] - TOP[i]) * t2) for i in range(3))
        d.line([(0, y), (size, y)], fill=c)
    return img


def rounded_mask(size, radius):
    m = Image.new("L", (size, size), 0)
    ImageDraw.Draw(m).rounded_rectangle([0, 0, size - 1, size - 1], radius=radius, fill=255)
    return m


def hexagon(cx, cy, r, rot=-90):
    import math
    return [(cx + r * math.cos(math.radians(rot + 60 * i)),
             cy + r * math.sin(math.radians(rot + 60 * i))) for i in range(6)]


def draw_molecule(d, size):
    cx = cy = size / 2
    r = size * 0.245
    pts = hexagon(cx, cy, r)
    lw = max(2, int(size * 0.030))
    # 外环
    d.polygon(pts, outline=(255, 255, 255), width=lw)
    # 内环 (芳香性提示)
    inner = hexagon(cx, cy, r * 0.60)
    d.polygon(inner, outline=(255, 255, 255, 170), width=max(2, int(lw * 0.62)))
    # 顶点节点
    nr = size * 0.042
    for i, (x, y) in enumerate(pts):
        col = (255, 255, 255) if i % 2 == 0 else (191, 219, 254)
        d.ellipse([x - nr, y - nr, x + nr, y + nr], fill=col)
    # 三条取代基键 + 末端节点
    import math
    for ang, rr in ((30, 1.62), (150, 1.62), (270, 1.70)):
        x1 = cx + r * math.cos(math.radians(ang))
        y1 = cy + r * math.sin(math.radians(ang))
        x2 = cx + r * rr * math.cos(math.radians(ang))
        y2 = cy + r * rr * math.sin(math.radians(ang))
        d.line([x1, y1, x2, y2], fill=(255, 255, 255), width=max(2, int(lw * 0.8)))
        d.ellipse([x2 - nr * 0.78, y2 - nr * 0.78, x2 + nr * 0.78, y2 + nr * 0.78],
                  fill=(147, 197, 253))


def main():
    img = gradient_bg(S)
    draw_molecule(ImageDraw.Draw(img), S)
    mask = rounded_mask(S, int(S * 0.22))
    out = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    out.paste(img, (0, 0), mask)
    # 顶部高光
    hl = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    ImageDraw.Draw(hl).rounded_rectangle([0, 0, S - 1, int(S * 0.46)], radius=int(S * 0.22),
                                         fill=(255, 255, 255, 26))
    out = Image.alpha_composite(out, hl)

    out.resize((512, 512), Image.LANCZOS).save(PNG)
    out.save(ICO, sizes=[(256, 256), (128, 128), (64, 64), (48, 48), (32, 32), (24, 24), (16, 16)])
    # 运行时窗口图标 (PhotoImage 只能用 PNG/GIF; 这里存 64px base64)
    import io
    buf = io.BytesIO()
    out.resize((64, 64), Image.LANCZOS).save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode()
    with open(B64PY, "w", encoding="utf-8") as f:
        f.write("# -*- coding: utf-8 -*-\n# 自动生成, 勿手改 (gui/make_icon.py)\n")
        f.write("ICON_B64 = (\n")
        for i in range(0, len(b64), 96):
            f.write('    "%s"\n' % b64[i:i + 96])
        f.write(")\n")
    print("icon.ico  %d bytes" % os.path.getsize(ICO))
    print("app_icon.png %d bytes" % os.path.getsize(PNG))
    print("icon_b64.py  %d bytes" % os.path.getsize(B64PY))


if __name__ == "__main__":
    main()
