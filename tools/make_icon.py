"""Draw the app icon (a shaded ball-and-stick molecule) and save it as assets/icon.ico.

Run with:  python tools/make_icon.py      (needs Pillow)
"""

import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter

ROOT = Path(__file__).resolve().parent.parent
S = 1024  # draw large, downsample for crisp small sizes


def shaded_ball(img, cx, cy, r, rgb):
    """Sphere with a soft highlight toward the upper-left, like the viewer's lighting."""
    ball = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    px = ball.load()
    lx, ly, lz = -0.45, -0.55, 0.70  # light direction
    ln = math.sqrt(lx * lx + ly * ly + lz * lz)
    lx, ly, lz = lx / ln, ly / ln, lz / ln
    for y in range(int(cy - r) - 1, int(cy + r) + 2):
        for x in range(int(cx - r) - 1, int(cx + r) + 2):
            dx, dy = (x - cx) / r, (y - cy) / r
            d2 = dx * dx + dy * dy
            if d2 > 1:
                continue
            dz = math.sqrt(1 - d2)
            diff = max(0.0, dx * lx + dy * ly + dz * lz)
            # specular (Blinn-ish, view along +z)
            hx, hy, hz = lx, ly, lz + 1
            hn = math.sqrt(hx * hx + hy * hy + hz * hz)
            spec = max(0.0, (dx * hx + dy * hy + dz * hz) / hn) ** 40
            shade = 0.28 + 0.72 * diff
            c = [min(255, int(ch * shade + 255 * 0.55 * spec)) for ch in rgb]
            edge = min(1.0, (1 - math.sqrt(d2)) * r / 1.5)  # anti-aliased rim
            px[x, y] = (*c, int(255 * edge))
    img.alpha_composite(ball)


def stick(draw, a, b, w, rgb):
    draw.line([a, b], fill=(*rgb, 255), width=w)


def main():
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))

    # rounded-square background with a subtle radial glow
    bg = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    g = bg.load()
    for y in range(S):
        for x in range(S):
            d = math.hypot(x - S * 0.55, y - S * 0.42) / S
            t = min(1.0, d * 1.5)
            g[x, y] = (int(36 - 26 * t), int(46 - 33 * t), int(66 - 46 * t), 255)
    mask = Image.new("L", (S, S), 0)
    ImageDraw.Draw(mask).rounded_rectangle([24, 24, S - 24, S - 24], radius=210, fill=255)
    img.paste(bg, (0, 0), mask)

    # molecule layout: one centre atom bonded to three neighbours
    C = (S * 0.50, S * 0.53)
    atoms = [
        ((S * 0.24, S * 0.30), 118, (255, 255, 255)),  # H-ish white
        ((S * 0.78, S * 0.27), 132, (70, 110, 250)),  # N blue
        ((S * 0.62, S * 0.82), 125, (255, 40, 40)),  # O red
    ]
    sticks = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    sd = ImageDraw.Draw(sticks)
    for (p, _, rgb) in atoms:
        mid = ((C[0] + p[0]) / 2, (C[1] + p[1]) / 2)
        stick(sd, C, mid, 64, (150, 150, 150))
        stick(sd, mid, p, 64, tuple(int(v * 0.85) for v in rgb))
    # soften stick edges a touch
    img.alpha_composite(sticks.filter(ImageFilter.GaussianBlur(1.2)))

    for (p, r, rgb) in atoms:
        shaded_ball(img, p[0], p[1], r, rgb)
    shaded_ball(img, C[0], C[1], 170, (150, 150, 150))  # carbon in the middle, drawn last (in front)

    out = ROOT / "assets"
    out.mkdir(exist_ok=True)
    img.save(out / "icon.png")
    # small copy for the page favicon (also becomes the taskbar icon of the app window)
    img.resize((128, 128), Image.LANCZOS).save(out / "favicon.png", optimize=True)
    img.save(out / "icon.ico", sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    print("wrote", out / "icon.ico")


if __name__ == "__main__":
    main()
