"""
Convert whole/ PNG snapshots → assets/images/{book_id}_pattern.jpg

For each pattern, loads the raw Ovito PNG, resizes to TARGET_W px wide,
appends a yarn-colour legend strip at the bottom, and saves as JPEG.

Run from repository root:
    python utils/build_mesh_images.py
"""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# ── Config ────────────────────────────────────────────────────────────────────
WHOLE_DIR = Path("whole")
OUT_DIR   = Path("assets/images")
TARGET_W  = 600      # output width in px
JPEG_Q    = 90

# Ovito default colours blended 60 % original + 40 % white → pastel
YARN_COLORS = {
    1: (255, 163, 163),   # soft red
    2: (163, 163, 255),   # soft blue
    3: (255, 255, 102),   # soft yellow
    4: (255, 163, 255),   # soft magenta
    5: (163, 255, 133),   # soft green
    6: (224, 255, 209),   # soft mint
    7: (209, 102, 255),   # soft purple
    8: (133, 255, 255),   # soft cyan
}

# book_id → filename stem in whole/   (no .png)
ID_TO_FILE = {
    "1084":     "1084",
    "1092":     "1092",
    "2001":     "2001",
    "2002":     "2002",
    "2009":     "2009",
    "7002":     "7002",
    "7001":     "7001",
    "7000":     "7000",
    "2023":     "2023",
    "2048":     "2048",
    "3023v2":   "3023",
    "3025":     "3025",
    "3021":     "3021",
    "3024":     "3024",
    "3026":     "3026",
    "3013":     "3013",
    "3018":     "3018",
    "3001":     "3001",
    "3003":     "3003",
    "3034":     "3034",
    "3043":     "3043",
    "3051_1t":  "3051",
    "3051v2":   "3051-03",
    "3052":     "3052",
    "3059":     "3059",
}

# ── Legend ────────────────────────────────────────────────────────────────────
LEGEND_H    = 52    # height of the legend strip in px
SWATCH_W    = 28    # coloured square width
SWATCH_H    = 24    # coloured square height
PAD         = 8     # padding between elements
FONT_SIZE   = 16

def make_legend(width: int) -> Image.Image:
    """Return a white PIL image containing the 8-yarn colour legend."""
    img  = Image.new("RGB", (width, LEGEND_H), (255, 255, 255))
    draw = ImageDraw.Draw(img)

    try:
        font = ImageFont.truetype("/usr/share/fonts/liberation/LiberationSans-Regular.ttf", FONT_SIZE)
    except OSError:
        font = ImageFont.load_default()

    n = len(YARN_COLORS)
    # total width needed per item: swatch + 2px gap + label + PAD
    label_w   = 14   # approximate char width for single digit
    item_w    = SWATCH_W + 4 + label_w + PAD
    block_w   = item_w * n - PAD
    x0        = (width - block_w) // 2
    y_swatch  = (LEGEND_H - SWATCH_H) // 2

    for i, (yarn, color) in enumerate(YARN_COLORS.items()):
        x = x0 + i * item_w
        # coloured square with thin border
        draw.rectangle([x, y_swatch, x + SWATCH_W, y_swatch + SWATCH_H],
                       fill=color, outline=(160, 160, 160))
        # yarn number
        tx = x + SWATCH_W + 4
        ty = y_swatch + (SWATCH_H - FONT_SIZE) // 2
        draw.text((tx, ty), str(yarn), fill=(50, 50, 50), font=font)

    return img


def convert(book_id: str, stem: str) -> bool:
    src = WHOLE_DIR / f"{stem}.png"
    if not src.exists():
        print(f"  WARN: {src} not found — skipping {book_id}")
        return False

    img = Image.open(src).convert("RGB")

    # Resize to TARGET_W, keep aspect ratio
    w, h = img.size
    new_h = round(h * TARGET_W / w)
    img   = img.resize((TARGET_W, new_h), Image.LANCZOS)

    # Append legend strip
    legend = make_legend(TARGET_W)
    combined = Image.new("RGB", (TARGET_W, new_h + LEGEND_H), (255, 255, 255))
    combined.paste(img,    (0, 0))
    combined.paste(legend, (0, new_h))

    out = OUT_DIR / f"{book_id}_pattern.jpg"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    combined.save(out, "JPEG", quality=JPEG_Q, optimize=True)
    print(f"  {book_id:10s} → {out}  ({combined.size[0]}×{combined.size[1]})")
    return True


if __name__ == "__main__":
    print(f"Converting {len(ID_TO_FILE)} patterns …")
    ok = sum(convert(bid, stem) for bid, stem in ID_TO_FILE.items())
    print(f"Done: {ok}/{len(ID_TO_FILE)} written to {OUT_DIR}/")
