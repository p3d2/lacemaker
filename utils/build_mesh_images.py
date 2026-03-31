"""
Convert whole/ PNG snapshots → assets/images/{book_id}_pattern.jpg

For each pattern, loads the raw Ovito PNG (RGBA, transparent background),
composites onto white, resizes to TARGET_W px wide, and saves as JPEG.

Run from repository root:
    python utils/build_mesh_images.py
"""

from pathlib import Path
from PIL import Image

# ── Config ────────────────────────────────────────────────────────────────────
WHOLE_DIR = Path("whole")
OUT_DIR   = Path("assets/images")
TARGET_W  = 600      # output width in px
JPEG_Q    = 90

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

def convert(book_id: str, stem: str) -> bool:
    src = WHOLE_DIR / f"{stem}.png"
    if not src.exists():
        print(f"  WARN: {src} not found — skipping {book_id}")
        return False

    img = Image.open(src).convert("RGBA")

    # Composite onto white background (replaces transparent → white)
    w, h = img.size
    white = Image.new("RGB", (w, h), (255, 255, 255))
    white.paste(img, mask=img.split()[3])  # alpha channel as mask
    img = white

    # Resize to TARGET_W, keep aspect ratio
    new_h = round(h * TARGET_W / w)
    img   = img.resize((TARGET_W, new_h), Image.LANCZOS)

    out = OUT_DIR / f"{book_id}_pattern.jpg"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img.save(out, "JPEG", quality=JPEG_Q, optimize=True)
    print(f"  {book_id:10s} → {out}  ({img.size[0]}×{img.size[1]})")
    return True


if __name__ == "__main__":
    print(f"Converting {len(ID_TO_FILE)} patterns …")
    ok = sum(convert(bid, stem) for bid, stem in ID_TO_FILE.items())
    print(f"Done: {ok}/{len(ID_TO_FILE)} written to {OUT_DIR}/")
