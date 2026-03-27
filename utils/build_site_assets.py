#!/usr/bin/env python3
"""
Build site assets for the Lacemaker Jekyll website.

For each pattern with analyse=1 in jobs.tsv this script:
  - Copies the initial-state JPG  →  assets/plots/{Label}/initial.jpg
  - Copies the GIF timelapse      →  assets/plots/{Label}/{combo}/timelapse.gif
  - Converts lengths PDF → JPG    →  assets/plots/{Label}/{combo}/lengths.jpg
  - Converts hist PDF → JPG       →  assets/plots/{Label}/{combo}/hist.jpg
  - Converts relArea PDF → JPG    →  assets/plots/{Label}/{combo}/relarea.jpg

Run from the repository root:
    python utils/build_site_assets.py

Requirements: Ghostscript (gs) must be on PATH.
"""

import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────
SIM_DIR  = Path("output/simulations")
ASSETS   = Path("assets/plots")
JOBS_TSV = Path("assets/data/jobs.tsv")
# Parameter string embedded in simulation filenames (from lace_maker.py call)
PARAMS   = "0.25_1.0_0.0_30.0_0.1"
DPI      = 300   # resolution for PDF → JPG
JPEGQ    = 85    # JPEG quality

# ── Pattern map: book_id → (label, full_name, family) ────────────────────────
PATTERN_MAP = {
    "1084":    ("S2",  "Square-2",  "Square"),
    "1092":    ("S1",  "Square-1",  "Square"),
    "2001":    ("D1",  "Diamond-1", "Diamond"),
    "2002":    ("D2",  "Diamond-2", "Diamond"),
    "2009":    ("D3",  "Diamond-3", "Diamond"),
    "7002":    ("D4",  "Diamond-4", "Diamond"),
    "7001":    ("D5",  "Diamond-5", "Diamond"),
    "7000":    ("D6",  "Diamond-6", "Diamond"),
    "2023":    ("D7",  "Diamond-7", "Diamond"),
    "2048":    ("H1",  "Hexa-1",    "Hexa"),
    "3023v2":  ("R1",  "Rose-1",    "Rose"),
    "3025":    ("R2",  "Rose-2",    "Rose"),
    "3021":    ("R3",  "Rose-3",    "Rose"),
    "3024":    ("R4",  "Rose-4",    "Rose"),
    "3026":    ("R5",  "Rose-5",    "Rose"),
    "3013":    ("R6",  "Rose-6",    "Rose"),
    "3018":    ("R7",  "Rose-7",    "Rose"),
    "3001":    ("R8",  "Rose-8",    "Rose"),
    "3003":    ("R9",  "Rose-9",    "Rose"),
    "3034":    ("R10", "Rose-10",   "Rose"),
    "3043":    ("R11", "Rose-11",   "Rose"),
    "3051_1t": ("R12", "Rose-12",   "Rose"),
    "3051v2":  ("R12", "Rose-12",   "Rose"),
    "3052":    ("R13", "Rose-13",   "Rose"),
    "3059":    ("Sp1", "Spruce-1",  "Spruce"),
}

# ── Helpers ───────────────────────────────────────────────────────────────────

def pdf_to_jpg(pdf_path: Path, jpg_path: Path) -> bool:
    """Convert first page of a PDF to JPG using Ghostscript. Returns True on success."""
    if not pdf_path.exists():
        return False
    jpg_path.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run([
        "gs", "-dNOPAUSE", "-dBATCH", "-sDEVICE=jpeg",
        f"-dJPEGQ={JPEGQ}", f"-r{DPI}",
        "-dFirstPage=1", "-dLastPage=1",
        f"-sOutputFile={jpg_path}", str(pdf_path)
    ], capture_output=True)
    return result.returncode == 0 and jpg_path.exists()


def copy_if_exists(src: Path, dst: Path) -> bool:
    """Copy src to dst if src exists. Returns True on success."""
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def combo_str(combo: list) -> str:
    return "_".join(str(y) for y in combo)


# ── Core processing ───────────────────────────────────────────────────────────

def process_pattern(pid: str, combos: list, label: str, initial_done: set) -> list:
    """
    Copy/convert all assets for one book ID.
    Returns a list of experiment dicts for use in the .md front matter.
    """
    plots_dir = SIM_DIR / f"Pattern_{pid}" / "plots"
    if not plots_dir.exists():
        print(f"  [WARN] {plots_dir} not found — skipping {pid}")
        return []

    prefix   = f"pattern{pid}_{PARAMS}"
    out_dir  = ASSETS / label
    out_dir.mkdir(parents=True, exist_ok=True)

    # Initial-state JPG (from stabilisation, no contract_fix)
    if label not in initial_done:
        copied = copy_if_exists(
            plots_dir / f"{prefix}_area_dist.jpg",
            out_dir / "initial.jpg"
        )
        if copied:
            print(f"    initial.jpg ✓")
        initial_done.add(label)

    results = []
    for combo in combos:
        cs  = combo_str(combo)
        cp  = f"{prefix}_contract_fix_{cs}"
        # Use variant-prefixed path when the label covers multiple book IDs (e.g. R12)
        # The caller sets is_variant=True in that case via the returned dicts.
        dst = out_dir / cs
        dst.mkdir(exist_ok=True)

        # Timelapse GIF
        has_gif = copy_if_exists(plots_dir / f"{cp}_area_dist.gif",  dst / "timelapse.gif")

        # Yarn lengths PDF → JPG
        has_lengths = pdf_to_jpg(plots_dir / f"{cp}_lengths_subplots.pdf", dst / "lengths.jpg")

        # Pore histogram PDF → JPG  (non-b first, then b-variant fallback)
        has_hist = pdf_to_jpg(plots_dir / f"{cp}_hist.pdf",  dst / "hist.jpg")
        if not has_hist:
            has_hist = pdf_to_jpg(plots_dir / f"{cp}b_hist.pdf", dst / "hist.jpg")

        # Relative area PDF → JPG  (non-b first, then b-variant fallback)
        has_relarea = pdf_to_jpg(plots_dir / f"{cp}_relArea.pdf",  dst / "relarea.jpg")
        if not has_relarea:
            has_relarea = pdf_to_jpg(plots_dir / f"{cp}b_relArea.pdf", dst / "relarea.jpg")

        status = " ".join(filter(None, [
            "gif" if has_gif else "",
            "lengths" if has_lengths else "",
            "hist" if has_hist else "",
            "relarea" if has_relarea else "",
        ])) or "none"
        print(f"    {cs}: [{status}]")

        results.append({
            "combo":       cs,
            "has_gif":     has_gif,
            "has_lengths": has_lengths,
            "has_hist":    has_hist,
            "has_relarea": has_relarea,
        })

    return results


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Reading jobs.tsv ...")
    rows = []
    with open(JOBS_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if str(r["analyse"]).strip() == "1":
                rows.append(r)

    # Group by label so R12's two book IDs are combined
    label_groups = {}
    for r in rows:
        pid = str(r["pattern_id"]).strip()
        if pid not in PATTERN_MAP:
            continue
        label, full_name, family = PATTERN_MAP[pid]
        exp_raw = str(r["exp"]).strip().strip('"')
        try:
            combos = json.loads(exp_raw)
        except Exception:
            combos = []
        if label not in label_groups:
            label_groups[label] = dict(full_name=full_name, family=family, pids=[])
        label_groups[label]["pids"].append((pid, combos))

    initial_done = set()
    summary = {}

    for label, info in label_groups.items():
        print(f"\n[{label}] {info['full_name']}")
        all_exps = []
        for pid, combos in info["pids"]:
            exps = process_pattern(pid, combos, label, initial_done)
            # When a label spans multiple book IDs (e.g. R12), move already-copied
            # assets into a variant subdirectory so paths don't collide.
            if len(info["pids"]) > 1:
                for e in exps:
                    e["variant"] = pid
                    src_dir = ASSETS / label / e["combo"]
                    var_dir = ASSETS / label / pid / e["combo"]
                    if src_dir.exists() and not var_dir.exists():
                        var_dir.parent.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(src_dir), str(var_dir))
            all_exps.extend(exps)
        summary[label] = all_exps

    print("\n\nSummary")
    print("-------")
    for label, exps in summary.items():
        n = len(exps)
        gifs  = sum(1 for e in exps if e["has_gif"])
        jpgs  = sum(1 for e in exps if e["has_lengths"] or e["has_hist"] or e["has_relarea"])
        print(f"  {label:4s}  {n} experiment(s)  {gifs} GIFs  {jpgs} plots converted")

    print("\nDone. Commit assets/plots/ to publish the images.")


if __name__ == "__main__":
    main()
