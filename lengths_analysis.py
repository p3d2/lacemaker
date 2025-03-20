# lengths_analysis.py
import os
import sys
import json
import subprocess
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 9})
FILTER_L0 = 48.0

color_map = [
    "#FF5C5C", "#5C5CFF", "#FFD966", "#FF5CFF",
    "#5CFF5C", "#A6FFCC", "#995CFF", "#5CFFFF"
]

def convert_fonts_to_outlines(input_pdf, output_pdf):
    gs_command = [
        "gs", "-o", output_pdf, "-sDEVICE=pdfwrite",
        "-dNoOutputFonts",  # convert fonts to outlines
        "-dNOPAUSE", "-dBATCH",
        "-dSubsetFonts=true", "-dEmbedAllFonts=true",
        "-dDownsampleColorImages=false",
        "-dDownsampleGrayImages=false",
        "-dDownsampleMonoImages=false",
        input_pdf
    ]
    subprocess.run(gs_command, check=True)

def plot_lengths(json_file, pattern_folder, dump_file):
    # 1) Read JSON
    with open(json_file, 'r') as f:
        data = json.load(f)  # list of { time, yarns: [ {mol, type, length}, ... ] }

    # 2) Build a dictionary keyed by (mol, type), storing times and lengths
    yarn_dict = {}  # (mol, type) -> { 'times': [], 'lengths': [] }

    for entry in data:
        t = entry['time']
        for y in entry['yarns']:
            mol_val = y['mol']
            typ_val = y['type']
            length  = y['length']
            key = (mol_val, typ_val)
            if key not in yarn_dict:
                yarn_dict[key] = {
                    'times': [],
                    'lengths': []
                }
            yarn_dict[key]['times'].append(t)
            yarn_dict[key]['lengths'].append(length)

    # 3) For each (mol, type) series, sort by time and compute (L(t)-L0)/L0
    type_groups = {}  # type -> list of dicts { 'times': [...], 'relvar': [...], 'mol': mol_val }
    for (mol_val, typ_val), info in yarn_dict.items():
        times_array = np.array(info['times'])
        lengths_array = np.array(info['lengths'])
        sort_idx = np.argsort(times_array)
        times_sorted = times_array[sort_idx]
        lengths_sorted = lengths_array[sort_idx]

        # Filter by baseline length
        L0 = lengths_sorted[0]
        if L0 < FILTER_L0:
            # Skip this yarn entirely
            continue

        # Compute rel. variation
        rel_var = (lengths_sorted - L0) / L0
        
        if typ_val not in type_groups:
            type_groups[typ_val] = []
        type_groups[typ_val].append({
            'times': times_sorted,
            'relvar': rel_var,
            'mol': mol_val
        })

    # 4) Create subplots: 2 rows x 4 columns, one subplot per type
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(5.0 / 2.54, 5.04 / 2.54),
                             sharex=True, sharey=True)
    axes = axes.flatten()

    all_types = [1,2,3,4,5,6,7,8]
    (ymin, ymax) = (-0.05, 0.05)
    for i, ttype in enumerate(all_types):
        ax = axes[i]
        if ttype not in type_groups:
            ax.text(0.5, 0.98, f"Type {ttype}", fontsize=7, ha="center", va="top", transform=ax.transAxes)
            continue

        group_list = type_groups[ttype]
        
        for yarn_entry in group_list:
            times = yarn_entry['times']

            nframes = 20
            contract = 0.25
            formula1_values = [1 - (t / (nframes/2)) * contract for t in range(1, int(nframes/2) + 1)]
            formula2_values = [1 - contract + ((t - (nframes/2)) / (nframes/2)) * contract for t in range(int(nframes/2) + 1, nframes + 1)]
            delta = [1] + formula1_values + formula2_values

            relv  = yarn_entry['relvar']
            c = color_map[(ttype - 1) % 8]
            ax.plot(times, relv, alpha=0.5, color=c,  zorder=10)
            if np.min(relv) < ymin: ymin = np.min(relv)
            if np.max(relv) > ymax: ymax = np.max(relv)        

        ax.text(0.5, 0.98, f"Type {ttype}", fontsize=7, ha="center", va="top", transform=ax.transAxes)
        ax.axhline(0.0, color='gray', lw=1, linestyle='--')     

    fig.supxlabel(r'Simulation Step', fontsize=9, y=-0.075)

    for ax in axes[0::4]:
        ax.set_ylabel(r'$\frac{L - L_0}{L_0}$')

    for ax in axes:
        ax.set_xlim(0, 20)
        if ymin < -0.05: ymin = math.floor(ymin / 0.05) * 0.05 
        if ymax > 0.05: ymax = math.ceil(ymax / 0.05) * 0.05
        
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05)) # Major ticks every 0.5 units
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.025))
        if ymin < -0.05 or ymax > 0.05:
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

        ax.set_ylim(ymin, ymax)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))   # Major ticks every 5 units
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))   # Minor ticks every 1 unit
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '' if x == 0 else str(int(x))))

    pdf_file = os.path.join(pattern_folder, 'plots', dump_file + "_lengths_subplots.pdf")
    fig.savefig(pdf_file, dpi=300, bbox_inches='tight')
    compressed_pdf = pdf_file.replace(".pdf", "_c.pdf")
    convert_fonts_to_outlines(pdf_file, compressed_pdf)
    os.replace(compressed_pdf, pdf_file)

def main(pattern, base_dir="output/simulations"):
    # 1) Find subdirectories named "Pattern_3025*"
    pattern_dirs = [
        d for d in os.listdir(base_dir)
        if d.startswith(pattern) and os.path.isdir(os.path.join(base_dir, d))
    ]

    # 2) For each pattern directory, look for the matching *_yarn_lengths.json
    for pattern_dir in pattern_dirs:
        folder_path = os.path.join(base_dir, pattern_dir)
        plots_folder = os.path.join(folder_path, "plots")
        if not os.path.exists(plots_folder):
            continue

        # 3) Gather all JSON files that match our naming convention.
        for file_name in os.listdir(plots_folder):
            if file_name.endswith("_yarn_lengths.json"):
                json_file = os.path.join(plots_folder, file_name)
                dump_file = file_name.replace("_yarn_lengths.json", "")
                # 4) Call plot_lengths with the discovered JSON.
                plot_lengths(json_file, folder_path, dump_file)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Error: Missing required argument.")
        print("Usage: python script.py <pattern_name>")
        sys.exit(1)  # Exit with error

    pattern_arg = sys.argv[1]
    main(pattern=pattern_arg)