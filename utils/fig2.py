import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import os

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 9})


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

# — Manual data entry: fill in your own values below —
data = [
    {'pattern': 'K', 'sim_id': 1, 'performance': -0.09, 'tag': 1, 'marker': '+'},
    {'pattern': 'W', 'sim_id': 1, 'performance': 0.04, 'tag': 1, 'marker': 'x'},
    {'pattern': 'D1', 'sim_id': 1, 'performance': -0.03,  'tag': 1, 'marker': 'D'},
    {'pattern': 'D2', 'sim_id': 1, 'performance': -0.04,  'tag': 1, 'marker': 'D'},
    {'pattern': 'D3', 'sim_id': 1, 'performance': 0.19,  'tag': 2, 'marker': 'D'},
    {'pattern': 'D4', 'sim_id': 1, 'performance': 0.11,  'tag': 2, 'marker': 'D'},
    {'pattern': 'D5', 'sim_id': 1, 'performance': 0.59, 'tag': 2, 'marker': 'D'},
    {'pattern': 'D6', 'sim_id': 1, 'performance': 0.02, 'tag': 2, 'marker': 'D'},
    {'pattern': 'D7', 'sim_id': 1, 'performance': -0.01, 'tag': 1, 'marker': 'D'},
    {'pattern': 'S1', 'sim_id': 1, 'performance': 0.09, 'tag': 1, 'marker': 's'},
    {'pattern': 'S2', 'sim_id': 1, 'performance': -0.04, 'tag': 1, 'marker': 's'},
    {'pattern': 'H1', 'sim_id': 1, 'performance': 0.62, 'tag': 1, 'marker': 'h'},
    {'pattern': 'R1', 'sim_id': 1, 'performance': 0.97, 'tag': 1, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 2, 'performance': 0.75, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 3, 'performance': 0.95, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 4, 'performance': 0.90, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 5, 'performance': 0.95, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 6, 'performance': 0.89, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 7, 'performance': 0.62, 'tag': 4, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 8, 'performance': 0.53, 'tag': 4, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 9, 'performance': 0.79, 'tag': 4, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 10, 'performance': 0.63, 'tag': 4, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 11, 'performance': 0.90, 'tag': 4, 'marker': 'o'},
    {'pattern': 'R1', 'sim_id': 12, 'performance': 0.35, 'tag': 8, 'marker': 'o'},
    {'pattern': 'R2', 'sim_id': 1, 'performance': 0.44, 'tag': 1, 'marker': 'o'},
    {'pattern': 'R3', 'sim_id': 1, 'performance': 0.30, 'tag': 1, 'marker': 'o'},
    {'pattern': 'R4', 'sim_id': 1, 'performance': 0.27, 'tag': 1, 'marker': 'o'},
    {'pattern': 'R4', 'sim_id': 2, 'performance': 0.22, 'tag': 8, 'marker': 'o'},
    {'pattern': 'R4', 'sim_id': 3, 'performance': 0.31, 'tag': 4, 'marker': 'o'},
    {'pattern': 'R4', 'sim_id': 4, 'performance': 0.45, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R4', 'sim_id': 5, 'performance': 0.43, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R5', 'sim_id': 1, 'performance': 0.31, 'tag': 1, 'marker': 'o'},
    {'pattern': 'R6', 'sim_id': 1, 'performance': 0.68, 'tag': 1, 'marker': 'o'},
    {'pattern': 'R7', 'sim_id': 1, 'performance': 0.63, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R8', 'sim_id': 1, 'performance': 0.13, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R8', 'sim_id': 2, 'performance': 0.13, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R9', 'sim_id': 1, 'performance': 0.10, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R9', 'sim_id': 2, 'performance': 0.07, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R9', 'sim_id': 3, 'performance': 0.07, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R9', 'sim_id': 4, 'performance': 0.13, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R10', 'sim_id': 1, 'performance': 0.43, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R11', 'sim_id': 1, 'performance': 0.13, 'tag': 2, 'marker': 'o'},
    {'pattern': 'R12', 'sim_id': 1, 'performance': 0.35, 'tag': 1, 'marker': 'o'},
    {'pattern': 'R12', 'sim_id': 2, 'performance': 0.26, 'tag': 4, 'marker': 'o'},
    {'pattern': 'R12', 'sim_id': 3, 'performance': -0.16, 'tag': 4, 'marker': 'o'},
    {'pattern': 'R13', 'sim_id': 1, 'performance': 0.32, 'tag': 1, 'marker': 'o'},
    {'pattern': 'Sp1', 'sim_id': 1, 'performance': 0.37, 'tag': 4, 'marker': '^'},
]

# Create DataFrame
df = pd.DataFrame(data)

# 1) performance range per pattern
perf_range = df.groupby('pattern')['performance'] \
               .agg(['min', 'max']) \
               .to_dict(orient='index')

# 2) count of sims per pattern
pattern_counts = df['pattern'].value_counts().to_dict()

# — Plotting —
fig, ax = plt.subplots(figsize=(9.5 / 2.54, 6.0 / 2.54))

# Determine category positions
categories = df['pattern'].unique()
x_locs = np.arange(len(categories))
pos_map = {cat: x for x, cat in enumerate(categories)}

# Build a color‐map from your unique patterns
unique_patterns = df['pattern'].unique()
palette = matplotlib.colormaps['tab20']
color_map = {p: palette(i % 20) for i, p in enumerate(unique_patterns)}

# Compute how many share the exact (pattern,performance)
df['dup_count'] = (
    df.groupby(['pattern','performance'])['performance']
      .transform('count')
)

xw = df['pattern'].nunique()
jitter_amount = 0.5

seen_patterns = set()
np.random.seed(0)
# Plot each simulation
for _, row in df.iterrows():
    pattern = row['pattern']
    x0 = pos_map[row['pattern']]
    y  = row['performance']
    c  = color_map[row['pattern']]
    if row['dup_count'] > 1:
        x = x0 + np.random.uniform(-jitter_amount, jitter_amount)
    else:
        x = x0

    # draw one full-range stem per pattern
    if pattern not in seen_patterns:
        # if only one point, start stem at 0; otherwise at the min performance
        if pattern_counts[pattern] == 1:
            y_min = 0.0
        else:
            if perf_range[pattern]['min'] < 0.0:
                y_min = perf_range[pattern]['min']
            else:
                y_min = 0.0

        y_max = perf_range[pattern]['max']

        ax.vlines(x0, y_min, y_max, color=color_map[row['pattern']], linewidth=0.5, zorder=-1)
        seen_patterns.add(pattern)

    y = row['performance']
    ax.scatter(x, y, color=color_map[row['pattern']],
               s=18, zorder=3, marker=row['marker'])
    ax.text(x+0.3, y+0.01, str(row['tag']),
            ha='center', va='bottom', fontsize=5, zorder=4)

# Formatting the axes
ax.set_xticks(x_locs)
ax.set_xticklabels(categories, fontsize=7, rotation=60, ha='center')
ax.set_ylabel(r'$\delta_\mathrm{{end}}$')
ax.grid(axis='y', linestyle='--', alpha=0.5, linewidth=0.5)
ax.set_xlim(-0.5, xw-0.5)
ax.set_ylim(-0.2, 1.0)
fig.savefig('fig.pdf', format='pdf', dpi=300, bbox_inches='tight')

convert_fonts_to_outlines("fig.pdf", "fig_c.pdf")
os.replace("fig_c.pdf", "fig.pdf")