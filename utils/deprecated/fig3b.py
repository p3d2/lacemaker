import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import os
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from scipy.spatial import ConvexHull
from scipy.interpolate import splprep, splev
from natsort import index_natsorted

def smooth_hull(pts, expand=1.05, n_pts=200, smoothness=0.01):
    hull = ConvexHull(pts)
    hp = pts[hull.vertices]
    hp = np.vstack([hp, hp[0]])  # close loop
    tck, u = splprep([hp[:,0], hp[:,1]], s=smoothness, per=True)
    u_fine = np.linspace(0, 1, n_pts)
    x_s, y_s = splev(u_fine, tck)
    cen = np.column_stack((x_s, y_s)).mean(axis=0)
    xp = cen[0] + (x_s - cen[0]) * expand
    yp = cen[1] + (y_s - cen[1]) * expand
    return xp, yp

# LaTeX labels, font size
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 9})
sample_name = 'sample_reduced'

def convert_fonts_to_outlines(input_pdf, output_pdf):
    gs_command = [
        "gs", "-o", output_pdf, "-sDEVICE=pdfwrite",
        "-dNoOutputFonts",
        "-dNOPAUSE", "-dBATCH",
        "-dSubsetFonts=true", "-dEmbedAllFonts=true",
        "-dDownsampleColorImages=false",
        "-dDownsampleGrayImages=false",
        "-dDownsampleMonoImages=false",
        input_pdf
    ]
    subprocess.run(gs_command, check=True)

# load data
df = pd.read_csv(r"/scratch/work/silvap1/lacemaker/all_results.csv")

# --- 1) sort by sample_full so legend/order follows alphabetically ---
df = df.iloc[index_natsorted(df[sample_name])]

# --- 2) build a color‐map per sample_full ---
unique_samples = df[sample_name].unique()
palette = matplotlib.colormaps['tab20']
color_map = {s: palette(i % 20) for i, s in enumerate(unique_samples)}

# count pieces in contract string & find highlight rows
df['contract_count'] = df['contract'].str.count('_') + 1
idx_max_contract = df.groupby(sample_name)['contract_count'].idxmax()
highlight_idx = set(idx_max_contract)

# --- 3) make the scatter plot, collecting handles for legend ---
fig, ax = plt.subplots(figsize=(5.0/2.54, 6.0/2.54))
handles, labels = [], []

for sample in unique_samples:
    sub = df[df[sample_name] == sample]
    mark = sub['marker'].iloc[0]
    col  = f"#{sub['color'].iloc[0]}" #color_map[sample]
    # plot all points for this sample at once
    h = ax.scatter(
        sub['area_change_t20'],
        sub['area_change_t10'] - sub['area_change_t20'],
        marker=mark,
        facecolor=col,
        edgecolor='black', #['red' if i in highlight_idx else 'black' for i in sub.index],
        linewidth=0.5,
        s=10,
        alpha=0.8,
        zorder=3
    )
    handles.append(h)
    labels.append(sample)

# --- 4) formatting ---
ax.set_xlabel(r'$\delta A_{20}$')
ax.set_ylabel(r'$\delta A_{10} - \delta A_{20}$')

ax.hlines([0.0], -1, 1, color='black', linewidth=0.5,
          linestyles='--', zorder=-1)
ax.vlines([0.0], -1, 1, color='black', linewidth=0.5,
          linestyles='--', zorder=-1)

# optional rectangles (uncomment to use)
# rects = [
#     patches.Rectangle((0.0, 0.0), 0.1, 0.1, facecolor='red',   alpha=0.2, zorder=0),
#     patches.Rectangle((0.1,0.0), 0.3, 0.1, facecolor='orange',alpha=0.2, zorder=0),
#     patches.Rectangle((0.0,0.1), 0.1, 0.3, facecolor='green', alpha=0.2, zorder=0),
# ]
# for r in rects: ax.add_patch(r)

ax.set_xlim(-0.025, 0.525)
ax.set_ylim(-0.025, 0.525)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

# --- legend outside (same as first code) ---
by_label = dict(zip(labels, handles))
# ax.legend(
#     by_label.values(),
#     by_label.keys(),
#     bbox_to_anchor=(1.02, 1),
#     loc='upper left',
#     fontsize=5,
#     frameon=False,
#     ncol=2
# )

fig.tight_layout()

# save only this combined figure
out_pdf = 'figb.pdf'
fig.savefig(out_pdf, format='pdf', dpi=300, bbox_inches='tight')
convert_fonts_to_outlines(out_pdf, 'figb_c.pdf')
os.replace('figb_c.pdf', out_pdf)

plt.close(fig)
