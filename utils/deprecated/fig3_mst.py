import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import os
import matplotlib.ticker as ticker
from natsort import index_natsorted
from scipy.spatial import distance_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

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

# natural sort and color‐map
df = df.iloc[index_natsorted(df[sample_name])]
unique_samples = df[sample_name].unique()
palette = matplotlib.colormaps['tab20']
color_map = {s: palette(i % 20) for i, s in enumerate(unique_samples)}

# mark highest‐contract point per sample (optional)
df['contract_count'] = df['contract'].str.count('_') + 1
idx_max_contract = df.groupby(sample_name)['contract_count'].idxmax()
df['highlight'] = df.index.isin(idx_max_contract)

# --- 1) scatter + MST overlay ---
fig, ax = plt.subplots(figsize=(5.5/2.54, 6.25/2.54))
for sample in unique_samples:
    sub = df[df[sample_name] == sample]
    # extract coordinates
    x = -sub['rel_max_contraction'].to_numpy()
    y = sub['area_change_t10'].to_numpy()
    col = f"#{sub['color'].iloc[0]}"
    marker = sub['marker'].iloc[0]

    # separate highlighted point
    highlight_mask = sub['highlight'].to_numpy()

    # normal points
    # ax.scatter(
    #     x[~highlight_mask], y[~highlight_mask],
    #     marker=marker,
    #     facecolor=col,
    #     edgecolor='black',
    #     linewidth=0.25,
    #     s=4,
    #     zorder=3
    # )

    # highlighted point (larger)
    ax.scatter(
        x[highlight_mask], y[highlight_mask],
        marker=marker,
        facecolor=col,
        edgecolor='black',
        linewidth=0.5,
        s=10,  # increase size here
        zorder=4
    )

    # compute MST for this sample's points
    if len(x) > 1:
        pts = np.column_stack((x, y))
        D = distance_matrix(pts, pts)
        Tcsr = minimum_spanning_tree(D)
        rows, cols = Tcsr.nonzero()
        # draw edges colored by sample
        for i, j in zip(rows, cols):
            xi, yi = x[i], y[i]
            xj, yj = x[j], y[j]
            line, = ax.plot([xi, xj], [yi, yj], color=col, linewidth=1.0, alpha=1.0, zorder=2)
            line.set_solid_capstyle('round')

# formatting axes
ax.set_xlabel(r'Rel. max contraction')
ax.set_ylabel(r'$\delta A_{10}$')
ax.set_xlim(-0.0125, 0.275)
ax.set_ylim(-0.025, 0.8)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.hlines(0.0, -1, 1, color='black', linewidth=0.5, linestyles='--', zorder=-1)
ax.vlines(0.0, -1, 1, color='black', linewidth=0.5, linestyles='--', zorder=-1)
fig.tight_layout()

# save main figure
fig.savefig('fig.pdf', format='pdf', dpi=300, bbox_inches='tight')
convert_fonts_to_outlines('fig.pdf', 'fig_c.pdf')
os.replace('fig_c.pdf', 'fig.pdf')
plt.close(fig)

# --- 2) legend-only ---
by_label = {sample: plt.Line2D([0], [0], marker=df[df[sample_name]==sample]['marker'].iloc[0],
                              color='w', markerfacecolor=f"#{df[df[sample_name]==sample]['color'].iloc[0]}",
                              markeredgecolor='k', markersize=4, markeredgewidth=0.25)
            for sample in unique_samples}

fig_leg = plt.figure(figsize=(2.5, 3.0))
ax_leg = fig_leg.add_subplot(111)
ax_leg.axis('off')
leg = ax_leg.legend(
    by_label.values(), by_label.keys(),
    loc='center', frameon=False, ncol=1, fontsize=5, 
)
fig_leg.tight_layout()
fig_leg.savefig('legend.pdf', format='pdf', dpi=300, bbox_inches='tight')
convert_fonts_to_outlines('legend.pdf', 'legend_c.pdf')
os.replace('legend_c.pdf', 'legend.pdf')
plt.close(fig_leg)
