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

# ----------------------------- helper ----------------------------------
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 9})

sample_name = 'sample_reduced'

def convert_fonts_to_outlines(input_pdf: str, output_pdf: str):
    """Run Ghostscript so every font is converted to outlines (no text).
    Parameters
    ----------
    input_pdf : str
        Path to input PDF.
    output_pdf : str
        Path where the converted PDF will be saved.
    """
    gs_command = [
        "gs", "-o", output_pdf, "-sDEVICE=pdfwrite",
        "-dNoOutputFonts", "-dNOPAUSE", "-dBATCH",
        "-dSubsetFonts=true", "-dEmbedAllFonts=true",
        "-dDownsampleColorImages=false",
        "-dDownsampleGrayImages=false",
        "-dDownsampleMonoImages=false",
        input_pdf,
    ]
    subprocess.run(gs_command, check=True)

# ----------------------------- data ------------------------------------

df = pd.read_csv(r"/scratch/work/silvap1/lacemaker/all_results.csv")

# natural sort so legend order is alphabetical
df = df.iloc[index_natsorted(df[sample_name])]
unique_samples = df[sample_name].unique()

# colour and marker come from the CSV itself
palette = matplotlib.colormaps['tab20']
color_map = {s: palette(i % 20) for i, s in enumerate(unique_samples)}

# -----------------------------------------------------------------------
# find highlight point (row with max contract parts) per sample ----------
# -----------------------------------------------------------------------

df['contract_count'] = df['contract'].str.count('_') + 1
highlight_idx = set(
    df.groupby(sample_name)['contract_count'].idxmax()
)

# ----------------------------- figure ----------------------------------

fig, ax = plt.subplots(figsize=(5.5 / 2.54, 6.25 / 2.54))
handles, labels = [], []

for sample in unique_samples:
    sub = df[df[sample_name] == sample]
    marker = sub['marker'].iloc[0]
    col = f"#{sub['color'].iloc[0]}"

    # coordinates for this pattern
    x = sub['area_change_t20'].to_numpy()
    y = (sub['area_change_t10'] - sub['area_change_t20']).to_numpy()

    # mask for the highlight point
    mask_high = np.isin(sub.index, list(highlight_idx))

    # normal points (small)
    ax.scatter(
        x[~mask_high], y[~mask_high],
        marker=marker,
        facecolor=col,
        edgecolor='black',
        linewidth=0.25,
        s=4,
        zorder=3
    )

    # highlight point (larger)
    ax.scatter(
        x[mask_high], y[mask_high],
        marker=marker,
        facecolor=col,
        edgecolor='black',
        linewidth=0.5,
        s=10,  
        zorder=4
    )

    # minimum‑spanning tree for this pattern (≥2 points)
    if len(x) > 1:
        pts = np.column_stack((x, y))
        D = distance_matrix(pts, pts)
        Tcsr = minimum_spanning_tree(D)
        rows, cols = Tcsr.nonzero()
        for i, j in zip(rows, cols):
            ax.plot(
                [x[i], x[j]], [y[i], y[j]],
                color=col,
                linewidth=1.0,
                alpha=0.8,
                zorder=2,
            )

    # store handle for legend (one per pattern)
    h = ax.scatter([], [], marker=marker, facecolor=col, edgecolor='black', s=20)
    handles.append(h)
    labels.append(sample)

# ----------------------------- axes style ------------------------------

ax.set_xlabel(r"$\delta A_{20}$")
ax.set_ylabel(r"$\delta A_{10} - \delta A_{20}$")

ax.hlines(0.0, -1, 1, color='black', linewidth=0.5, linestyles='--', zorder=-1)
ax.vlines(0.0, -1, 1, color='black', linewidth=0.5, linestyles='--', zorder=-1)

ax.set_xlim(-0.025, 0.525)
ax.set_ylim(-0.025, 0.475)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

fig.tight_layout()

# ----------------------------- save ------------------------------------

out_pdf = 'figb.pdf'
fig.savefig(out_pdf, format='pdf', dpi=300, bbox_inches='tight')
convert_fonts_to_outlines(out_pdf, 'figb_c.pdf')
os.replace('figb_c.pdf', out_pdf)
plt.close(fig)

# ----------------------------- legend ----------------------------------

by_label = dict(zip(labels, handles))
fig_leg = plt.figure(figsize=(2.5, 3.0))
ax_leg = fig_leg.add_subplot(111)
ax_leg.axis('off')
ax_leg.legend(
    by_label.values(), by_label.keys(),
    loc='center', frameon=False, ncol=1, fontsize=5,
)
fig_leg.tight_layout()

fig_leg.savefig('legendb.pdf', format='pdf', dpi=300, bbox_inches='tight')
convert_fonts_to_outlines('legendb.pdf', 'legendb_c.pdf')
os.replace('legendb_c.pdf', 'legendb.pdf')
plt.close(fig_leg)
