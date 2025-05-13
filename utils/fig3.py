import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import os
import matplotlib.ticker as ticker
from natsort import index_natsorted

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

# load
df = pd.read_csv(r"/scratch/work/silvap1/lacemaker/all_results.csv")

# sort & color‐map
df = df.iloc[index_natsorted(df[sample_name])]
unique_samples = df[sample_name].unique()
palette = matplotlib.colormaps['tab20']
color_map = {s: palette(i % 20) for i, s in enumerate(unique_samples)}

# mark highest‐contract point per sample
df['contract_count'] = df['contract'].str.count('_') + 1
idx_max_contract = df.groupby(sample_name)['contract_count'].idxmax()
df['highlight'] = df.index.isin(idx_max_contract)

# collect handles/labels for legend
handles, labels = [], []

# --- 1) create the main scatter plot (no legend) ---
fig, ax = plt.subplots(figsize=(5.0/2.54, 6.0/2.54))
for sample in unique_samples:
    sub = df[df[sample_name] == sample]
    mark = sub['marker'].iloc[0]
    col  = f"#{sub['color'].iloc[0]}" #color_map[sample]
    h = ax.scatter(
        -sub['rel_max_contraction'],
        sub['area_change_t10'],
        marker=mark,
        facecolor=col,
        edgecolor='black',
        linewidth=0.5,
        s=10,
        alpha=0.8,
        zorder=3
    )
    # only add one handle per sample
    handles.append(h)
    labels.append(sample)

# formatting
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

# save data‐only plot
fig.savefig('fig.pdf', format='pdf', dpi=300, bbox_inches='tight')
convert_fonts_to_outlines("fig.pdf", "fig_c.pdf")
os.replace("fig_c.pdf", "fig.pdf")
plt.close(fig)


# --- 2) create & save the legend by itself ---
# build ordered mapping to remove duplicates
by_label = dict(zip(labels, handles))

# figure size tuned to legend contents
fig_leg = plt.figure(figsize=(2.5, 3.0))
ax_leg = fig_leg.add_subplot(111)
ax_leg.axis('off')  # no axes

leg = ax_leg.legend(
    by_label.values(),
    by_label.keys(),
    loc='center',
    frameon=False,
    ncol=1,
    fontsize=5
)
fig_leg.tight_layout()

fig_leg.savefig('legend.pdf', format='pdf', dpi=300, bbox_inches='tight')
convert_fonts_to_outlines("legend.pdf", "legend_c.pdf")
os.replace("legend_c.pdf", "legend.pdf")
plt.close(fig_leg)