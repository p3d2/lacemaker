import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import os
import matplotlib.colors as mcolors
from natsort import index_natsorted
import matplotlib.ticker as ticker

# ----------------------------- settings --------------------------------
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 9})
sample_name = 'sample_reduced'
bar_width = 0.15/2.54          # width of each bar in inches
group_gap = bar_width * 1.5     # gap between sample groups
left_margin = 0.0             # inches
right_margin = 0.0            # inches
height_full = 6.0/2.54        # original full height (in inches)
height_half = height_full / 2 # half height for split plots

# Turn off axis 2
axes2 = False

# Ghostscript outline conversion helper
def convert_fonts_to_outlines(input_pdf, output_pdf):
    gs = [
        'gs', '-o', output_pdf, '-sDEVICE=pdfwrite',
        '-dNoOutputFonts', '-dNOPAUSE', '-dBATCH',
        '-dSubsetFonts=true', '-dEmbedAllFonts=true',
        '-dDownsampleColorImages=false',
        '-dDownsampleGrayImages=false',
        '-dDownsampleMonoImages=false',
        input_pdf
    ]
    subprocess.run(gs, check=True)

# ---------------------------- load & sort ------------------------------
df = pd.read_csv(r"/scratch/work/silvap1/lacemaker/all_results.csv")
df = df.iloc[index_natsorted(df[sample_name])]
unique_samples = list(df[sample_name].unique())
df['contract_count'] = df['contract'].str.count('_') + 1
# find split index (inclusive 'R5')
split_idx = unique_samples.index('R5') + 1
# two sample lists
groups = [unique_samples[:split_idx], unique_samples[split_idx:]]

for gi, samples in enumerate(groups, start=1):
    # count bars per sample and compute x positions
    bar_counts = [len(df[df[sample_name] == s]) for s in samples]
    x_positions = []
    xticks = []
    xticklabels = []
    curr_x = 0.0
    for s, count in zip(samples, bar_counts):
        xs = curr_x + np.arange(count) * bar_width
        x_positions.append((s, xs))
        xticks.append(xs.mean())
        xticklabels.append(s)
        curr_x += count * bar_width + group_gap
    total_width = curr_x - group_gap
    fig_width = left_margin + total_width + right_margin

    # create half-height figure
    fig, ax = plt.subplots(figsize=(fig_width, height_half))
    if axes2: ax2 = ax.twinx()

    # plotting bars and points
    for s, xs in x_positions:
        sub = df[df[sample_name] == s]
        A10 = sub['area_change_t10'].to_numpy()
        A20 = sub['area_change_t20'].to_numpy()
        contraction = -sub['rel_max_contraction'].to_numpy()
        base_color = mcolors.to_rgba(f"#{sub['color'].iloc[0]}")
        darker = tuple(c * 0.8 for c in base_color[:3]) + (1,)
        bar_width_adj = bar_width * 0.8
        ax.bar(xs, A20, width=bar_width_adj, color=darker, edgecolor='none', linewidth=0.2)
        ax.bar(xs, A10 - A20, bottom=A20, width=bar_width_adj, color=base_color, edgecolor='none', linewidth=0.2)
        if axes2: ax2.scatter(xs, contraction, marker='o', s=2, color='gray', zorder=3)

        max_cnt = sub['contract_count'].max()          # denom for this pattern
        for x_val, top, cnt, marker in zip(xs, A10, sub['contract_count'], sub['marker']):
            inv = max_cnt / cnt                        # e.g. 8/2 → 4
            label = f"{int(inv)}" if inv.is_integer() else f"{inv:.2f}"
            #if np.floor(float(label)) == 1:
            #    ax.bar(xs, A20, width=bar_width, color=darker, edgecolor='black', linewidth=0.2)
            #    ax.bar(xs, A10 - A20, bottom=A20, width=bar_width, color=base_color, edgecolor='black', linewidth=0.2)
            if cnt == max_cnt:                         # full-active case
                ax.scatter(
                    x_val, top+0.025,
                    marker=marker, color=base_color, edgecolor='black', s=5, zorder=4,  linewidth=0.5
                )

    # axis limits and dashed zero line
    ax.set_xlim(-bar_width, total_width)
    ax.set_ylim(-0.025, 0.8)
    ax.hlines([0.0, 0.2, 0.4, 0.6], -bar_width, total_width, color='black', linewidth=0.25, linestyles='--', zorder=-1)
    if axes2: ax2.set_ylim(-0.025*0.25/0.8, 0.25)

    # tick locators
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    if axes2: ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.25/2))
    if axes2: ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.25/4))

    # xticks and labels
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=60, ha='center', fontsize=7)

    # axis labels
    if gi == 1:
        ax.set_ylabel('Rel. area change')
        if axes2: ax2.set_ylabel('Rel. max contraction', color='gray')
    else:
        ax.set_ylabel('Rel. area change')
        # hide repeat labels if desired

    if axes2: ax2.tick_params(axis='y', colors='gray')

    fig.tight_layout(pad=0.2)

    # save
    out = f'barplots_part{gi}.pdf'
    fig.savefig(out, format='pdf', dpi=300, bbox_inches='tight')
    convert_fonts_to_outlines(out, f'{out.replace(".pdf","_c.pdf")}')
    os.replace(f'{out.replace(".pdf","_c.pdf")}', out)
    plt.close(fig)
