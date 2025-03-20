import os
import sys
import json
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.neighbors import KernelDensity
from scipy.interpolate import interp1d, PchipInterpolator
from scipy.signal import find_peaks
import matplotlib.cm as cm
import imageio
import matplotlib.ticker as ticker

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 9})

color1 = "#000000"  # Black
color2 = "#8D8D8D"  # Gray

max_area = 200
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

    subprocess.run(gs_command, check=True)

def plot_from_json(json_file, pattern_folder, dump_file, split_analysis=False, bw=0.1):
    # Create directories if they don't exist
    plots_dir = os.path.join(pattern_folder, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    with open(json_file, 'r') as f:
        data = json.load(f)

    x_list = []
    density_list = []
    cmap = plt.get_cmap('turbo')

    # Prepare figures for KDE and histogram plots
    if False:
        fig, ax = plt.subplots(figsize=(8.5 / 2.54, 5.5 / 2.54))
        fig2, ax2 = plt.subplots(figsize=(8.5 / 2.54, 5.5 / 2.54))
        fig3, ax3 = plt.subplots(figsize=(8.5 / 2.54, 5.5 / 2.54))
    else:
        fig, ax = plt.subplots(figsize=(4.5 / 2.54, 5 / 2.54))
        fig2, ax2 = plt.subplots(figsize=(4.5 / 2.54, 5 / 2.54))
        fig3, ax3 = plt.subplots(figsize=(4.5 / 2.54, 5 / 2.54))

    times = []
    tot_areas = []
    if split_analysis:
        tot_areas2_lower = []  # for the lower (smaller) peak
        tot_areas2_higher = []  # for the higher (larger) peak
    else:
        tot_areas2 = []  # single representative peak

    for idx, time_data in enumerate(data):
        t = time_data['time']
        times.append(t)
        areas = np.array(time_data['areas'])

        # Plot the KDE
        areas_array = np.array(areas)
        min_area_threshold = 0.1

        # Filter areas
        filtered_areas = areas_array[areas_array >= min_area_threshold]
        filtered_areas = filtered_areas.reshape(-1, 1)
        
        proc_areas = np.log10(filtered_areas + 1)
        
        len_t = len(areas_array[areas_array >= 10.0])
        tot_areas.append(np.median(areas_array[areas_array >= 20.0]))

        # Weighted average area (weight = area)
        weighted_avg = np.sum(filtered_areas * filtered_areas) / np.sum(filtered_areas)
    
        # Convert to same log-scale your plot uses: y = log10(area + 1)
        wavg_plot_val = np.log10(weighted_avg + 1)

        # KDE setup
        x_min, x_max = 0.0, np.log10(max_area)
        x_values = np.linspace(x_min, x_max, 1000).reshape(-1, 1)

        kde = KernelDensity(kernel='epanechnikov', bandwidth=bw)
        kde.fit(proc_areas, sample_weight=proc_areas.flatten())

        # Score samples (log density)
        log_density = kde.score_samples(x_values)

        # Convert log density to density
        density = np.exp(log_density).flatten()

        peaks_idx, props = find_peaks(density, height=np.max(density)/24, prominence=np.max(density)/32)
        #print(peaks_idx)
        x_list.append(x_values.flatten())
        density_list.append(density)

        shift = 1.0 * t
        
        # Plot shaded regions
        num_bins = 40
        counts, bin_edges = np.histogram(proc_areas, bins=num_bins, density=True, weights=proc_areas)
        x_hist = np.linspace(0.0, np.log10(max_area), 1000)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        norm = mcolors.Normalize(vmin=0, vmax=np.log10(max_area))
        colors = cmap(norm(bin_centers))
        peak_size = data[-1]['time']/50

        alpha3 = 1.0

        mask = (x_hist < bin_edges[0])
        ax2.fill_betweenx(x_hist[mask], shift, shift + density[mask] * peak_size, lw=0.25, color=colors[0], alpha=alpha3, zorder=-idx, rasterized=True)
        for i, (left, right, color) in enumerate(zip(bin_edges[:-1], bin_edges[1:], colors)):
            mask = (x_hist >= left) & (x_hist < right)
            ax2.fill_betweenx(x_hist[mask], shift, shift + density[mask] * peak_size, lw=0.25, color=color, alpha=alpha3, zorder=-idx, rasterized=True)
        mask = x_hist >= bin_edges[-1]
        ax2.fill_betweenx(x_hist[mask], shift, shift + density[mask] * peak_size, lw=0.25, color=colors[-1], alpha=alpha3, zorder=-idx, rasterized=True)
        ax2.plot(shift + density * peak_size, x_values, color='black', lw=0.5, label='Epanechnikov KDE', zorder=-idx)
        
        # Handle peak extraction based on analysis mode:
        if split_analysis:
            if len(peaks_idx) == 0:
                # Fallback: if no peaks, use first index
                lower_idx = higher_idx = 0
            elif len(peaks_idx) == 1:
                # Only one peak found: duplicate it
                lower_idx = higher_idx = peaks_idx[0]
            else:
                selected_peaks = peaks_idx[-2:]
                lower_idx, higher_idx = selected_peaks[0], selected_peaks[1]
                if lower_idx < low_cut: # 675:
                    lower_idx = np.nan #higher_idx
            
            # Save the peaks converted back to original scale.
            #print(lower_idx)
            if not np.isnan(lower_idx):
                tot_areas2_lower.append(10 ** x_values[lower_idx] - 1)
                ax2.scatter(shift + density[lower_idx] * peak_size, x_values[lower_idx],
                        color=color2, marker='x', s=6, lw=0.75)
            else:
                tot_areas2_lower.append(np.nan) 
            tot_areas2_higher.append(10 ** x_values[higher_idx] - 1)
            
            # Plot each peak marker with different colors:
            ax2.scatter(shift + density[higher_idx] * peak_size, x_values[higher_idx], color=color1, marker='o', s=3, lw=0.5)
        else:
            # Default: use the last peak
            tot_areas2.append(10 ** x_values[peaks_idx[-1]] - 1)
            #ax2.scatter(shift + density[peaks_idx[-1]] * peak_size, x_values[peaks_idx[-1]], color='black', marker='o', s=1)

        #ax2.scatter(shift, np.log10(tot_areas[-1] + 1), color='red', marker='x', s=2)
        #ax2.scatter(shift, wavg_plot_val, color='black', marker='o', s=2)

    ax3.plot(times, np.zeros(len(times)), 'k--', lw=0.5, dashes=(5,10))
    rel_val = (tot_areas-tot_areas[0]) / tot_areas[0]
    if split_analysis:
        # Compute relative changes separately for the lower and higher peaks:
        
        tot_areas2_lower = np.array([x if not isinstance(x, np.ndarray) else x.item() for x in tot_areas2_lower])
        tot_areas2_higher = np.array(tot_areas2_higher)
        if True:
            rel_val2_lower = (tot_areas2_lower - tot_areas2_higher[0]) / tot_areas2_higher[0]
        else:
            rel_val2_lower = (tot_areas2_lower - tot_areas2_lower[0]) / tot_areas2_lower[0]
        rel_val2_higher = (tot_areas2_higher - tot_areas2_higher[0]) / tot_areas2_higher[0]

        ax3.plot(times, rel_val2_lower, color=color2, lw=1, zorder=-2)
        ax3.scatter(times, rel_val2_lower, marker='x', s=6, lw=0.75, color=color2, zorder=0)
        ax3.plot(times, rel_val2_higher, color=color1, lw=1, zorder=-2)
        ax3.scatter(times, rel_val2_higher, marker='o', s=3, lw=0.5, color=color1, edgecolors=color1, zorder=1)
        maxH = rel_val2_higher[-1][0] #np.max(rel_val2_higher)
        #maxL = rel_val2_lower[-1] #np.max(rel_val2_lower)
        idxL = np.where(~np.isnan(rel_val2_lower))[0][-1]
        maxL = rel_val2_lower[idxL]
        ax3.annotate(r'$\delta_\mathrm{{end}}={:.2f}$'.format(maxH),  
             xy=(times[-2], maxH+0.05),    # Point where annotation appears
             xytext=(times[-2], maxH+0.05),  # Offset text slightly above point
             fontsize=7, color=color1, ha='right', va='bottom')
        ax3.annotate(r'$\delta_\mathrm{{end}}={:.2f}$'.format(maxL),  
             xy=(times[idxL-1], maxL-0.05),    # Point where annotation appears
             xytext=(times[idxL-1], maxL-0.05),  # Offset text slightly above point
             fontsize=7, color=color2, ha='right', va='top')
    else:
        rel_val2 = (np.array(tot_areas2) - tot_areas2[0]) / tot_areas2[0]
        maxH = rel_val2[-1][0] #np.max(rel_val2_higher)
        ax3.plot(times, rel_val2, 'k', lw=1)
        ax3.scatter(times, rel_val2, color='black', marker='o', s=2)
        ax3.annotate(r'$\delta_\mathrm{{end}}={:.2f}$'.format(maxH),  
             xy=(times[-2], maxH+0.05),    # Point where annotation appears
             xytext=(times[-2], maxH+0.05),  # Offset text slightly above point
             fontsize=7, color=color1, ha='right', va='bottom')

    ax3.set_ylim(min_y, 1.1)
    #ax3.set_ylim(-0.35, 1.1)
    # Finalize histogram plot
    ax2.set_xlim(0.0, data[-1]['time'] + 1.1*np.max(density) * peak_size)
    ax2.set_ylim(0, np.log10(max_area))
    area_labels = [0, 1, 10, 100, 200]
    ytick_positions = [np.log10(a + 1) for a in area_labels]
    ytick_labels = [f"${a}$" for a in area_labels] 

    #ax3.yaxis.set_major_locator(ticker.MultipleLocator(0.05)) # Major ticks every 0.5 units
    #ax3.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
    # Then apply them to your axis:
    ax2.set_yticks(ytick_positions)
    ax2.set_yticklabels(ytick_labels)
    ax2.set_ylabel(r'Hole area $(\mathrm{mm}^2)$')

    # Prepare data for the KDE contour plot
    global_x_min = min([xv.min() for xv in x_list])
    global_x_max = max([xv.max() for xv in x_list])
    common_x_values = np.linspace(global_x_min, global_x_max, 1000)

    kde_data_interpolated = []

    for x_vals, dens in zip(x_list, density_list):
        interp_func = interp1d(x_vals, dens, bounds_error=False, fill_value=0)
        density_common = interp_func(common_x_values)
        kde_data_interpolated.append(density_common)

    kde_data = np.array(kde_data_interpolated)
    extra_steps = 5
    actual_times = np.array([time_data['time'] for time_data in data])
    t_old = actual_times
    t_new = np.linspace(actual_times.min(), actual_times.max(), len(density_list)*extra_steps - 1)

    density_interpolated = np.zeros((len(t_new), kde_data.shape[1]))

    # Perform PCHIP interpolation
    for i in range(kde_data.shape[1]):
        interpolator = PchipInterpolator(t_old, kde_data[:, i])
        density_interpolated[:, i] = interpolator(t_new)

    T_grid, X_grid = np.meshgrid(t_new, common_x_values)

    density_interpolated = np.nan_to_num(
        density_interpolated,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    # Plot KDE contour
    ax.pcolormesh(T_grid, X_grid, density_interpolated.T, shading='auto', cmap='magma')
    if False:
        ax.set_xlabel(r'Time $(\mathrm{s})$')
        ax2.set_xlabel(r'Time $(\mathrm{s})$')
        ax3.set_xlabel(r'Time $(\mathrm{s})$')
    else:
        ax.set_xlabel(r'Simulation Step')
        ax2.set_xlabel(r'Simulation Step')
        ax3.set_xlabel(r'Simulation Step')
    ax.set_ylabel(r'Hole area $(\mathrm{mm}^2)$')
    ax.set_xlim(t_new.min(), t_new.max())
    ax.set_ylim(0, np.log10(max_area))
    ax.set_yticks([0, 1, 2, np.log10(max_area)])
    ax.set_yticklabels(['1', '10', '100', str(max_area)])
    
    
    ax3.set_ylabel(r'$\delta = \frac{A - A_0}{A_0}$')
    ax3.set_xlim(0.0, data[-1]['time'])

    # Save the plots
    fig_filename = os.path.join(plots_dir, dump_file + '_kde.png')
    fig2_filename = os.path.join(plots_dir, dump_file + '_hist.pdf')
    fig3_filename = os.path.join(plots_dir, dump_file + '_relArea.pdf')
    fig4_filename = os.path.join(plots_dir, dump_file + '_relArea.png')
    fig.savefig(fig_filename, format='png', dpi=300, bbox_inches='tight')
    fig2.savefig(fig2_filename, format='pdf', dpi=300, bbox_inches='tight')
    fig3.savefig(fig3_filename, format='pdf', dpi=300, bbox_inches='tight')
    ax3.plot(times, rel_val, 'r', lw = 1, zorder=-1)
    fig3.savefig(fig4_filename, format='png', dpi=300, bbox_inches='tight')

    # Compress and convert fonts to outlines for each PDF
    for pdf_file in [fig2_filename, fig3_filename]:
        compressed_pdf = pdf_file.replace(".pdf", "_c.pdf")
        convert_fonts_to_outlines(pdf_file, compressed_pdf)

        # Optional: Replace original file with the compressed version
        os.replace(compressed_pdf, pdf_file)

if __name__ == '__main__':
    split_analysis = False
    bandwidth = 0.1  # Default bandwidth
    min_y = -0.1
    if len(sys.argv) < 3 or len(sys.argv) > 7:
        print("Usage: python holes_analysis.py <pattern_folder> <dump_file> [split] [bandwidth] [lowcut] [plot y min]")
        sys.exit(1)

    if len(sys.argv) >= 4:
        split_analysis = sys.argv[3].lower() in ("true", "1", "yes", "0", "-0.1")
    
    # Check if the fifth argument exists, otherwise use default 0.1
    if len(sys.argv) >= 5:
        try:
            bandwidth = float(sys.argv[4])  # Convert to float
            low_cut = int(sys.argv[5])
            min_y = float(sys.argv[6])
        except ValueError:
            print(f"Error: Invalid bandwidth '{sys.argv[4]}'. Must be a number.")
            sys.exit(1)

    pattern_folder = sys.argv[1]
    dump_file = sys.argv[2]

    # Construct the JSON file path
    json_file = os.path.join(pattern_folder, 'plots', dump_file + '_holes.json')

    # Check if the JSON file exists
    if not os.path.exists(json_file):
        print(f"Error: JSON file '{json_file}' does not exist.")
        sys.exit(1)

    # Call the plotting function
    plot_from_json(json_file, pattern_folder, dump_file, split_analysis, bandwidth)