import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon as ShapelyPolygon
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors
from sklearn.neighbors import KernelDensity
from scipy.interpolate import interp1d, PchipInterpolator
import matplotlib.cm as cm
import imageio

plt.rcParams.update({'font.size': 20})

def plot_from_json(json_file, pattern_folder, dump_file, method='log10'):
    # Create directories if they don't exist
    plots_dir = os.path.join(pattern_folder, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    with open(json_file, 'r') as f:
        data = json.load(f)

    frames = []
    x_list = []
    density_list = []
    cmap = plt.get_cmap('jet')

    # Prepare figures for KDE and histogram plots
    fig, ax = plt.subplots(figsize=(8, 8))
    fig2, ax2 = plt.subplots(figsize=(8, 8))

    for idx, time_data in enumerate(data):
        t = time_data['time']
        areas = np.array(time_data['areas'])

        if t == 0:
            max_area = 1000  # Or set to max(areas) if preferred

        # Plot the KDE
        areas_array = np.array(areas)
        min_area_threshold = 1.0

        # Filter areas
        filtered_areas = areas_array[areas_array >= min_area_threshold]
        filtered_areas = filtered_areas.reshape(-1, 1)

        if method == 'log10':
            proc_areas = np.log10(filtered_areas)
        elif method == 'sqrt':
            proc_areas = np.sqrt(filtered_areas)
        else:
            proc_areas = filtered_areas

        # KDE setup
        x_min, x_max = 0.0, 3.0
        x_values = np.linspace(x_min, x_max, 1000).reshape(-1, 1)

        if method == 'log10':
            kde = KernelDensity(kernel='epanechnikov', bandwidth=0.05)
        elif method == 'sqrt':
            kde = KernelDensity(kernel='epanechnikov', bandwidth=1.0)
        else:
            kde = KernelDensity(kernel='epanechnikov', bandwidth=0.1)

        kde.fit(proc_areas, sample_weight=proc_areas.flatten())

        # Score samples (log density)
        log_density = kde.score_samples(x_values)

        # Convert log density to density
        density = np.exp(log_density).flatten()
        x_list.append(x_values.flatten())
        density_list.append(density)

        if method == 'log10':
            shift = 1.0 * t
        elif method == 'sqrt':
            shift = 0.1 * t
        else:
            shift = t

        # Plot shaded regions
        num_bins = 40
        counts, bin_edges = np.histogram(proc_areas, bins=num_bins, density=True, weights=proc_areas)
        x_hist = np.linspace(0.0, 3.0, 1000)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        norm = mcolors.Normalize(vmin=0, vmax=3.0)
        colors = cmap(norm(bin_centers))
        peak_size = 0.5
        alpha3 = 1.0

        mask = (x_hist < bin_edges[0])
        ax2.fill_betweenx(x_hist[mask], shift, shift + density[mask] * peak_size, color=colors[0], alpha=alpha3, zorder=-int(t))
        for i, (left, right, color) in enumerate(zip(bin_edges[:-1], bin_edges[1:], colors)):
            mask = (x_hist >= left) & (x_hist < right)
            ax2.fill_betweenx(x_hist[mask], shift, shift + density[mask] * peak_size, color=color, alpha=alpha3, zorder=-int(t))
        mask = x_hist >= bin_edges[-1]
        ax2.fill_betweenx(x_hist[mask], shift, shift + density[mask] * peak_size, color=colors[-1], alpha=alpha3, zorder=-int(t))
        ax2.plot(shift + density * peak_size, x_values, color='black', lw=1, label='Epanechnikov KDE', zorder=-int(t-1))

    # Finalize histogram plot
    ax2.set_xlim(0.0, len(density_list) + 0.5)
    ax2.set_ylim(0, 3.0)
    ax2.set_yticks([0, 1, 2, 3])
    ax2.set_yticklabels(['1', '10', '100', '1000'])
    ax2.set_xlabel('Time Step')
    ax2.set_ylabel('Area')

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
    t_old = np.arange(len(density_list))
    extra_steps = 5
    t_new = np.linspace(t_old[0], t_old[-1], len(density_list) * extra_steps - 1)

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
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Area')
    ax.set_xlim(t_new.min(), t_new.max())
    ax.set_ylim(0, 3.0)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(['1', '10', '100', '1000'])

    # Save the plots
    fig_filename = os.path.join(plots_dir, dump_file + '_kde.png')
    fig2_filename = os.path.join(plots_dir, dump_file + '_hist.png')
    fig.savefig(fig_filename, format='png', dpi=300, bbox_inches='tight')
    fig2.savefig(fig2_filename, format='png', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python holes_analysis.py <pattern_folder> <dump_file>")
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
    plot_from_json(json_file, pattern_folder, dump_file, method='log10')