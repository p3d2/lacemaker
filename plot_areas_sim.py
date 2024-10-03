# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 17:11:00 2024

@author: silvap1
"""

import argparse
import os

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import imageio.v2 as imageio
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d, PchipInterpolator
from shapely.geometry import Point
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union
from sklearn.neighbors import KernelDensity

fg_color = 'black'
bg_color = 'none'

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans']
mpl.rcParams['font.size'] = 20

mpl.rcParams['text.color'] = fg_color
mpl.rcParams['axes.labelcolor'] = fg_color
mpl.rcParams['xtick.color'] = fg_color
mpl.rcParams['ytick.color'] = fg_color
mpl.rcParams['axes.edgecolor'] = fg_color  # Changes the edge color of the axes

mpl.rcParams['figure.facecolor'] = bg_color
mpl.rcParams['axes.facecolor'] = bg_color

# Create the colormap
cmap = plt.cm.turbo

# Function to normalize and get color
def get_color(value, vmin=0, vmax=0.675):
    # Normalize value
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    return cmap(norm(value))

def read_lammps_dump(filename, xmin, xmax, ymin, ymax):
    with open(filename, 'r') as file:
        lines = file.readlines()

    forces = []
    pe_atom = []
    pos_atom = []
    num_atoms = 0
    coord_atoms = []
    id_atom = []
    atom_type = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if 'ITEM: NUMBER OF ATOMS' in line:
            num_atoms = int(lines[i + 1])
        elif 'ITEM: ATOMS' in line:
            # ITEM: ATOMS id mol type x y z fx fy fz c_perAtomPE
            i += 1  # Start reading atom data from the next line
            timestep_data = []
            pe = []
            pos = []
            ids = []
            a_type = []
            coord = []
            for _ in range(num_atoms):
                if i >= len(lines):
                    break
                parts = lines[i].split()
                if len(parts) >= 10:  # Ensure there are enough parts in the line
                    ids.append(int(parts[0]))
                    mol = int(parts[1])
                    a_type.append(int(parts[2]))
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    fx, fy, fz = float(parts[6]), float(parts[7]), float(parts[8])
                    force_magnitude = np.sqrt(fx**2 + fy**2 + fz**2)
                    timestep_data.append(force_magnitude)

                    if xmin < x < xmax and ymin < y < ymax:
                        pe_value = float(parts[9])
                    else:
                        pe_value = 0  # Assign zero if outside ROI

                    pe.append(pe_value)
                    coord.append([x, y, z])
                i += 1
            pos_atom.append(ids)
            id_atom.append(ids)
            atom_type.append(a_type)
            forces.append(timestep_data)
            pe_atom.append(pe)
            coord_atoms.append(coord)
            continue
        i += 1

    return forces, np.array(id_atom), np.array(atom_type), np.array(pe_atom), np.array(pos_atom), coord_atoms

# Function to extract polygons from the merged geometry
def get_polygons(geom):
    if geom.geom_type == 'Polygon':
        return [geom]
    elif geom.geom_type == 'MultiPolygon':
        return list(geom.geoms)
    else:
        return []

def save_plot(positions):
    fig, ax = plt.subplots(figsize=(8, 8))
    fig3, ax3 = plt.subplots(figsize=(8, 8))

    method = 'log10'
    density_list, x_list = [], []
    frames = []
    for t in range(len(positions)):
        fig2, ax2 = plt.subplots(figsize=(8, 8))

        pos_0 = np.array(positions[t]) 
        x = pos_0[:, 0]
        y = pos_0[:, 1]
        points = np.column_stack((x, y))
        radius = 0.5

        circles = [Point(px, py).buffer(radius) for px, py in points]

        # Merge circles with union
        merged = unary_union(circles)
        polygons = get_polygons(merged)

        # Prepare matplotlib patches
        patches = []
        for poly in polygons:
            exterior_coords = np.array(poly.exterior.coords)
            patches.append(Polygon(exterior_coords))

        # Create a PatchCollection
        collection = PatchCollection(patches, facecolor='grey', edgecolor='black', alpha=0.5)

        # Plotting
        ax2.add_collection(collection)

        areas = []
        for poly in polygons:
            for interior in poly.interiors:
                interior_polygon = ShapelyPolygon(interior)
                areas.append(interior_polygon.area)

        if t == 0:
            max_area = 1000 #max(areas)

        for poly in polygons:
            for (interior, area) in zip(poly.interiors, areas):
                interior_coords = np.array(interior.coords)
                # Use the same color or modify as needed
                if method == 'log10': color = cmap(np.log10(area + 1) / np.log10(max_area + 1))
                if method == 'sqrt': color = cmap(np.sqrt(area/max_area))
                hole = Polygon(
                    interior_coords,
                    facecolor=color,
                    edgecolor='gray',
                    linewidth=1,
                    alpha=0.7
                )
                ax2.add_patch(hole)

        # Set plot limits
        ax2.set_xlim(x.min() - radius, x.max() + radius)
        ax2.set_ylim(y.min() - radius, y.max() + radius)
        ax2.set_aspect('equal')  # Ensure the aspect ratio is equal
        frame_filename = os.path.join('output', 'simulations', pattern_folder, 'plots', f'frame_{t}.png')
        fig2.savefig(frame_filename, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
        frames.append(frame_filename)
        plt.close(fig2)
        # Plot the KDE
        areas_array = np.array(areas)
        min_area_threshold = 1.0  # Adjust this value as needed

        # Filter areas to include only those above the threshold
        filtered_areas = areas_array[areas_array >= min_area_threshold]
        
        filtered_areas = filtered_areas.reshape(-1, 1)  # Reshape for sklearn apply log to the data
        if method == 'log10': proc_areas = np.log10(filtered_areas) # Use log data
        if method == 'sqrt': proc_areas = np.sqrt(filtered_areas)
        
        # Define the range for the KDE plot
        x_min, x_max = 0.0, 3.0
        x_values = np.linspace(x_min, x_max, 1000).reshape(-1, 1)

        # Initialize KernelDensity
        if method == 'log10': kde = KernelDensity(kernel='epanechnikov', bandwidth=0.05)
        if method == 'sqrt': kde = KernelDensity(kernel='epanechnikov', bandwidth=1.0)  # Adjust bandwidth as needed
        
        kde.fit(proc_areas, sample_weight=proc_areas.flatten())
        
        # Score samples (log density)
        log_density = kde.score_samples(x_values)

        # Convert log density to density
        density = np.exp(log_density).flatten()
        x_list.append(x_values.flatten())
        density_list.append(density)

        if method == 'log10': shift = 1.0 * t
        if method == 'sqrt': shift = 0.1 * t  # To make a waterfall plot
        
        # Plot shaded regions:
        num_bins = 40
        counts, bin_edges = np.histogram(proc_areas, bins=num_bins, density=True, weights=proc_areas)
        x = np.linspace(0.0, 3.0, 1000)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        norm = mcolors.Normalize(vmin=0, vmax=3.0)
        colors = cmap(norm(bin_centers))
        peak_size = 0.5
        #ax3.bar(bin_centers, counts, width=(bin_edges[1] - bin_edges[0]), bottom=-shift,
        #       color=colors, edgecolor='black', alpha=0.5)

        mask = (x < bin_edges[0])
        alpha3 = 1.0
        ax3.fill_betweenx(x[mask], shift,  shift + density[mask]*peak_size,  color=colors[0], alpha=alpha3)
        for i, (left, right, color) in enumerate(zip(bin_edges[:-1], bin_edges[1:], colors)):
            mask = (x >= left) & (x < right)
            ax3.fill_betweenx(x[mask], shift,  shift + density[mask]*peak_size,  color=color, alpha=alpha3)
        mask = x >= bin_edges[-1]
        ax3.fill_betweenx(x[mask], shift,  shift + density[mask]*peak_size,  color=colors[-1], alpha=alpha3)
        ax3.plot(shift + density*peak_size, x_values,  color='black', lw=1, label='Epanechnikov KDE')
    ax3.set_xlim(0.0, len(density_list)+0.5)
    ax3.set_yticks([0, 1, 2, 3], ['1', '10', '100', '1000'])
    ax3.set_xlabel('Time Step')
    ax3.set_ylabel('Area')
    # Define the common x_values grid
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
    t_new = np.linspace(t_old[0], t_old[-1], len(density_list)*extra_steps-1)

    density_interpolated = np.zeros((len(t_new), kde_data.shape[1]))

    # Perform PCHIP interpolation for each x_value
    for i in range(kde_data.shape[1]):
        interpolator = PchipInterpolator(t_old, kde_data[:, i])
        density_interpolated[:, i] = interpolator(t_new)

    num_time_steps = density_interpolated.shape[0]
    num_x_values = density_interpolated.shape[1]
    num_x_ticks = 10  
    num_y_ticks = 10

    T_grid, X_grid = np.meshgrid(t_new, common_x_values)

    density_interpolated = np.nan_to_num(
        density_interpolated,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    # Use contourf to plot the data
    ax.pcolormesh(T_grid, X_grid, density_interpolated.T, shading='auto', cmap='magma')
    ax.invert_yaxis()
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Area')

    # Optionally, set the limits of the axes
    ax.set_xlim(t_new.min(), t_new.max())
    ax3.set_yticks([0, 1, 2, 3], ['1', '10', '100', '1000'])

    fig.savefig(os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_kde.png'),
                format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    fig3.savefig(os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_hist.png'),
                format='png', dpi=300, bbox_inches='tight', pad_inches=0)

    gif_filename = os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_area_dist.gif')
    with imageio.get_writer(gif_filename, mode='I', fps=10, loop=0) as writer:
        for filename in frames:
            image = imageio.imread(filename)
            writer.append_data(image)

    for filename in frames:
        os.remove(filename)

parser = argparse.ArgumentParser(description='Processing dump file name')
parser.add_argument('folder', type=str, help='The name of the dump file to be processed')
parser.add_argument('trj_file', type=str, help='The name of the dump file to be processed')
parser.add_argument('--mol_r', nargs="+", type=int, help='Select the yarns (mol id)')
parser.add_argument('--roi', nargs="+", type=float, help='ROI to consider PE')
parser.add_argument('--yrange', nargs="+", type=float, help='Range for y plots')

args = parser.parse_args()
pattern_folder = args.folder
dump_file = args.trj_file
mol_r = args.mol_r
xmin, xmax, ymin, ymax = args.roi

(force_data, id_data, atom_type, pe_data, pos, coord_atoms) = read_lammps_dump(
    os.path.join('output', 'simulations', pattern_folder, dump_file + '.lammpstrj'),
    xmin, xmax, ymin, ymax)

for k in range(len(pe_data)):
    pe = pe_data[k]
    id_atom = id_data[k]
    type_atom = atom_type[k]
    id_flat = pos[k]
    c_flat = np.array(coord_atoms[k])
    # Get the sort indices based on the flattened ID matrix
    sort_indices = np.argsort(id_flat)

    # Use the indices to sort arrays
    pe_data[k] = pe[sort_indices]
    id_data[k] = id_atom[sort_indices]
    atom_type[k] = type_atom[sort_indices]
    pos[k] = id_flat[sort_indices]
    coord_atoms[k] = c_flat[sort_indices]

len_data, pe_sort, id_new = [], [], []
for k in range(len(pe_data)):
    len_temp, len_yarn = [], []
    pe_temp, pe_yarn = [], []
    id_temp, id_yarn = [], []

    for l in range(len(pe_data[k]) - 1):
        if id_data[k][l] != id_data[k][l + 1]:
            len_yarn.append(len_temp)
            len_temp = []
            pe_yarn.append(pe_temp)
            pe_temp = []
            id_yarn.append(id_temp)
            id_temp = []
        else:
            p1 = coord_atoms[k][l]
            p2 = coord_atoms[k][l + 1]
            len_temp.append(np.linalg.norm(p2 - p1))
            pe_temp.append(pe_data[k][l])
            id_temp.append(id_data[k][l])
    len_data.append(len_yarn)
    pe_sort.append(pe_yarn)
    id_new.append(id_yarn)

save_plot(coord_atoms)