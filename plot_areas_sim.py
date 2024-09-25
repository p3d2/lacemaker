# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 17:11:00 2024

@author: silvap1
"""

import argparse, json, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from scipy.fft import fft
from scipy.interpolate import interp1d, RectBivariateSpline, RBFInterpolator
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#from KDEpy import FFTKDE
from scipy.stats import norm
from scipy.spatial import Delaunay
from shapely.geometry import Point
from shapely import geometry, ops
from shapely.ops import unary_union
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm
import matplotlib.colors as mcolors
from shapely.strtree import STRtree
from shapely import vectorized
from shapely.geometry import Polygon as ShapelyPolygon
import seaborn as sns
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity

###patch start###
from mpl_toolkits.mplot3d.axis3d import Axis
if not hasattr(Axis, "_get_coord_info_old"):
    def _get_coord_info_new(self, renderer):
        mins, maxs, centers, deltas, tc, highs = self._get_coord_info_old(renderer)
        mins += deltas / 4
        maxs -= deltas / 4
        return mins, maxs, centers, deltas, tc, highs
    Axis._get_coord_info_old = Axis._get_coord_info  
    Axis._get_coord_info = _get_coord_info_new
###patch end###

fg_color = 'black'
bg_color = 'none'

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans']
mpl.rcParams['font.size'] = 20

mpl.rcParams['text.color'] = fg_color
mpl.rcParams['axes.labelcolor'] = fg_color
mpl.rcParams['xtick.color'] = fg_color
mpl.rcParams['ytick.color'] = fg_color
mpl.rcParams['axes.labelcolor'] = fg_color
mpl.rcParams['axes.edgecolor'] = fg_color  # Changes the edge color of the axes

# Optionally, set a dark background for better contrast
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
            coord = []
            id = []
            a_type = []
            for _ in range(num_atoms):
                if i >= len(lines):
                    break
                parts = lines[i].split()
                if len(parts) >= 3:  # Ensure there are enough parts in the line
                    pos.append(float(parts[0]))
                    id.append(int(parts[1]))
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
            pos_atom.append(pos)
            id_atom.append(id)
            atom_type.append(a_type)
            forces.append(timestep_data)
            pe_atom.append(pe)
            coord_atoms.append(coord)
            continue
        i += 1

    return (forces, np.array(id_atom), np.array(atom_type), np.array(pe_atom), np.array(pos_atom), coord_atoms)

# Function to extract polygons from the merged geometry
def get_polygons(geom):
    if geom.geom_type == 'Polygon':
        return [geom]
    elif geom.geom_type == 'MultiPolygon':
        return list(geom.geoms)
    else:
        return []

def save_plot(positions):

    fig, ax = plt.subplots(figsize=(12, 8))
    fig2, ax2 = plt.subplots(figsize=(12, 8))
    for t in range(len(pos)):
    #for t in range(15,16):
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
        ##fig, ax = plt.subplots(figsize=(12, 12))
        ax2.add_collection(collection)

        areas = []
        for poly in polygons:
            for interior in poly.interiors:
                interior_polygon = ShapelyPolygon(interior)
                areas.append(interior_polygon.area)

        if t == 0:
            max_area = max(areas) 
        cmap = mpl.colormaps['turbo']

        for poly in polygons:
            for (interior, area) in zip(poly.interiors, areas):
                interior_coords = np.array(interior.coords)
                # Use the same color or modify as needed
                color = cmap(area/max_area)
                hole = Polygon(
                    interior_coords,
                    facecolor=color,  # Or assign a different mapping if desired
                    edgecolor='white',
                    linestyle='dashed',
                    linewidth=1,
                    alpha=0.7
                )
                ax2.add_patch(hole)

        # Set plot limits
        ax2.set_xlim(x.min() - radius, x.max() + radius)
        ax2.set_ylim(y.min() - radius, y.max() + radius)
        ax2.set_aspect('equal')  # Ensure the aspect ratio is equal
        if t ==0:
            fig2.savefig(os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_area_dist.png'), format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    
    
        # Plot the KDE using Seaborn with Gaussian kernel and Silverman's bandwidth
        areas_array = np.array(areas)
        min_area_threshold = 1.0  # Adjust this value as needed

        # Filter areas to include only those above the threshold
        filtered_areas = areas_array[areas_array >= min_area_threshold]
        # areas_array = filtered_areas
        
        kde = gaussian_kde(areas_array, bw_method='silverman')
        x_min, x_max = 10.0, areas_array.max()
        x_values = np.linspace(x_min, x_max, 1000)
        density = kde(x_values)

        # Plot the KDE
        #ax.plot(x_values, density - 0.002*t, color='blue')
        #ax.fill_between(x_values, density - 0.001*t, alpha=0.5, color='skyblue')
        
        filtered_areas = np.array(filtered_areas).reshape(-1, 1)  # Reshape for sklearn

        # Define the range for the KDE plot
        x_min, x_max = filtered_areas.min(), filtered_areas.max()
        x_values = np.linspace(x_min, x_max, 1000).reshape(-1, 1)

        # Initialize KernelDensity with a non-Gaussian kernel, e.g., 'epanechnikov'
        kde = KernelDensity(kernel='epanechnikov', bandwidth=5.0)  # Adjust bandwidth as needed

        kde.fit(filtered_areas, sample_weight=filtered_areas.flatten())

        # Score samples (log density)
        log_density = kde.score_samples(x_values)

        # Convert log density to density
        density = np.exp(log_density)

        shift = 0.01 * t # To make a waterfall plot
        # Plot the KDE
        ax.plot(x_values, density-shift, color='black', lw=1, label='Epanechnikov KDE')
        
        counts, bin_edges = np.histogram(filtered_areas, bins=100, density=True, weights=filtered_areas)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        norm = mcolors.Normalize(vmin=0, vmax=max_area)
        colors = cmap(norm(bin_centers))
        ax.bar(bin_centers, counts, width=(bin_edges[1] - bin_edges[0]), bottom=-shift,
                color=colors, edgecolor='black', alpha=0.5)
        #ax.scatter(bin_centers, counts - 0.0025*t, alpha=0.5)

        # Add titles and labels
    plt.title('KDE of Interior Hole Areas')
    plt.xlabel('Area')
    plt.ylabel('Density')
    
    fig.savefig(os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_kde.png'), format='png', dpi=300, bbox_inches='tight', pad_inches=0)

parser = argparse.ArgumentParser(description='Processing dump file name')
parser.add_argument('folder', type=str, help='The name of the dump file to be processed')
parser.add_argument('trj_file', type=str, help='The name of the dump file to be processed')
parser.add_argument('--mol_r', nargs="+", type=int, help='Select the yarns (mol id)')
parser.add_argument('--roi', nargs="+", type=float, help='Roi to consider PE')
parser.add_argument('--yrange', nargs="+", type=float, help='Range for y plots')

args = parser.parse_args()
pattern_folder = args.folder
dump_file = args.trj_file
mol_r = args.mol_r
xmin, xmax, ymin, ymax = args.roi

(force_data, id_data, atom_type, pe_data, pos, coord_atoms) = read_lammps_dump(os.path.join('output','simulations', pattern_folder, dump_file + '.lammpstrj'), xmin, xmax, ymin, ymax)

for k in range(len(pe_data)):
    pe = pe_data[k]
    id_atom = id_data[k]
    type_atom = atom_type[k]
    id_flat = pos[k]
    c_flat = np.array(coord_atoms[k])
    # Get the sort indices based on the flattened ID matrix
    sort_indices = np.argsort(id_flat)

    # Use the indices to sort both arrays
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

    for l in range(len(pe_data[k])-1):
        if id_data[k][l] != id_data[k][l+1]:
            len_yarn.append(len_temp)
            len_temp = []
            pe_yarn.append(pe_temp)
            pe_temp = []
            id_yarn.append(id_temp)
            id_temp = []
        else:
            p1 = coord_atoms[k][l]
            p2 = coord_atoms[k][l+1]
            len_temp.append(np.linalg.norm(p2-p1))
            pe_temp.append(pe_data[k][l])
            id_temp.append(id_data[k][l])
    len_data.append(len_yarn)
    pe_sort.append(pe_yarn)  
    id_new.append(id_yarn)

save_plot(coord_atoms)