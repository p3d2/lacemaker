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
import json
import imageio.v2 as imageio
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from shapely.geometry import Point
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union

fg_color = 'black'
bg_color = 'none'

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans']
mpl.rcParams['font.size'] = 20

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

def save_holes(positions):
    
    data = []
    frames = []
    method = 'log10'
    for t in range(len(positions)):

        pos_0 = np.array(positions[t]) 
        x = pos_0[:, 0]
        y = pos_0[:, 1]
        points = np.column_stack((x, y))
        radius = 0.5

        circles = [Point(px, py).buffer(radius) for px, py in points]

        # Merge circles with union
        merged = unary_union(circles)
        polygons = get_polygons(merged)

        # Draw gif
        patches = []
        for poly in polygons:
            exterior_coords = np.array(poly.exterior.coords)
            patches.append(Polygon(exterior_coords))

        collection = PatchCollection(patches, facecolor='grey', edgecolor='black', alpha=0.5)
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_facecolor('white')
        ax.add_collection(collection)
        max_area = 1000
        areas = []
        for poly in polygons:
            exterior_coords = np.array(poly.exterior.coords)
            patches.append(Polygon(exterior_coords))
            for interior in poly.interiors:
                interior_coords = np.array(interior.coords)
                area = ShapelyPolygon(interior_coords).area
                areas.append(area)
                if method == 'log10':
                    color = cmap(np.log10(area + 1) / np.log10(max_area))
                elif method == 'sqrt':
                    color = cmap(np.sqrt(area / max_area))
                else:
                    color = 'grey'
                hole = Polygon(interior_coords, facecolor=color, edgecolor='grey', linewidth=1, alpha=1.0)
                ax.add_patch(hole)

        ax.set_xlim(x.min() - radius, x.max() + radius)
        ax.set_ylim(y.min() - radius, y.max() + radius)

        frame_filename = os.path.join(pattern_folder, 'plots', dump_file + f'frame_{t}.png')
        fig.savefig(frame_filename, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
        frames.append(frame_filename)
        plt.close(fig)

        # Save time and areas only
        time_data = {
            'time': t,
            'areas': areas  # List of areas of polygons at time t
        }
        data.append(time_data)
        
    # Write data to a JSON file
    output_file = os.path.join(pattern_folder, 'plots', dump_file + '_holes.json')
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    # Generate GIF from frames
    gif_filename = os.path.join(pattern_folder, 'plots', dump_file + '_area_dist.gif')
    with imageio.get_writer(gif_filename, mode='I', fps=10, loop=0) as writer:
        for filename in frames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Clean up temporary frame files
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
    os.path.join(pattern_folder, dump_file + '.lammpstrj'),
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

save_holes(coord_atoms)