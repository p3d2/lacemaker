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
from shapely.geometry import Point, LineString
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union
from scipy.optimize import linear_sum_assignment

fg_color = 'black'
bg_color = 'none'

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans']
mpl.rcParams['font.size'] = 20

MIN_HOLE_AREA = 1.0  # area threshold
MAX_CENTROID_DIST = 2.5  # tweak this distance threshold to match holes from frame to frame
radius = 0.1# yarn radiues

# Create the colormap
cmap = plt.cm.turbo

def read_lammps_dump(filename, xmin, xmax, ymin, ymax):
    """
    Reads a LAMMPS dump file that has lines of format:
        ITEM: ATOMS id mol type x y z fx fy fz c_perAtomPE
    Returns arrays for forces, id, mol, atom_type, pe, coord_atoms.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    forces = []
    pe_atom = []
    num_atoms = 0
    coord_atoms = []
    id_atom = []
    mol_atom = [] 
    atom_type = []

    i = 0
    while i < len(lines):
        line = lines[i]
        if 'ITEM: NUMBER OF ATOMS' in line:
            num_atoms = int(lines[i + 1])
        elif 'ITEM: ATOMS' in line:
            i += 1  # Start reading atom data
            timestep_forces = []
            timestep_pe = []
            ids = []
            mol_frame = []
            a_type = []
            coords = []

            for _ in range(num_atoms):
                if i >= len(lines):
                    break
                parts = lines[i].split()
                if len(parts) >= 10:
                    # Read atom data
                    aid = int(parts[0])
                    mol_id = int(parts[1])
                    typ = int(parts[2])
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    fx, fy, fz = float(parts[6]), float(parts[7]), float(parts[8])
                    force_magnitude = np.sqrt(fx**2 + fy**2 + fz**2)

                    # Potential energy only if within ROI
                    if xmin < x < xmax and ymin < y < ymax:
                        pe_value = float(parts[9])
                    else:
                        pe_value = 0

                    # Store in arrays
                    ids.append(aid)
                    mol_frame.append(mol_id)
                    a_type.append(typ)
                    timestep_forces.append(force_magnitude)
                    timestep_pe.append(pe_value)
                    coords.append([x, y, z])

                i += 1

            # Done reading one timestep's atoms
            id_atom.append(ids)
            mol_atom.append(mol_frame)
            atom_type.append(a_type)
            forces.append(timestep_forces)
            pe_atom.append(timestep_pe)
            coord_atoms.append(coords)
            continue
        i += 1

    return (
        np.array(forces),
        np.array(id_atom),
        np.array(mol_atom),
        np.array(atom_type),
        np.array(pe_atom),
        np.array(coord_atoms),
    )

def get_polygons(geom):
    """Handle single or multi-polygon merges in shapely."""
    if geom.geom_type == 'Polygon':
        return [geom]
    elif geom.geom_type == 'MultiPolygon':
        return list(geom.geoms)
    else:
        return []

def save_holes(positions, pattern_folder, dump_file, roi_area):
    """
    Builds circle-buffer polygons around each atom coordinate,
    merges them, and extracts hole areas. Plots them and writes:
      - a GIF showing holes
      - a JSON with hole (normalized) areas at each frame
    """

    hole_labels = {}
    next_label_id = 1

    data = []
    frames = []
    method = 'log10'
    max_area = 1000.0  # for color scale in plot

    for t in range(len(positions)):
        pos_0 = np.array(positions[t])
        x = pos_0[:, 0]
        y = pos_0[:, 1]

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_facecolor('white')

        #circles = [Point(px, py).buffer(radius, resolution=64) for px, py in pos_0[:, :2]]
        #merged = unary_union(circles)
        line_segments = []
        for i in range(len(pos_0) - 1):
            p1 = (pos_0[i][0], pos_0[i][1])
            p2 = (pos_0[i+1][0], pos_0[i+1][1])
            if np.linalg.norm(np.array(p2) - np.array(p1)) < 1:
                segment = LineString([p1, p2]).buffer(radius)
                line_segments.append(segment)
        merged = unary_union(line_segments)
        polygons = get_polygons(merged)

        areas = []
        holes_this_frame = []  # each hole => {'area': val, 'centroid': (cx, cy)}
        patches = []
        for poly in polygons:
            exterior_coords = np.array(poly.exterior.coords)
            # Possibly add a patch for the exterior if you wish
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

                # Always use the same edge color (remove any thresholding on area)
                edge_color = 'grey'
                hole_patch = Polygon(interior_coords, facecolor=color, edgecolor=edge_color, linewidth=1, alpha=1.0)
                ax.add_patch(hole_patch)

                # Always add every hole (remove filtering by MIN_HOLE_AREA)
                centroid = ShapelyPolygon(interior_coords).centroid
                hole_info = {
                    'area': area,
                    'centroid': (centroid.x, centroid.y),
                    'coords': interior_coords
                }
                holes_this_frame.append(hole_info)

        # If you want the merged region (the big grey region) displayed:
        collection = PatchCollection(patches, facecolor='grey', edgecolor='black', alpha=0.5, zorder = -1)
        ax.add_collection(collection)

        # Now that we have all holes >= MIN_HOLE_AREA for this frame:
        # (A) If t == 0, assign new labels to all
        # (B) If t > 0, match holes to existing labels if centroids are close enough
        frame_holes_with_labels = []  # store { 'label': X, 'area': Y }

        if t == 0 or len(hole_labels) == 0:
            # First frame (or no previous labels): assign new labels to all holes
            for hole_info in holes_this_frame:
                hole_labels[next_label_id] = hole_info['centroid']
                frame_holes_with_labels.append({
                    'label': next_label_id,
                    'area': hole_info['area'],
                    'centroid': hole_info['centroid']
                })
                next_label_id += 1
        else:
            # For subsequent frames, build a cost matrix between previous holes and current holes
            prev_labels = list(hole_labels.items())  # each item is (label, centroid)
            num_prev = len(prev_labels)
            num_curr = len(holes_this_frame)
            cost_matrix = np.zeros((num_prev, num_curr))
            for i, (label, prev_centroid) in enumerate(prev_labels):
                for j, hole_info in enumerate(holes_this_frame):
                    curr_centroid = hole_info['centroid']
                    cost_matrix[i, j] = np.sqrt((prev_centroid[0] - curr_centroid[0])**2 +
                                                (prev_centroid[1] - curr_centroid[1])**2)

            # Solve the assignment problem using the Hungarian algorithm
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            
            new_labels = {}
            assigned_curr = set()
            # For each matched pair, use the old label if the distance is acceptable
            for i, j in zip(row_ind, col_ind):
                if cost_matrix[i, j] <= MAX_CENTROID_DIST:
                    label = prev_labels[i][0]
                    new_labels[label] = holes_this_frame[j]['centroid']
                    frame_holes_with_labels.append({
                        'label': label,
                        'area': holes_this_frame[j]['area'],
                        'centroid': holes_this_frame[j]['centroid']
                    })
                    assigned_curr.add(j)
            # For any current hole not matched, assign a new label
            for j, hole_info in enumerate(holes_this_frame):
                if j not in assigned_curr:
                    label = next_label_id
                    new_labels[label] = hole_info['centroid']
                    frame_holes_with_labels.append({
                        'label': label,
                        'area': hole_info['area'],
                        'centroid': hole_info['centroid']
                    })
                    next_label_id += 1
            hole_labels = new_labels

        for hole_dict in frame_holes_with_labels:
            cx, cy = hole_dict['centroid']
            lbl = str(hole_dict['label'])
            ax.text(cx, cy, lbl,
                    fontsize=6,       # small font
                    color='black',    # text color (adjust if needed)
                    ha='center', 
                    va='center')

        # Store data in JSON
        time_data = {
            'time': t,
            'areas': areas,  # all holes
            'areas_label': frame_holes_with_labels
        }
        data.append(time_data)

        # Set axis limits & save figure
        ax.set_xlim(x.min() - radius, x.max() + radius)
        ax.set_ylim(y.min() - radius, y.max() + radius)

        frame_filename = os.path.join(pattern_folder, 'plots', f"{dump_file}_frame_{t}.png")
        fig.savefig(frame_filename, dpi=300, bbox_inches='tight', pad_inches=0)
        frames.append(frame_filename)
        plt.close(fig)


    # Write holes to JSON
    output_file = os.path.join(pattern_folder, 'plots', f"{dump_file}_holes.json")
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)

    # Make a GIF
    gif_filename = os.path.join(pattern_folder, 'plots', f"{dump_file}_area_dist.gif")
    with imageio.get_writer(gif_filename, mode='I', fps=10, loop=0) as writer:
        for filename in frames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Clean up frames
    for filename in frames:
        os.remove(filename)

def main():
    parser = argparse.ArgumentParser(description='Processing dump file name')
    parser.add_argument('folder', type=str, help='Simulation folder name')
    parser.add_argument('trj_file', type=str, help='Dump file name')
    parser.add_argument('--mol_r', nargs="+", type=int, help='Select yarns (mol id)')
    parser.add_argument('--roi', nargs="+", type=float, help='ROI [xmin xmax ymin ymax]')
    parser.add_argument('--yrange', nargs="+", type=float, help='Range for y plots')

    args = parser.parse_args()
    pattern_folder = args.folder #os.path.join('output', 'simulations', args.folder)
    dump_file = args.trj_file
    mol_r = args.mol_r
    xmin, xmax, ymin, ymax = args.roi

    roi_area = (xmax - xmin) * (ymax - ymin)
    roi_width = (xmax - xmin)  # if you ever need it

    # Read dump
    force_data, id_data, mol_data, atom_type, pe_data, coord_atoms = read_lammps_dump(
        os.path.join(pattern_folder, dump_file + '.lammpstrj'),
        xmin, xmax, ymin, ymax
    )

    # We'll build the final list for yarn-length JSON
    detailed_length_data = []

    # Sort each frame by ID (user says consecutive ID means consecutive chain)
    for k in range(len(pe_data)):
        id_flat = id_data[k]
        idx_sort = np.argsort(id_flat)

        # Reorder all arrays by ID
        id_data[k]       = id_data[k][idx_sort]
        mol_data[k]      = mol_data[k][idx_sort]
        atom_type[k]     = atom_type[k][idx_sort]
        pe_data[k]       = pe_data[k][idx_sort]
        force_data[k]    = force_data[k][idx_sort]
        coord_atoms[k]   = coord_atoms[k][idx_sort]

    # Compute yarn lengths at each timestep
    for k in range(len(pe_data)):
        coords_k = coord_atoms[k]
        mol_k    = mol_data[k]
        type_k   = atom_type[k]

        # Group atoms by their 'mol'
        yarn_dict = {}
        for i_atom, m in enumerate(mol_k):
            if m not in yarn_dict:
                yarn_dict[m] = {
                    'coords': [],
                    'types': []
                }
            yarn_dict[m]['coords'].append(coords_k[i_atom])
            yarn_dict[m]['types'].append(type_k[i_atom])

        # For each yarn, sum distances
        yarn_info_list = []
        for mol_val, dct in yarn_dict.items():
            c_list = dct['coords']
            # Summation of adjacent distances
            yarn_length = 0.0
            for i_c in range(len(c_list) - 1):
                p1 = c_list[i_c]
                p2 = c_list[i_c + 1]
                yarn_length += np.linalg.norm(np.array(p2) - np.array(p1))

            # If all atoms have the same type, just take the first
            first_type = dct['types'][0]

            yarn_info_list.append({
                'mol': int(mol_val),
                'type': int(first_type),
                'length': yarn_length
            })

        # Append to global structure
        detailed_length_data.append({
            'time': k,
            'yarns': yarn_info_list
        })

    # Write the yarn-length data
    length_file = os.path.join(pattern_folder, 'plots', f"{dump_file}_yarn_lengths.json")
    with open(length_file, 'w') as f:
        json.dump(detailed_length_data, f, indent=2)

    # Finally, compute hole areas + plot + GIF
    save_holes(coord_atoms, pattern_folder, dump_file, roi_area)

if __name__ == "__main__":
    main()