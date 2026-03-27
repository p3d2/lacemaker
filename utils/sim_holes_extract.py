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
import cv2

fg_color = 'black'
bg_color = 'none'

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['DejaVu Sans']
mpl.rcParams['font.size'] = 20

MIN_HOLE_AREA = 0.1  # area threshold
MAX_CENTROID_DIST = 2.5  # tweak this distance threshold to match holes from frame to frame
radius = 0.1# yarn radiues

CSTEP = 10

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
    frames2 = []
    method = 'log10'
    max_area = 200.0  # for color scale in plot

    for t in range(len(positions)):
        pos_0 = np.array(positions[t])
        x = pos_0[:, 0]
        y = pos_0[:, 1]

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(8, 8))
        fig2, ax2 = plt.subplots(figsize=(8, 8))
        ax.set_facecolor('white')
        ax2.set_facecolor('white')

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
                sim2exp = (50 / 128) ** 2
                area = ShapelyPolygon(interior_coords).area * sim2exp

                # Collect in "areas" list (all holes, big or small)        
                areas.append(area)

                # Decide the patch color (you can keep the same color scale):
                if method == 'log10':
                    color = cmap(np.log10(area + 1) / np.log10(max_area))
                elif method == 'sqrt':
                    color = cmap(np.sqrt(area / max_area))
                else:
                    color = 'grey'

                # Plot this hole (regardless of its size)
                edge_color = 'red' if area < MIN_HOLE_AREA else 'grey'
                edge_color2 = 'grey'
                hole_patch = Polygon(interior_coords, facecolor=color, edgecolor=edge_color, linewidth=1, alpha=1.0)
                hole_patch2 = Polygon(interior_coords, facecolor=color, edgecolor=edge_color2, linewidth=1, alpha=1.0)
                ax.add_patch(hole_patch)
                ax2.add_patch(hole_patch2)

                # Now do the labeling logic, but ONLY store big holes in holes_this_frame.
                if area >= MIN_HOLE_AREA:
                    centroid = ShapelyPolygon(interior_coords).centroid
                    hole_info = {
                        'area': area,
                        'centroid': (centroid.x, centroid.y),
                        'coords': interior_coords
                    }
                    holes_this_frame.append(hole_info)

        # If you want the merged region (the big grey region) displayed:
        collection = PatchCollection(patches, facecolor='grey', edgecolor='black', alpha=0.5, zorder = -1)
        collection2 = PatchCollection(patches, facecolor='grey', edgecolor='black', alpha=0.5, zorder = -1)
        ax.add_collection(collection)
        ax2.add_collection(collection2)

        # Now that we have all holes >= MIN_HOLE_AREA for this frame:
        # (A) If t == 0, assign new labels to all
        # (B) If t > 0, match holes to existing labels if centroids are close enough
        frame_holes_with_labels = []  # store { 'label': X, 'area': Y }

        if t == 0:
            # First frame: assign new labels to all holes above threshold
            for hole_info in holes_this_frame:
                hole_labels[next_label_id] = hole_info['centroid']
                frame_holes_with_labels.append({
                    'label': next_label_id,
                    'area': hole_info['area'],
                    'centroid': hole_info['centroid']
                })
                next_label_id += 1
        # else:
        #     # Subsequent frames: try to match holes to existing labels
        #     used_labels = set()
        #     for hole_info in holes_this_frame:
        #         cx, cy = hole_info['centroid']
        #         best_label = None
        #         best_dist = float('inf')

        #         for label_id, old_centroid in hole_labels.items():
        #             ox, oy = old_centroid
        #             dist = np.sqrt((cx - ox)**2 + (cy - oy)**2)
        #             if dist < best_dist and dist < MAX_CENTROID_DIST:
        #                 best_label = label_id
        #                 best_dist = dist

        #         if best_label is not None:
        #             # update that label’s centroid
        #             hole_labels[best_label] = (cx, cy)
        #             frame_holes_with_labels.append({
        #               'label': best_label,
        #               'area': hole_info['area'],
        #               'centroid': (cx, cy)
        #           })
        #             used_labels.add(best_label)
        #         else:
        #             # no good match => new label
        #             hole_labels[next_label_id] = (cx, cy)
        #             frame_holes_with_labels.append({
        #                 'label': next_label_id,
        #                 'area': hole_info['area'],
        #                 'centroid': (cx, cy)
        #             })
        #             used_labels.add(next_label_id)
        #             next_label_id += 1
        else:
            # Subsequent frames: one-to-one matching between existing labels and current holes.
            candidate_matches = []
            for idx, hole_info in enumerate(holes_this_frame):
                cx, cy = hole_info['centroid']
                for label_id, old_centroid in hole_labels.items():
                    ox, oy = old_centroid
                    dist = np.sqrt((cx - ox)**2 + (cy - oy)**2)
                    if dist < MAX_CENTROID_DIST:
                        candidate_matches.append((dist, label_id, idx))
            
            # Sort all candidate matches by distance (smallest first)
            candidate_matches.sort(key=lambda x: x[0])
            
            assigned_labels = set()
            assigned_holes = set()
            for dist, label_id, idx in candidate_matches:
                # Assign the label if it hasn't been used and this hole hasn't been matched yet.
                if label_id not in assigned_labels and idx not in assigned_holes:
                    hole_info = holes_this_frame[idx]
                    cx, cy = hole_info['centroid']
                    # Update the centroid for the matched label.
                    hole_labels[label_id] = (cx, cy)
                    frame_holes_with_labels.append({
                        'label': label_id,
                        'area': hole_info['area'],
                        'centroid': (cx, cy)
                    })
                    assigned_labels.add(label_id)
                    assigned_holes.add(idx)

            # (Optional) remove labels that didn't appear at all this frame
            # hole_labels = {
            #     lid: c for lid, c in hole_labels.items() if lid in used_labels
            # }

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

        ax2.set_xlim(x.min() - radius, x.max() + radius)
        ax2.set_ylim(y.min() - radius, y.max() + radius)
        ax2.set_xlabel('')
        ax2.set_ylabel('')
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_xticks([])
        ax2.set_yticks([])

        if t % CSTEP == 0:
            frame_filename2 = os.path.join(pattern_folder, 'plots', f"{dump_file}_frame2_{t}.png")
            fig2.savefig(frame_filename2, dpi=96, bbox_inches='tight', pad_inches=0.01)
            frames2.append(frame_filename2)
            plt.close(fig2)


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

    png_filename = os.path.join(pattern_folder, 'plots', f"{dump_file}_area_dist.png")
    images = []
    for k in range(len(frames2)):
        image = imageio.imread(frames2[k])
        cv2.putText(image, f"Step {CSTEP*k}", (10, 560), cv2.FONT_HERSHEY_SIMPLEX, 1.5, (0, 0, 0), 15)
        cv2.putText(image, f"Step {CSTEP*k}", (10, 560), cv2.FONT_HERSHEY_SIMPLEX, 1.5, (255, 255, 255), 5)
        images.append(image)
    # Combine images side by side (horizontally)
    thumbnail = np.hstack(images)
    
    #ROWS
    #row1 = np.hstack(images[:2]) # First row: Stack first two images horizontally
    # Second row: The third image needs to be aligned with the first column
    #empty_space = np.ones_like(images[0]) * 255  # Create a blank white image as a placeholder
    #row2 = np.hstack([images[2], empty_space])  # Align third image to the first column
    # Stack both rows vertically
    #thumbnail = np.vstack([row1, row2])
    
    # Save the combined image
    thumbnail_rgb = cv2.cvtColor(thumbnail, cv2.COLOR_BGR2RGB)
    jpeg_filename = png_filename.replace('.png', '.jpg')
    cv2.imwrite(jpeg_filename, thumbnail_rgb, [cv2.IMWRITE_JPEG_QUALITY, 80]) 
    imageio.imwrite(png_filename, cv2.cvtColor(thumbnail_rgb, cv2.COLOR_RGB2BGR))

    # Clean up frames
    for filename in frames:
        os.remove(filename)

    for filename in frames2:
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