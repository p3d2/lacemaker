# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 17:11:00 2024

@author: silvap1
"""

import argparse, json, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.fft import fft
from scipy.interpolate import interp1d
from scipy.signal import find_peaks

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

def save_plot_1(data1, data2, data3, mol_range, N, bin):

    all_x = []
    all_y = []

    lq = []
    med = []
    uq = []
    maxq = []

    for t in range(len(data1)):
        for k in range(len(data1[0])):
            if data3[t][k][0] in mol_range:

                f_l = np.cumsum(np.array(data1[t][k]).ravel())
                f_p = np.array(data2[t][k]).ravel()

                lq.append((data3[t][k][0],np.percentile(f_p,25)))
                med.append((data3[t][k][0],np.percentile(f_p,50)))
                uq.append((data3[t][k][0],np.percentile(f_p,75)))
                maxq.append((data3[t][k][0],np.percentile(f_p,99)))

    labels = sorted(set(x for x, _ in lq))

    n_labels = len(labels)
    fig, axes = plt.subplots(nrows=1, ncols=n_labels, sharex=False, sharey=True, figsize=(5 * n_labels, 5))

    # Check if axes is an array or a single AxesSubplot object (in case there's only one label)
    if n_labels == 1:
        axes = [axes]  # make it a list for consistent handling below
    # Iterate over each label to create its subplot
    for ax, label in zip(axes, labels):
        # Extract data for the current label
        x_lq, y_lq = zip(*[(x, y) for x, y in lq if x == label])
        x_med, y_med = zip(*[(x, y) for x, y in med if x == label])
        x_uq, y_uq = zip(*[(x, y) for x, y in uq if x == label])
        x_max, y_max = zip(*[(x, y) for x, y in maxq if x == label])

        # Plot the data in the current axis
        ax.plot(np.arange(len(data1))+1, y_max, color='red', zorder=2)
        ax.plot(np.arange(len(data1))+1, y_med, color='orange', zorder=2)
        ax.fill_between(np.arange(len(data1))+1, y_lq, y_uq, alpha=0.5, zorder=1)

        # Optional: set a title or a legend for each subplot
        ax.set_title(f"MolID: {label}")
        ax.set_xlim(0, len(data1))

        if label == labels[0]:
            ax.set_ylabel('Potential Energy')

    fig.text(0.5, 0.0,'Position (normalized)', va='top')  # Set x-axis label for each subplot or you can move this to only the bottom subplot        
    
    # Save the entire figure
    plt.savefig(os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_combined.png'), format='png', dpi=600, bbox_inches='tight', pad_inches=0)

    # Close the figure to free up memory
    plt.close(fig)

def save_plot_2(data1, data2, data3, mol_range, N):

    fig, ax = plt.subplots()
    all_x = []
    all_y = []
    all_t = []
    
    lq = []
    med = []
    uq = []

    for t in range(len(data1)):
        for k in range(len(data1[0])):
            if data3[t][k][0] in mol_range:

                f_l = np.cumsum(np.array(data1[t][k]).ravel())
                f_p = np.array(data2[t][k]).ravel()

                f_linear = interp1d(f_l, f_p)
                xnew = np.linspace(f_l[0], f_l[-1], int(N+1))
                ynew = f_linear(xnew)

                all_x.append((xnew-xnew[0])/(xnew[-1]-xnew[0]))
                all_y.append(ynew + 0.01*t)

                lq.append((data3[t][k][0],np.percentile(f_p,25)))
                med.append((data3[t][k][0],np.percentile(f_p,50)))
                uq.append((data3[t][k][0],np.percentile(f_p,75)))

    all_x_flat = np.array(all_x).ravel()
    all_y_flat = np.array(all_y).ravel()

    #H, xedges, yedges = np.histogram2d(all_x_flat, all_y_flat, bins=bin)
    #ax.pcolormesh(xedges, yedges, H.T, shading='auto', cmap='turbo')#, vmax=100.0)

    labels = set(x for x, _ in lq) 
    for label in labels:
        x_lq, y_lq = zip(*[(x, y) for x, y in lq if x == label])
        x_med, y_med = zip(*[(x, y) for x, y in med if x == label])
        x_uq, y_uq = zip(*[(x, y) for x, y in uq if x == label])

        ax.plot(np.arange(len(data1)), y_med, color='black', zorder=2)  # Plotting the median line
        ax.fill_between(np.arange(len(data1)), y_lq, y_uq, label = str(x_lq), alpha=0.5, zorder=1) 

    plt.xlabel('Time (iterations)')
    plt.ylabel('PE')
    plt.savefig(os.path.join(dump_file + '_' + '_'.join(map(str, mol_range)) + '.png'), format='png', dpi=300, bbox_inches='tight',pad_inches=0)
    plt.close('all')

parser = argparse.ArgumentParser(description='Processing dump file name')
parser.add_argument('folder', type=str, help='The name of the dump file to be processed')
parser.add_argument('trj_file', type=str, help='The name of the dump file to be processed')
parser.add_argument('--mol_r', nargs="+", type=int, help='Select the yarns (mol id)')
parser.add_argument('--roi', nargs="+", type=float, help='Roi to consider PE')

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

if len(mol_r) == 1:
    save_plot_2(len_data, pe_sort, id_new, mol_r, 1e4)

if len(mol_r) > 1:
    if len(mol_r) == 2:
        mol_r = range(mol_r[0], mol_r[1]+1)
    save_plot_1(len_data, pe_sort, id_new, mol_r, 1e4, 1000)