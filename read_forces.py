# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 17:11:00 2024

@author: silvap1
"""

import argparse, json, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

fg_color = 'white'
bg_color = 'none'

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']
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

def read_lammps_dump(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    forces = []
    pe_atom = []
    pos_atom = []
    num_atoms = 0
    coord_atoms = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if 'ITEM: NUMBER OF ATOMS' in line:
            num_atoms = int(lines[i + 1])
        elif 'ITEM: ATOMS' in line:
            i += 1  # Start reading atom data from the next line
            timestep_data = []
            pe = []
            pos = []
            coord = []
            for _ in range(num_atoms):
                if i >= len(lines):
                    break
                parts = lines[i].split()
                if len(parts) >= 3:  # Ensure there are enough parts in the line
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    fx, fy, fz = float(parts[5]), float(parts[6]), float(parts[7])
                    force_magnitude = np.sqrt(fx**2 + fy**2 + fz**2)
                    timestep_data.append(force_magnitude)
                    pe.append(float(parts[8]))
                    pos.append(float(parts[0]))
                    coord.append([x, y, z])
                i += 1
            forces.append(timestep_data)
            pe_atom.append(pe)
            pos_atom.append(pos)
            coord_atoms.append(coord)
            continue
        i += 1

    return (forces, np.array(pe_atom), np.array(pos_atom), coord_atoms)

# Replace 'your_dump_file.dump' with the path to your LAMMPS dump file
dump_file = 'output_contract.lammpstrj'
(force_data, pe_data, pos, coord_atoms) = read_lammps_dump(dump_file)

for k in range(len(pe_data)):
    pe = pe_data[k]
    id_flat = pos[k]
    c_flat = np.array(coord_atoms[k])
    # Get the sort indices based on the flattened ID matrix
    sort_indices = np.argsort(id_flat)

    # Use the indices to sort both arrays
    pe_data[k] = pe[sort_indices]
    pos[k] = id_flat[sort_indices]
    coord_atoms[k] = c_flat[sort_indices]

len_data = []
pe_yarn = []
for k in range(len(pe_data)):
    len_yarn = 0
    for l in range(6125,7260):
        p1 = coord_atoms[k][l]
        p2 = coord_atoms[k][l+1]
        len_yarn += np.linalg.norm(p2-p1)
            
    len_data.append(len_yarn)
    pe_yarn.append(sum(pe_data[k][(6125+40):(7260-40)]))
    

fig, ax = plt.subplots()
ax.plot(len_data/len_data[0], pe_yarn,'-o')
ax.set_xlabel('L/L_0')
ax.set_ylabel('Potential Energy')
#plt.savefig(r'C:\Users\silvap1\Google Drive\Scripts\lammps_lace\images\plotmol14.png', format='png', dpi=300, bbox_inches='tight',pad_inches=0)
#plt.close('all')

# For a surface plot, the time and position values need to form a grid
X, Y = np.meshgrid(np.arange(len(force_data[0])), np.arange(0,len(force_data)))

# Interpolating force values onto the grid
# This requires the force_values array to be the same shape as X and Y
# If this is not the case, you will need to interpolate your data
#Z = np.array(force_data[10:])  # Replace this with actual interpolated data
Z = pe_data
#Z_clipped = np.clip(Z[:,3149:3738], 0, 1.0)
Z_clipped = np.clip(Z[:,(6125+40):(7260-40)], 0, 1.0)

# Contour plot
fig, ax = plt.subplots()

ax.contourf(Z_clipped, cmap='turbo', levels=100, vmin=0.0, vmax=0.4)  # Set color limits with vmin and vmax
plt.title('Potential Energy')
plt.xlabel('Position')
plt.ylabel('Time')
#plt.savefig(r'C:\Users\silvap1\Google Drive\Scripts\lammps_lace\images\surf.png', format='png', dpi=300, bbox_inches='tight',pad_inches=0)
#plt.close('all')