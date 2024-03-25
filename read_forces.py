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

fg_color = 'white'
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

def read_lammps_dump(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    forces = []
    pe_atom = []
    pos_atom = []
    num_atoms = 0
    coord_atoms = []
    id_atom = []
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
            id = []
            for _ in range(num_atoms):
                if i >= len(lines):
                    break
                parts = lines[i].split()
                if len(parts) >= 3:  # Ensure there are enough parts in the line
                    pos.append(float(parts[0]))
                    id.append(int(parts[1]))
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    fx, fy, fz = float(parts[6]), float(parts[7]), float(parts[8])
                    force_magnitude = np.sqrt(fx**2 + fy**2 + fz**2)
                    timestep_data.append(force_magnitude)
                    pe.append(float(parts[9]))
                    coord.append([x, y, z])
                i += 1
            pos_atom.append(pos)
            id_atom.append(id)
            forces.append(timestep_data)
            pe_atom.append(pe)
            coord_atoms.append(coord)
            continue
        i += 1

    return (forces, np.array(id_atom), np.array(pe_atom), np.array(pos_atom), coord_atoms)

parser = argparse.ArgumentParser(description='Processing dump file name')
parser.add_argument('folder', type=str, help='The name of the dump file to be processed')
parser.add_argument('trj_file', type=str, help='The name of the dump file to be processed')

args = parser.parse_args()
pattern_folder = args.folder
dump_file = args.trj_file

(force_data, id_data, pe_data, pos, coord_atoms) = read_lammps_dump(os.path.join('output','simulations', pattern_folder, dump_file + '.lammpstrj'))

for k in range(len(pe_data)):
    pe = pe_data[k]
    id_atom = id_data[k]
    id_flat = pos[k]
    c_flat = np.array(coord_atoms[k])
    # Get the sort indices based on the flattened ID matrix
    sort_indices = np.argsort(id_flat)

    # Use the indices to sort both arrays
    pe_data[k] = pe[sort_indices]
    id_data[k] = id_atom[sort_indices]
    pos[k] = id_flat[sort_indices]
    coord_atoms[k] = c_flat[sort_indices]

len_data = []
pe_sort = []
for k in range(len(pe_data)):
    len_temp = []
    len_yarn =[]
    pe_temp = []
    pe_yarn = []
    for l in range(len(pe_data[k])-1):
        if id_data[k][l] != id_data[k][l+1]:
            len_yarn.append(len_temp)
            len_temp = []
            pe_yarn.append(pe_temp)
            pe_temp = []
        else:
            p1 = coord_atoms[k][l]
            p2 = coord_atoms[k][l+1]
            len_temp.append(np.linalg.norm(p2-p1))
            pe_temp.append(pe_data[k][l])
    len_data.append(len_yarn)
    pe_sort.append(pe_yarn)  


f_l = np.cumsum([item for sublist in len_data[0][20:21] for item in sublist])
f_p = [item for sublist in pe_sort[0][20:21] for item in sublist]

f_linear = interp1d(f_l, f_p)
xnew = np.linspace(f_l[0], f_l[-1], int(10e5+1))
ynew = f_linear(xnew)

fft_result = fft(ynew)
# Get the power spectrum (magnitude of the FFT)
power_spectrum = np.abs(fft_result)
frequencies = np.fft.fftfreq(len(ynew), d=(f_l[-1]-f_l[0])/10e5)

# Contour plot
#for q in range(20,30):
fig, ax = plt.subplots()

#ax.plot(np.abs(ft_flat))
#ax.plot(frequencies[1:], power_spectrum[1:])
ax.plot(xnew, ynew)
ax.plot(f_l, f_p, 'r.')
#ax.contourf([lst[q] for lst in pe_sort], cmap='turbo', levels=100, vmin=0.0, vmax=0.4)  # Set color limits with vmin and vmax
#ax.plot(pe_sort[:][30])  # Set color limits with vmin and vmax
#ax.plot(frequencies, np.abs(fourier_transform))
#plt.xlim(0,5)
#plt.ylim(0,1000.0)
plt.title('Potential Energy')
plt.xlabel('Position')
plt.ylabel('Time')
#plt.savefig(os.path.join('output','simulations', pattern_folder, dump_file + '-' + str(q) + '.png'), format='png', dpi=300, bbox_inches='tight',pad_inches=0)
plt.savefig(os.path.join('output','simulations', pattern_folder, dump_file + '.png'), format='png', dpi=300, bbox_inches='tight',pad_inches=0)

plt.close('all')