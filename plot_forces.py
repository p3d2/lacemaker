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
from scipy.fft import fft
from scipy.interpolate import interp1d, RectBivariateSpline, RBFInterpolator
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

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

def save_plot_1(data1, data2, data3, mol_range, N, bin, yax_min, y_axmax):

    ylen = []

    lq = []
    med = []
    uq = []
    maxq = []

    for t in range(len(data1)):
        for k in range(len(data1[0])):
            if data3[t][k][0] in mol_range:


                f_l0 = np.sum(np.array(data1[0][k]).ravel())
                f_l = np.sum(np.array(data1[t][k]).ravel())

                ylen.append((data3[t][k][0],(f_l-f_l0)/f_l0))

                f_p = np.array(data2[t][k]).ravel()
                f_p0 = np.percentile(np.array(data2[0][k]).ravel(),50)

                lq.append((data3[t][k][0],np.percentile(f_p,25)-f_p0))
                med.append((data3[t][k][0],np.percentile(f_p,50)-f_p0))
                uq.append((data3[t][k][0],np.percentile(f_p,75)-f_p0))
                maxq.append((data3[t][k][0],np.percentile(f_p,99)-f_p0))

    labels = sorted(set(x for x, _ in lq))

    n_labels = len(labels)
    fig, axes = plt.subplots(nrows=1, ncols=n_labels, sharex=False, sharey=True, figsize=(5 * n_labels, 5))

    # Check if axes is an array or a single AxesSubplot object (in case there's only one label)
    if n_labels == 1:
        axes = [axes]  # make it a list for consistent handling below
    # Iterate over each label to create its subplot
    for ax, label in zip(axes, labels):
        # Extract data for the current label
        
        x_ylen, y_ylen = zip(*[(x, y) for x, y in ylen if x == label])

        x_lq, y_lq = zip(*[(x, y) for x, y in lq if x == label])
        x_med, y_med = zip(*[(x, y) for x, y in med if x == label])
        x_uq, y_uq = zip(*[(x, y) for x, y in uq if x == label])
        x_max, y_max = zip(*[(x, y) for x, y in maxq if x == label])

        t = np.arange(1,len(y_ylen)+1)
        # Plot the data in the current axis
        ax.plot(t, y_ylen, '--', color='blue', zorder=2)
        #ax.plot(t, y_max, color='red', zorder=2)
        ax.plot(t, y_med, color='black', zorder=2)
        ax.fill_between(t, y_lq, y_uq, alpha=0.5, zorder=1, color='gray')

        # Optional: set a title or a legend for each subplot
        ax.set_title(f"MolID: {label}")
        ax.set_xlim(1.0, np.max(t))
        ax.set_ylim(yax_min, yax_max)
        
        if label == labels[0]:
            ax.set_ylabel('PE-PE(0)')

    fig.text(0.5, 0.0,'Iterations', va='top')  # Set x-axis label for each subplot or you can move this to only the bottom subplot        
    
    plt.subplots_adjust(wspace=0.0)
    # Save the entire figure
    plt.savefig(os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_combined.png'), format='png', dpi=150, bbox_inches='tight', pad_inches=0)

    # Close the figure to free up memory
    plt.close(fig)

def fill_between_3d(ax, X, Y, Z1, Z2, where=None, color='grey', alpha=0.5):
    verts = []
    for i in range(len(X)-1):
        if where is None or (where[i] and where[i+1]):  # Check condition at each step
            verts.append([
                (X[i], Y[i], Z1[i]), (X[i], Y[i], Z2[i]),
                (X[i+1], Y[i+1], Z2[i+1]), (X[i+1], Y[i+1], Z1[i+1])
            ])
    
    # Create the PolyCollection object with the specified vertices
    poly = Poly3DCollection(verts, facecolors=color, alpha=alpha)
    ax.add_collection3d(poly)

def save_plot_2(data1, data2, data3, mol_range, N, bin, yax_min, y_axmax):

    Nbins = 100
    bin_edges = np.linspace(0, 0.675, Nbins-1)  # 99 bins between 0 and 0.675
    bin_edges = np.append(bin_edges, np.inf)

    n_labels = len(mol_range)
    #fig, axes = plt.subplots(nrows=1, ncols=n_labels, sharex=False, sharey=True, figsize=(5 * n_labels, 5))
    fig = plt.figure(figsize=(5 * n_labels, 10))
    all_x = []
    all_y = []
    all_t = []
    
    ylen = []

    lq = []
    med = []
    uq = []
    maxq = []

    pe_data = []
    e_vec = np.linspace(0,0.675,1000)
    for t in range(len(data1)):
        for k in range(len(data1[0])):
            if data3[t][k][0] in mol_range:
                f_l0 = np.sum(np.array(data1[0][k]).ravel())
                f_l = np.sum(np.array(data1[t][k]).ravel())
                f_p = np.array(data2[t][k]).ravel()
                f_p = f_p[f_p > 0]

                #f_linear = interp1d(f_l, f_p)
                #xnew = np.linspace(f_l[0], f_l[-1], int(N+1))
                #ynew = f_linear(xnew)

                #all_x.append((xnew-xnew[0])/(xnew[-1]-xnew[0]))
                #all_y.append(ynew + 0.01*t)
                ylen.append((data3[t][k][0],(np.mean(np.array(data1[t][k]).ravel())-0.25)/0.25))
                lq.append((data3[t][k][0],np.percentile(f_p,25)))
                med.append((data3[t][k][0],np.percentile(f_p,50)))
                uq.append((data3[t][k][0],np.percentile(f_p,75)))

                #hist, edges = np.histogram(np.array(data2[t][k]).ravel(), bins=bin_edges)
                #pe_hist.append((data3[t][k][0], hist))

                #density = gaussian_kde(np.array(data2[t][k]).ravel())
                #density.covariance_factor = lambda: .1  # Smaller bandwidth makes the curve smoother
                #density._compute_covariance()

                density = gaussian_kde(f_p)
                density.covariance_factor = lambda: .25  # Smaller bandwidth makes the curve smoother
                density._compute_covariance()
                
                pe_data.append((data3[t][k][0].T, t+1, density(e_vec)))
            
                #mol_index = mol_range.index(data3[t][k][0])
                #for q in range(len(f_p)):
                #    axes[mol_index].plot(t, f_p[q], 'ko', alpha=0.05, markeredgecolor='none', markersize=2)

    #all_x_flat = np.array(all_x).ravel()
    #all_y_flat = np.array(all_y).ravel()
    #H, xedges, yedges = np.histogram2d(all_x_flat, all_y_flat, bins=bin)
    #ax.pcolormesh(xedges, yedges, H.T, shading='auto', cmap='turbo')#, vmax=100.0)

    labels = sorted(set(x for x, _ in lq))

    if n_labels == 1:
        axes = [axes]  # make it a list for consistent handling below

    # # Iterate over each label to create its subplot
    time_range = np.linspace(1,len(data1), len(data1))
   
    # for ax, label in zip(axes, labels):
    # #     x_lq, y_lq = zip(*[(x, y) for x, y in lq if x == label])
    #     x_med, y_med = zip(*[(x, y) for x, y in med if x == label])
    #     x_uq, y_uq = zip(*[(x, y) for x, y in uq if x == label])

    #     x_ylen, y_ylen = zip(*[(x, y) for x, y in ylen if x == label])

    #     subsets = [(t, pes) for l, t, pes in pe_data if l == label]

    #     t_pe = np.concatenate([[t] * len(pes) for t, pes in subsets])  # Time values repeated
    #     z_pe = np.concatenate([pes for _, pes in subsets])
    #     data = np.vstack([t_pe, z_pe])
    #     kde = gaussian_kde(data, bw_method=[0.5, 0.1])
    #     time_grid, energy_grid = np.meshgrid(time_range, pe_range)
    #     kde_values = kde(np.vstack([time_grid.ravel(), energy_grid.ravel()])).reshape(time_grid.shape)

    #     ax.pcolormesh(time_grid, energy_grid, kde_values, shading='nearest', cmap='turbo')
    #     ax.plot(np.linspace(1,len(data1)+1,len(data1)), y_med, 'w--', alpha=0.1)
    # #     ax.plot(np.arange(len(data1)), y_med, color='black', zorder=2)  # Plotting the median line
    # #     ax.fill_between(np.arange(len(data1)), y_lq, y_uq, label = str(x_lq), alpha=0.5, zorder=1) 
    #     ax.set_title(f"MolID: {label}")
    #     #ax.set_xlim(1.0, np.max(t))
    #     #ax.set_ylim(yax_min, yax_max)
        
    #     if label == labels[0]:
    #         ax.set_ylabel('PE')

    ###############################################################################################################
    # Calculate KDE for each time and plot
    new_t = np.linspace(time_range.min(), time_range.max(), 100)
    interpolated_kdes = np.zeros((len(new_t), len(e_vec)))
    gs = gridspec.GridSpec(2, n_labels, figure=fig)

    #for ax, label in zip(axes, labels):
    for i in range(n_labels):

        _, _, kde_values = zip(*[(x, y, z) for x, y, z in pe_data if x == labels[i]])
        _, y_med = zip(*[(x, y) for x, y in med if x == labels[i]])
        _, y_uq = zip(*[(x, y) for x, y in uq if x == labels[i]])
        _, y_ylen = zip(*[(x, y) for x, y in ylen if x == labels[i]])

        kdes = np.array(kde_values)
        #spline = RectBivariateSpline(time_range, e_vec, kdes, kx=1, ky=1)
        #kde_interp = spline(new_t, e_vec).T
        
        #rbf = RBFInterpolator(np.array(np.meshgrid(time_range, e_vec)).reshape(2, -1).T, kdes.ravel(), kernel='gaussian', epsilon=3.0)
        #kde_interp = rbf(np.array(np.meshgrid(new_t, e_vec)).reshape(2, -1).T).reshape(len(new_t), len(e_vec))

        #for i, kde_values_at_pe in enumerate(np.array(kde_values)):  # Transpose to work along potential energy dimension
        #    interp_func = interp1d(new_t, kde_values_at_pe, kind='cubic')  # Cubic interpolation
        #    interpolated_kdes[i] = interp_func(e_vec)
        
        #ax.remove()  # Remove the 2D axis
        #ax = fig.add_subplot(1, n_labels, np.where(axes == ax)[0][0] + 1, projection='3d')

        ax = fig.add_subplot(gs[0,i], projection='3d')

        ax.patch.set_facecolor('none')
        for k in range(len(time_range)):
            xl = np.full(len(e_vec), time_range[k])
            yl = e_vec
            zl = kdes[k]/np.amax(kdes[k])
            ax.plot(xl, yl, zl, color='k', linewidth=0.5)
            fill_between_3d(ax, xl, yl, zl*0, zl, color=get_color(y_med[k]), alpha=0.5)
        ax.view_init(86, -4, 86)
        ax.tick_params(axis='z', labelleft=False)

        xmin, xmax = [1.0, max(time_range)]
        ymin, ymax = [0.0, 0.675]
        zmin, zmax = [0.0, 1.0]
        xb = np.array([xmin, xmax, xmax, xmin, xmin])
        yb = np.array([ymin, ymin, ymax, ymax, ymin]) 
        zb = np.array([zmin, zmin, zmin, zmin, zmax])
        edges_kw = dict(color='0.0', linewidth=0.5, linestyle='--', zorder=1e3)
                
        # for j in range(len(xb)-1):
        #     ax.plot([xb[j], xb[j+1]], [yb[j], yb[j+1]], [zmax, zmax], **edges_kw)
        #     ax.plot([xb[j], xb[j+1]], [yb[j], yb[j+1]], [zmin, zmin], **edges_kw)
        #     ax.plot([xb[j], xb[j]], [yb[j], yb[j]], [zmin, zmax], **edges_kw)

        # ax.plot([xb[-2], xb[-2]], [yb[-2], yb[-2]], [zmin, zmax], **edges_kw)
        # Bottom face
        ax.plot_surface(np.array([[xmin, xmin], [xmax, xmax]]),
                np.array([[ymin, ymax], [ymin, ymax]]),
                np.array([[zmin, zmin], [zmin, zmin]]),
                color='0.0', alpha=0.1)

        ax.plot_surface(np.array([[xmin, xmax], [xmin, xmax]]),
                np.array([[ymax, ymax], [ymax, ymax]]),
                np.array([[zmin, zmin], [zmax, zmax]]),
                color='0.0', alpha=0.5)

        ax.plot_surface(np.array([[xmax, xmax], [xmax, xmax]]),
                np.array([[ymin, ymin], [ymax, ymax]]),
                np.array([[zmin, zmax], [zmin, zmax]]),
                color='0.0', alpha=0.3)

        # Plot
        #X, Y = np.meshgrid(time_range, e_vec)
        #kde_interp = bisplev(X, Y, tck)
        #ax.pcolormesh(X, Y, kdes.T, shading='nearest', cmap='turbo')
        
        #ax.plot(time_range, y_med, 'w', alpha=0.5, linewidth=2, zorder=2)
        #x_positions = time_range[1:]-.5
        #for xpos in x_positions:
        #    ax.axvline(x=xpos, color='k', linewidth=1, zorder=1)  # Customize the line style as needed
        ax.set_title(f"MolID: {labels[i]}")
        ax.set_xlabel('')  # Label only the last subplot
        #ax.invert_xaxis()  # Optional: Invert y-axis if you prefer that format
        # Turn off the grid
        ax.grid(False)
        ax.set_proj_type('ortho')
        ax.autoscale(False)
        # Set the limits of the plot
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_zlim(zmin, zmax)

        ax.set_box_aspect([1, 1, 1], zoom=1.25)

        ax = fig.add_subplot(gs[1,i])
        norm = mcolors.Normalize(vmin=0, vmax=0.675)
        ax.plot(time_range, y_ylen, 'k')
        ax.plot((xmin-1, xmax+1), (0,0), linestyle = '--', linewidth=0.5, color='gray')

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(yax_min, yax_max)

        if i==0: ax.set_ylabel('$\Delta L$')
        else:
            ax.set_ylabel('')
            ax.set_yticks([], [])
    # Label the shared y-axis
    #axes[0].set_ylabel('Potential Energy')
   
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.0, hspace=0.0)
    plt.savefig(os.path.join('output', 'simulations', pattern_folder, 'plots', dump_file + '_combined.png'), format='png', dpi=150, bbox_inches='tight', pad_inches=0)
    plt.close('all')

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
yax_min, yax_max = args.yrange

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

#if len(mol_r) == 1:
save_plot_2(len_data, pe_sort, id_new, mol_r, 1e4, 1000, yax_min, yax_max)

#if len(mol_r) > 1:
#    if len(mol_r) == 2:
#        mol_r = range(mol_r[0], mol_r[1]+1)
#
#    save_plot_1(len_data, pe_sort, id_new, mol_r, 1e4, 1000, yax_min, yax_max)