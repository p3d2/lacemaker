# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 17:25:02 2024

@author: silvap1
"""

import numpy as np
import argparse, json, os
from scipy.interpolate import splprep, splev
import matplotlib.pyplot as plt

# Function to load unit pattern from json file
def load_data(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    return data

# Calculate spatial translations between nodes
def calc_translations(nodes, path_list):
    path = path_list['path']
    shift_keys = [list(map(int, key.strip("[]").split(", "))) for key in path_list['shifts'].keys()]
    
    spatial_shifts = []
    for i in range(len(path) - 1):
        current_crossing = path[i]
        next_crossing = path[i + 1]
        
        # Calculate the basic spatial shift
        shift_x = nodes[str(next_crossing)][0] - nodes[str(current_crossing)][0]
        shift_y = nodes[str(next_crossing)][1] - nodes[str(current_crossing)][1]
        
        # Check for special translations at the current crossing
        pair_shift = [current_crossing, next_crossing]
        if pair_shift in shift_keys:
            # Apply the special translation
            shift_x += path_list['shifts'][str(pair_shift)][0]
            shift_y += path_list['shifts'][str(pair_shift)][1]

        spatial_shifts.append((shift_x, shift_y))

    return spatial_shifts

def generate_yarns(nodes, trs, u_yarns):
    # Initialize the list of points with the starting node coordinates
    path_t, start_n, start_l, crossing, z0 = u_yarns

    points = []
    translations = []
    pt_x = nodes[str(start_n)][0]
    pt_y = nodes[str(start_n)][1]
    pt_z = crossing * z0
    crossing = -crossing
    points.append((pt_x, pt_y, pt_z))

    for i in range(len(trs)-1):
        if start_l > (len(trs)-1): start_l = 0
        pt_x += trs[start_l][0]
        pt_y += trs[start_l][1]
        pt_z = crossing * z0
        crossing = -crossing
        points.append((pt_x, pt_y, pt_z))
        translations.append((trs[start_l][0], trs[start_l][1]))
        start_l += 1
    
    return points, translations

def extend_points(points, n, vx, vy):
    extended_points = []
    for k in range(n):
        for point in points:
            # Calculate new x and y with the respective shifts
            new_x = point[0] + k * vx
            new_y = point[1] + k * vy
            extended_points.append((new_x, new_y, point[2]))
    return extended_points

def smooth_yarn(points, arc_length=1.0, smoothness=3, num_points=10000):
    # Convert points to numpy array
    points = np.array(points)
    x, y, z = points.T
    
    # Create the spline
    tck, u = splprep([x, y, z], s=smoothness)

    # Evaluate the spline at evenly spaced points
    u_even = np.linspace(0, 1, num_points)
    x_even, y_even, z_even = splev(u_even, tck)

    # Function to compute the euclidean distance between two points
    def euclidean(p1, p2):
        return np.sqrt(np.sum((p1 - p2)**2))

    # Initialize parameters
    result_list = [tuple(splev(0, tck))]  # Start with the first point
    last_point = result_list[-1]

    # Iterate over the evenly spaced spline points
    for i in range(1, num_points):
        current_point = (x_even[i], y_even[i], z_even[i])
        if euclidean(np.array(last_point), np.array(current_point)) >= arc_length:
            result_list.append(current_point)
            last_point = current_point

    return result_list

def replicate_points(points, nmol, vx, vy):
    replicated_points = []
    for point in points:
       # Calculate new x and y with the respective shifts
       new_x = point[0] + vx
       new_y = point[1] + vy
       replicated_points.append((nmol, new_x, new_y, point[2]))
    return replicated_points

def point_in_roi(point, roi_bounds):
    x, y, z = point[1:]
    return (roi_bounds['x_min'] <= x <= roi_bounds['x_max'] and
            roi_bounds['y_min'] <= y <= roi_bounds['y_max'] and
            roi_bounds['z_min'] <= z <= roi_bounds['z_max'])

def filter_points(curve_points, roi_bounds):
    return [point for point in curve_points if point_in_roi(point, roi_bounds)]

def is_close_to_boundary(x, y, x_min, y_min, x_max, y_max, margin):
    if (np.abs(x-x_min) <= margin or np.abs(x-x_max) <= margin): return True
    if (np.abs(y-y_min) <= margin or np.abs(y-y_max) <= margin): return True
    return False

def save_pattern(p0, trs, filename):
    
    colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999', '#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494']
    coords = []
    for (x0, y0), trans in zip(p0, trs):
        x, y = [x0], [y0]  # Initialize with the starting point
        for dx, dy in trans:
            x.append(x[-1] + dx)
            y.append(y[-1] + dy)
        coords.append((x, y))

    for idx, (x, y) in enumerate(coords):
        # Select color in a cyclic manner
        color = colors[idx % len(colors)]
        plt.plot(x, y, marker='o', color=color, lw = 2)
        
    # Setting the x and y axis labels
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.axis('equal')
    
    # Saving the plot to a file
    plt.savefig(filename)

def calc_special(yb, threshold):
    
    special = []
    
    # Create the list with the ends of all yarns (indexes of yarns are in order! otherwise it does not work)
    m_id = []
    c = 0
    for k in range(1, len(yb)): 
        if yb[k, 1] != yb[k-1, 1]:
            m_id.append((c, k-1))
            c = k
    m_id.append((c, len(yb)-1))
    
    # For the 2 ends of a yarn, check what are the closest points with other yarns that are at least closer than the threshold
    for i in range(len(m_id)-1):
        start_i, end_i = m_id[i]
        p0min = yb[start_i, 2:5]
        p0max = yb[end_i, 2:5]
        
        for j in range(len(m_id)):
            if i == j:
                continue
            dmin = threshold
            dmax = threshold
            cmin = None
            cmax = None
            for m in range(m_id[j][0], m_id[j][1]):
                p1 = yb[m, 2:5]
                dist_min = np.linalg.norm(p1 - p0min)
                dist_max = np.linalg.norm(p1 - p0max)
                if dist_min < dmin:
                    dmin = dist_min
                    cmin = int(yb[m,0])
                if dist_max < dmax:
                    dmax = dist_max
                    cmax = int(yb[m,0])
                    
            if cmin is not None: special.append((int(yb[start_i,0]), cmin))
            if cmax is not None: special.append((int(yb[end_i,0]), cmax))
            
    return special

def write_lammps_data(yarns, yarn_id, yarns_b, roi_bounds, dist, mass, units, threshold, filename="yarns.txt"):
    with open(filename, 'w') as file:
        
        # Calculate bonds
        bonds = []
        for k in range(1,len(yarns)):
            if yarn_id[k] == yarn_id[k-1]:
                bonds.append((yarns[k][0], k-1, k))
        
        # Calculate special bonds -> connect extreme particles with other extreme particles
        special_bonds = calc_special(yarns_b, threshold)
        
        # Calculate angles
        angles = []
        for k in range(2,len(yarns)):
            if yarn_id[k] - yarn_id[k-1] - yarn_id[k-2] == 0:
                angles.append((yarns[k][0], k-2, k-1, k))
        
        # Write headers
        num_atoms = len(yarns)
        num_bonds = len(bonds) + len(special_bonds)
        num_angles = len(angles)
            
        file.write("LAMMPS data file via Python script\n\n")
        file.write(f"{num_atoms} atoms\n")
        file.write(f"{num_bonds} bonds\n")
        file.write(f"{num_angles} angles\n\n")
        file.write("3 atom types\n")  # Assuming two types of atoms
        file.write("3 bond types\n")
        file.write("3 angle types\n\n")
        
        # Box bounds based on ROI
        file.write(f"{(roi_bounds['x_min']-dist)*units} {(roi_bounds['x_max']+dist)*units} xlo xhi\n")
        file.write(f"{(roi_bounds['y_min']-dist)*units} {(roi_bounds['y_max']+dist)*units} ylo yhi\n")
        file.write(f"{(roi_bounds['z_min']-dist)*units} {(roi_bounds['z_max']+dist)*units} zlo zhi\n\n")
        
        # Masses section
        file.write("Masses\n\n")
        file.write(f"1 {mass}  # Mass for atom type 1\n")
        file.write(f"2 {mass}  # Mass for atom type 2\n")
        file.write(f"3 {mass}  # Mass for atom type 3\n\n")
        

        # Atoms section
        file.write("Atoms\n\n")
        for k in range(len(yarns)):
            atom_type, x, y, z = yarns[k]
            file.write(f"{k+1} {yarn_id[k]+1} {atom_type} {x*units} {y*units} {z*units}\n")

        # Bonds section
        nb = len(bonds)
        file.write("\nBonds\n\n")
        for k in range(nb):
            file.write(f"{k+1} {bonds[k][0]} {bonds[k][1]+1} {bonds[k][2]+1}\n")  # Bond type 1
        
        for k in range(len(special_bonds)):
            file.write(f"{nb+k+1} {3} {special_bonds[k][0]+1} {special_bonds[k][1]+1}\n")  # Bond type 1
                
        # Angles section
        file.write("\nAngles\n\n")
        for k in range(len(angles)):
            file.write(f"{k+1} {angles[k][0]} {angles[k][1]+1} {angles[k][2]+1} {angles[k][3]+1}\n")

def main():
    
    # Create argument parser
    parser = argparse.ArgumentParser(description="Load unit pattern from a JSON file.")
    # Add argument for the JSON file path
    parser.add_argument('json_file', type=str, help='Path to the JSON file containing nodes and paths data.')
    parser.add_argument('--dist_particles', type=float, default=1.0, help='Distance between particles. Default: 1.0')
    parser.add_argument('--units', type=float, default=1.0, help='Units. Default: 1.0')
    parser.add_argument('--mass', type=float, default=1.0, help='Mass. Default: 1.0')
    parser.add_argument('--threshold', type=float, default=2.5, help='Threshold. Default: 2.5')

    # Parse arguments
    args = parser.parse_args()
    
    # Load data from the specified JSON file
    data = load_data(args.json_file)
    nodes = data['nodes']
    paths = data['paths']
    unit_yarns = data['unit_yarns']
    unit_rep = data['unit_repetion']
    roi_bounds = data['roi_bounds']
    
    dist_particles = args.dist_particles
    units = args.units
    mass = args.mass
    threshold = args.threshold
    filename = args.json_file.split('/')[-1].split('.')[0] + '_' + str(dist_particles) + '_' + str(units) + '_' + str(threshold) + '.png'
    
    # Calculate translations between nodes of each path
    path_translations = []
    for k in range(len(paths)):
        path_translations.append(calc_translations(nodes, paths[k]))
    
    # Create yarns
    yarns = []
    path_trs = []
    for k in range(len(unit_yarns)):
        # Generate unit yarns
        yarn, translations = generate_yarns(nodes, path_translations[unit_yarns[str(k)][0]], unit_yarns[str(k)])
        path_trs.append(translations)

        # Extend unit yarns
        n1, vx1, vy1, n2, vx2, vy2, mol = unit_rep[str(k)]
        yarn_ext = extend_points(yarn, n1, vx1, vy1)
        
        # Smooth and create fixed point to point distances
        yarn_sm = smooth_yarn(yarn_ext, arc_length=dist_particles/units)
        
        # Replicate yarns
        for l in range(n2+1):
            yarns.append(replicate_points(yarn_sm, mol, l*vx2, l*vy2))

    # Save pattern
    p0 = [nodes[str(unit_yarns[str(k)][1])] for k in range(len(unit_yarns))]
    save_pattern(p0, path_trs, os.path.join('output','patterns_data',filename))

    # Crop yarns -> Cut rectangular section
    yarns = [filter_points(yarn, roi_bounds) for yarn in yarns]
    
    # Flatten the list of yarns
    yarns_flat = [item for sublist in yarns for item in sublist]
    
    # Create a tag (yarn_id) for each independent yarn
    # Create a list of particles close to the boundaries so we can glue them
    yarn_id = []    
    yarns_bound = []
    c = 0
    distance = 0
    threshod_bound = 3.0 # this threshold selects all the particles that are close to the boundary to save some computation time later when calculating points closer to ends
    for k in range(len(yarns_flat)):
        if k > 0: distance = np.linalg.norm(np.array(yarns_flat[k]) - np.array(yarns_flat[k-1]))
        if distance > 2*dist_particles:
            c += 1
        yarn_id.append(c)
        
        if is_close_to_boundary(yarns_flat[k][1], yarns_flat[k][2], roi_bounds['x_min'], roi_bounds['y_min'], roi_bounds['x_max'], roi_bounds['y_max'], threshod_bound):
            yarns_bound.append((k, c, yarns_flat[k][1], yarns_flat[k][2], yarns_flat[k][3]))
    yarn_id = np.array(yarn_id)
    yarns_bound = np.array(yarns_bound)

    write_lammps_data(yarns_flat, yarn_id, yarns_bound, roi_bounds, dist_particles, mass, units, threshold, os.path.join('output','lammps_data',filename))
    
if __name__ == "__main__":
    main()