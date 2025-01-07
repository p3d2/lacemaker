# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 17:25:02 2024

@author: silvap1
"""

import numpy as np
import argparse, json, os, math
from scipy.interpolate import splprep, splev
import matplotlib.pyplot as plt

# Function to convert path nodes that have 'l' or 'r' to integer value
def to_int(value):
    # Remove last character if it's 'l' or 'r'
    val = str(value)
    if val[-1] in 'lr':
        return int(val[:-1])
    return int(val)

# Function to check for 'l' or 'r' at the end
def lr(value):
    val = str(value)
    return {'l': -1, 'r': 1}.get(val[-1], 0)

# Function to generate auxiliary points for twisting
def twist_points(tw_val, ang):
    points = []
    twist = np.abs(tw_val)
    
    theta = ang - math.pi/4

    x = -twist
    y = -twist

    for k in range(2*twist+1):
        points.append((x, y))
        if int( ((k-np.sign(tw_val))%4)/2 ) < 1:
            y += 2
        else:
            x += 2

    rotated_points = [(x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta)) for x, y in points]
    return rotated_points

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
        current_crossing = to_int(path[i])
        next_crossing = to_int(path[i + 1])
        
        # Calculate the basic spatial shift
        shift_x = nodes[str(next_crossing)][0] - nodes[str(current_crossing)][0]
        shift_y = nodes[str(next_crossing)][1] - nodes[str(current_crossing)][1]
        
        # Check for special translations at the current crossing
        pair_shift = [current_crossing, next_crossing]
        if pair_shift in shift_keys:
            # Apply the special translation
            shift_x += path_list['shifts'][str(pair_shift)][0]
            shift_y += path_list['shifts'][str(pair_shift)][1]

        # Calculate twist value
        tw_val = lr(path[i + 1]) * nodes[str(next_crossing)][2]

        spatial_shifts.append((shift_x, shift_y, tw_val, 0))

    return spatial_shifts

# Function to generate auxiliary points for twisting
def twist_points(tw_val, ang):
    points = []
    twist = np.abs(tw_val)
    
    theta = ang - math.pi/4

    x = -twist
    y = -twist

    for k in range(2*twist+1):
        points.append((x, y))
        if int( ((k-np.sign(tw_val))%4)/2 ) < 1:
            y += 2
        else:
            x += 2

    rotated_points = [(x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta)) for x, y in points]
    return rotated_points

# Calculate angles for special nodes
def node_angle(trs, paths):
    
    special_node_indices = {}
    seen_nodes = set()

    # Iterate over each path and its content
    for path_index, path_dict in enumerate(paths):
        path = path_dict['path']
        for node_index, node in enumerate(path):
            if isinstance(node, str) and (node.endswith('l') or node.endswith('r')):  # Check for special node
                numeric_key = node[:-1]  # Strip the last character to get the numeric prefix
                if node in seen_nodes:
                    raise ValueError(f"Duplicate node '{node}' found in paths, which is not allowed.")
                seen_nodes.add(node)  # Add node to the set to track that it's been processed
                
                if numeric_key not in special_node_indices:
                    special_node_indices[numeric_key] = []
                special_node_indices[numeric_key].append((path_index, node_index))
    
    # Validate the count of indices per key
    for key, indices in special_node_indices.items():
        if len(indices) != 2:
            raise ValueError(f"Special node '{key}' does not have exactly one 'l' and one 'r' entries, which is required to represent the twist between two yarns.")

    for key, indices in special_node_indices.items():
        (l1, u1), (l2, u2) = indices
        size1 = len(trs[l1])
        size2 = len(trs[l2])
        # Angle in respect to (0,1) vector
        angle1 = math.atan2(
            trs[l1][u1-1][1] + trs[l2][u2-1][1],
            trs[l1][u1-1][0] + trs[l2][u2-1][0]
        )
        angle2 = math.atan2(
            trs[l1][u1 % size1][1] + trs[l2][u2 % size2][1],
            trs[l1][u1 % size1][0] + trs[l2][u2 % size2][0]
        ) # the mod operator ensures cyclic calculation of the 2 vectors
        trs[l1][u1-1] = trs[l1][u1-1][:3] + ((angle1+angle2)/2,)
        trs[l2][u2-1] = trs[l2][u2-1][:3] + ((angle1+angle2)/2,)
        #print((indices, angle1, angle2))

    return trs

def generate_yarns(nodes, trs, u_yarns, rad):
    # Initialize the list of points with the starting node coordinates
    _, start_n, start_l, crossing, z0 = u_yarns
    points = []
    translations = []
    pt_x = nodes[str(start_n)][0]
    pt_y = nodes[str(start_n)][1]
    pt_z = crossing * z0
    crossing = -crossing
    #points.append((pt_x, pt_y, pt_z))

    for i in range(len(trs)):
        if start_l > (len(trs)-1): start_l = 0
        
        dx = trs[start_l][0]
        dy = trs[start_l][1]
        twists = trs[start_l][2]

        # Check if the fist node is a special node. If it is a special % 2 == 1, then there is no cross inversion
        if i == 0 and abs(nodes[str(start_n)][2]) > 0 and nodes[str(start_n)][2] % 2 != 0:
            crossing = -crossing

        pt_x += dx
        pt_y += dy
        pt_z = crossing * z0
        
        if twists != 0:
            angle = trs[start_l][3]
            aux_pts = twist_points(twists, angle)
            for k in range(len(aux_pts)):
                px = pt_x + aux_pts[k][0]*rad
                py = pt_y + aux_pts[k][1]*rad
                pz = crossing * z0 * math.cos(k * math.pi / 2) # the cosine cycles over 1, 0, -1, 0, and so on
                points.append((px, py, pz))        
        else:
            points.append((pt_x, pt_y, pt_z))

        if twists % 2 == 0: crossing = -crossing

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

def smooth_yarn(points, arc_length=1.0, smoothness=0.1, num_points=int(2e5)):
    # Convert points to numpy array
    points = np.array(points)
    x, y, z = points.T
    
    # Create the spline
    tck, _ = splprep([x, y, z], s=smoothness)

    # Evaluate the spline at evenly spaced points
    u_even = np.linspace(0, 1, num_points)
    x_even, y_even, z_even = splev(u_even, tck)

    # Initialize parameters
    result_list = [tuple(splev(0, tck))]  # Start with the first point
    last_point = result_list[-1]

    # Iterate over the evenly spaced spline points
    for i in range(1, num_points):
        current_point = (x_even[i], y_even[i], z_even[i])
        if np.linalg.norm(np.array(last_point) - np.array(current_point)) >= arc_length:
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
        plt.plot(x, y, marker='o', color=color, lw = 2, alpha=0.5)
        
    # Setting the x and y axis labels
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.axis('equal')
    
    # Saving the plot to a file
    plt.savefig(filename)

def calculate_bonds(yarns, yarn_id, yarns_b, threshold, m_max):

    # Calculate bonds
    bonds = []
    for k in range(1,len(yarns)):
        if yarn_id[k] == yarn_id[k-1]:
            bonds.append((yarns[k][0], k-1, k))

    # Calculate special bonds -> connect extreme particles with other extreme particles
    special_bonds = calc_special(yarns_b, threshold, m_max)
    return bonds + special_bonds

def calc_special(yb, threshold, mol):
    
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
    mol += 1
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
                    
            if cmin is not None: 
                special.append((mol, int(yb[start_i,0]), cmin, dmin))
                mol += 1
            if cmax is not None:
                special.append((mol, int(yb[end_i,0]), cmax, dmax))
                mol += 1
            
    return special

def calculate_angles(yarns, yarn_id):
    # Calculate angles
    angles = []
    for k in range(2,len(yarns)):
        if yarn_id[k] == yarn_id[k-1] == yarn_id[k-2]:
            angles.append((yarns[k][0], k-2, k-1, k))
    return angles

def write_lammps_data(yarns, yarn_id, bonds, angles, roi_bounds, dist, mass, units, ks1, ks2, kb, filename="yarns.txt"):
    with open(filename, 'w') as file:
        
        if 1 == 1: # bring origin to 0
            xmin = roi_bounds['x_min']
            ymin = roi_bounds['y_min']

        # Write headers
        num_atoms = len(yarns)
        num_bonds = len(bonds)
        num_angles = len(angles)
        
        num_atom_types = max(yarns, key=lambda x: x[0])[0]
        num_bond_types = max(bonds, key=lambda x: x[0])[0]
        num_angle_types = max(angles, key=lambda x: x[0])[0]

        file.write("LAMMPS data file via Python script\n\n")
        file.write(f"{num_atoms} atoms\n")
        file.write(f"{num_bonds} bonds\n")
        file.write(f"{num_angles} angles\n\n")
        file.write(f"{num_atom_types} atom types\n")  # Assuming two types of atoms
        file.write(f"{num_bond_types} bond types\n")
        file.write(f"{num_angle_types} angle types\n\n")
        
        # Box bounds based on ROI
        file.write(f"{(roi_bounds['x_min']-dist - xmin)*units} {(roi_bounds['x_max']+dist - xmin)*units} xlo xhi\n")
        file.write(f"{(roi_bounds['y_min']-dist - ymin)*units} {(roi_bounds['y_max']+dist - ymin)*units} ylo yhi\n")
        file.write(f"{(roi_bounds['z_min']-dist)*units} {(roi_bounds['z_max']+dist)*units} zlo zhi\n\n")
        
        # Masses section
        file.write("Masses\n\n")
        for k in range(num_atom_types):
            file.write(f"{k+1} {mass}  # Mass for atom type {k+1}\n")
        
        # Atoms section
        file.write("\nAtoms\n\n")
        for k, yarn in enumerate(yarns):
            atom_type, x, y, z = yarn
            file.write(f"{k+1} {yarn_id[k]+1} {atom_type} {(x - xmin)*units} {(y - ymin)*units} {z*units}\n")

        # Bonds section
        file.write("\nBonds\n\n")
        for k, bond in enumerate(bonds):
            file.write(f"{k+1} {bond[0]} {bond[1]+1} {bond[2]+1}\n") 
            
        # Angles section
        file.write("\nAngles\n\n")
        for k, angle in enumerate(angles):
            file.write(f"{k+1} {angle[0]} {angle[1]+1} {angle[2]+1} {angle[3]+1}\n")

        # Bond Coeffs section
        file.write("\nBond Coeffs\n\n")
        nb = 1
        for k, bond in enumerate(bonds):
            if bond[0] == nb:
                if len(bond) < 4: file.write(f"{nb} {ks1} {dist}\n")
                else: file.write(f"{nb} {ks2} {bond[3]*units}\n")
                nb += 1

        # Angle Coeffs section
        file.write("\nAngle Coeffs\n\n")
        na = 1
        for k, angle in enumerate(angles):
            if angle[0] == na:
                file.write(f"{na} {kb} {180.0}\n")
                na += 1
def main():
    
    # Create argument parser
    parser = argparse.ArgumentParser(description="Load unit pattern from a JSON file.")
    # Add argument for the JSON file path
    parser.add_argument('json_file', type=str, help='Path to the JSON file containing nodes and paths data.')
    parser.add_argument('--dist_particles', type=float, default=1.0, help='Distance between particles. Default: 1.0')
    parser.add_argument('--units', type=float, default=1.0, help='Units. Default: 1.0')
    parser.add_argument('--mass', type=float, default=1.0, help='Mass. Default: 1.0')
    parser.add_argument('--threshold', type=float, default=2.5, help='Threshold. Default: 2.5')
    parser.add_argument('--ks1', type=float, default=30.0, help='Yarns stretching constant. Default: 30.0')
    parser.add_argument('--ks2', type=float, default=10.0, help='Special yarns stretching constant. Default: 10.0')
    parser.add_argument('--kb', type=float, default=5.0, help='Yarns bending constant. Default: 50.0')
    parser.add_argument('--twist_distance', type=float, default=0.5, help='Distance between auxiliary points when twisting two yarns. Default: 0.5')
    parser.add_argument('--invert_y', type=float, default=-1.0, help='Invert the y axis. Default: -1.0 (invert)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Load data from the specified JSON file (note json structure changed and this file is still trying to keep the old structure while reading the new format)
    data = load_data(args.json_file)

    new_nodes = data["nodes"]  # new style
    old_nodes = {}
    for node_id_str, node_data in new_nodes.items():
        node_id = str(node_id_str)
        x = node_data["x"]
        y = node_data["y"]
        twist = node_data["twist"]
        old_nodes[node_id] = [x, y, twist]
    data["nodes"] = old_nodes

    # 2) Convert each path's "shifts" from a list of objects to a dict with keys like "[from, to]"
    #    so that data["paths"][i]["shifts"] is consistent with your old usage.
    for path_dict in data["paths"]:
        new_shifts_list = path_dict.get("shifts", [])
        old_shifts_dict = {}
        for shift_obj in new_shifts_list:
            key = f"[{shift_obj['from']}, {shift_obj['to']}]"
            old_shifts_dict[key] = [shift_obj["dx"], shift_obj["dy"]]
        path_dict["shifts"] = old_shifts_dict

    # 3) If your old code used 'unit_rep = data["unit_repetion"]', then reconstruct
    new_unit_yarns = data["unit_yarns"]
    old_unit_yarns = {}
    for yarn_id_str, yarn_info in new_unit_yarns.items():
        # Build the old 5-element list
        arr = [
            yarn_info["path_id"],          # index 0
            yarn_info["starting_node"],    # index 1
            yarn_info["start_path_index"], # index 2
            yarn_info["z_sign"],           # index 3
            yarn_info["z_height"]          # index 4
        ]
        old_unit_yarns[yarn_id_str] = arr
    data["unit_yarns"] = old_unit_yarns

    old_unit_rep = {}
    # We need to grab the "repetitions" from the original new_unit_yarns dict,
    # because we just overwrote data["unit_yarns"] with the old arrays.
    for yarn_id_str, yarn_info in new_unit_yarns.items():
        reps = yarn_info["repetitions"]
        name = yarn_info["name"]  # or any other field your old code expected
        # Example: build a 7-element array [count1, dx1, dy1, count2, dx2, dy2, name]
        # If you have more or fewer repetitions, adapt accordingly.
        if len(reps) >= 2:
            rep_array = [
                reps[0]["count"], reps[0]["dx"], reps[0]["dy"],
                reps[1]["count"], reps[1]["dx"], reps[1]["dy"],
                name
            ]
        elif len(reps) == 1:
            # If there's only one repetition object
            rep_array = [
                reps[0]["count"], reps[0]["dx"], reps[0]["dy"],
                0, 0.0, 0.0, name
            ]
        else:
            # No repetitions? Just fill with zeros
            rep_array = [0, 0.0, 0.0, 0, 0.0, 0.0, name]

        old_unit_rep[yarn_id_str] = rep_array

    data["unit_repetion"] = old_unit_rep

    # 4) For convenience, rename some data keys into your old variables,
    #    so the rest of your code can use them just as before:
    nodes        = data["nodes"]
    paths        = data["paths"]
    unit_yarns   = data["unit_yarns"]
    unit_rep     = data["unit_repetion"]
    roi_bounds   = data["roi_bounds"]

    dist_particles = args.dist_particles
    units = args.units
    mass = args.mass
    threshold = args.threshold
    ks1 = args.ks1
    ks2 = args.ks2
    kb = args.kb
    r_tw = args.twist_distance
    invert = args.invert_y
    filename = args.json_file.split('/')[-1].split('.')[0] + '_' + str(dist_particles) + '_' + str(units) + '_' + str(threshold) + '_' + str(ks1) + '_' + str(kb)

    # Invert y-axis
    if invert < 0:
        for key, values in nodes.items():
            nodes[key] = [values[0], -values[1], values[2]]

        for path in paths:
            updated_path = []
            for item in path["path"]:
                if isinstance(item, str) and 'l' in item:
                    updated_path.append(item.replace('l', 'r'))
                elif isinstance(item, str) and 'r' in item:
                    updated_path.append(item.replace('r', 'l'))
                else:
                    updated_path.append(item)
            path["path"] = updated_path

            for key, values in path["shifts"].items():
                path["shifts"][key] = [values[0], -values[1]]
        
        for key, values in unit_rep.items():
            values[2] = -values[2]  # Invert the 3rd element
            values[5] = -values[5]  # Invert the 6th element
            unit_rep[key] = values  # Assign the modified list back to the dictionary

        roi_bounds['y_min'], roi_bounds['y_max'] = -roi_bounds['y_max'], -roi_bounds['y_min']

    # Calculate translations between nodes of each path
    path_translations = []
    for k in range(len(paths)):
        path_translations.append(calc_translations(nodes, paths[k]))
    
    path_translations = node_angle(path_translations, paths)
    # Create yarns
    yarns = []
    path_trs = []
    mol_max = 0
    for k in range(len(unit_yarns)):
        # Generate unit yarns
        yarn, translations = generate_yarns(nodes, path_translations[unit_yarns[str(k)][0]], unit_yarns[str(k)], r_tw)
        path_trs.append(translations)
        
        # Extend unit yarns
        n1, vx1, vy1, n2, vx2, vy2, mol = unit_rep[str(k)]
        if mol > mol_max: mol_max = mol
        yarn_ext = extend_points(yarn, n1, vx1, vy1)
        
        # Smooth and create fixed point to point distances
        yarn_sm = smooth_yarn(yarn_ext, arc_length=dist_particles/units)
        
        # Replicate yarns
        for l in range(n2+1):
            yarns.append(replicate_points(yarn_sm, mol, l*vx2, l*vy2))

    # Save pattern
    p0 = [nodes[str(unit_yarns[str(k)][1])] for k in range(len(unit_yarns))]
    #save_pattern(p0, path_trs, os.path.join('output','patterns_data',filename + '.png'))

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

    # Calculate bonds (connect consecutive particles in the same yarn and ends of yarns within a threshold with other yarns)
    bonds = calculate_bonds(yarns_flat, yarn_id, yarns_bound, threshold, mol_max)

    # Calculate angles
    angles = calculate_angles(yarns_flat, yarn_id)

    write_lammps_data(yarns_flat, yarn_id, bonds, angles, roi_bounds, dist_particles, mass, units, ks1, ks2, kb, os.path.join('output','lammps_data', filename + '.data'))
    
if __name__ == "__main__":
    main()