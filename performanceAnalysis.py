import json
import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def analyze_files(holes_file, yarn_file, pattern_folder, dump_file):

    # 1) LOAD holes.json
    with open(holes_file, 'r') as f:
        holes_data = json.load(f)  # [ { "time": t, "areas_label": [...], ... }, ... ]

    # Build a structure: label_areas[label][t] = area
    label_areas = {}
    for frame in holes_data:
        t = frame["time"]
        for hole_dict in frame["areas_label"]:
            label = hole_dict["label"]
            area  = hole_dict["area"]
            if label not in label_areas:
                label_areas[label] = {}
            label_areas[label][t] = area

    label_rel_changes = {}  # { label: { "A0": value, "rel_changes": {t: rel_change}, "max_rel": value, "mean_rel": value } }
    for label, areas_over_time in label_areas.items():
        if 0 in areas_over_time and areas_over_time[0] > 0:
            A0 = areas_over_time[0]
            rel_changes = { t: abs(area - A0) / A0 for t, area in areas_over_time.items() }
            label_rel_changes[label] = {
                "A0": A0,
                "rel_changes": rel_changes,
                "max_rel": max(rel_changes.values()),
                "mean_rel": np.mean(list(rel_changes.values()))
            }

    # Find which label has the largest area overall
    best_label = None
    best_area = -1.0
    best_time = None

    for label, times_dict in label_areas.items():
        for t, area in times_dict.items():
            if area > best_area:
                best_label = label
                best_area  = area
                best_time  = t

    # For that label, get the area at t=0 (if it exists)
    area_at_t0 = label_areas[best_label].get(0, None)

    # 2) LOAD yarn_lengths.json
    with open(yarn_file, 'r') as f:
        yarn_data = json.load(f)  # [ { "time": t, "yarns": [{ "mol": X, "length": L }, ...] }, ... ]

    # Build a structure: yarn_lengths[mol][t] = length
    yarn_lengths = {}
    for frame in yarn_data:
        t = frame["time"]
        for yinfo in frame["yarns"]:
            mol_id = yinfo["mol"]
            length = yinfo["length"]
            if mol_id not in yarn_lengths:
                yarn_lengths[mol_id] = {}
            yarn_lengths[mol_id][t] = length

    # Find the minimum (l(t) - l(0)) / l(0)
    best_ratio = None
    best_mol   = None
    best_time_yarn = None

    for mol, times_dict in yarn_lengths.items():
        if 0 not in times_dict:
            # Skip if no time=0
            continue
        l0 = times_dict[0]
        if l0 == 0:
            continue

        # For each time, compute ratio
        for t, length_at_t in times_dict.items():
            ratio = (length_at_t - l0) / l0
            if (best_ratio is None) or (ratio < best_ratio):
                best_ratio = ratio
                best_mol   = mol
                best_time_yarn = t

    # Store the analysis in a dictionary
    analysis_results = {
        "largest_hole_label": best_label,
        "area_at_t0": area_at_t0,
        "largest_area": best_area,
        "time_of_largest_area": best_time,
        "best_ratio": abs(best_ratio) if best_ratio else None,
    }

    # Compute final metrics for your top-level scatter:
    #  - relative max contraction = absolute value of (l(t)-l(0))/l(0)
    #  - relative max area = (a_max - a0)/a_max
    rel_max_contraction = abs(best_ratio) if best_ratio is not None else None
    if area_at_t0 is not None and best_area > 0:
        rel_max_area = (best_area - area_at_t0) / best_area
    else:
        rel_max_area = None

    analysis_results["rel_max_contraction"] = rel_max_contraction
    analysis_results["rel_max_area"] = rel_max_area

    overall_max_modification = None
    if rel_max_contraction is not None and rel_max_area is not None:
        overall_max_modification = max(rel_max_contraction, rel_max_area)
    analysis_results["overall_max_modification"] = overall_max_modification

    # Save everything to a JSON file in the 'plots' folder
    out_json = os.path.join(pattern_folder, 'plots', dump_file + '_analysis.json')
    with open(out_json, 'w') as f:
        json.dump(analysis_results, f, indent=2)

    # 3) MAKE A PLOT: HOLES vs TIME
    fig1, ax1 = plt.subplots()
    for label, times_dict in label_areas.items():
        times_sorted = sorted(times_dict.keys())
        areas_sorted = [times_dict[t] for t in times_sorted]
        ax1.plot(times_sorted, areas_sorted, '-', linewidth=1, alpha=0.2) 
        # The 'o' marker is optional. 
        # alpha=0.2 => lines are quite transparent

    ax1.set_xlabel('Time')
    ax1.set_ylabel('Hole area')
    ax1.set_title('Holes vs. Time')

    # Highlight the hole with the largest modification (best_label)
    if best_label in label_areas:
        best_times = sorted(label_areas[best_label].keys())
        best_areas = [label_areas[best_label][t] for t in best_times]
        ax1.plot(best_times, best_areas, '-', linewidth=2, color='red',
                label=f'Hole {best_label} (max mod)')
        ax1.legend()

    fig1.tight_layout()
    fig1.savefig(os.path.join(pattern_folder, 'plots', dump_file + '_holes_vs_time.png'), dpi=300)

    # 4) MAKE A PLOT: YARN LENGTHS vs TIME
    fig2, ax2 = plt.subplots()
    for mol, times_dict in yarn_lengths.items():
        times_sorted = sorted(times_dict.keys())
        lengths_sorted = np.array([times_dict[t] for t in times_sorted])
        ax2.plot(times_sorted, (lengths_sorted-lengths_sorted[0])/lengths_sorted[0], '-', linewidth=1, alpha=0.2)

    ax2.set_xlabel('Time')
    ax2.set_ylabel('Length')
    ax2.set_title('Yarn Lengths vs. Time')
    fig2.tight_layout()
    fig2.savefig(os.path.join(pattern_folder, 'plots', dump_file +'_yarns_vs_time.png'), dpi=300)

    # Highlight the yarn with maximum contraction (best_mol)
    if best_mol in yarn_lengths:
        times = sorted(yarn_lengths[best_mol].keys())
        lengths = np.array([yarn_lengths[best_mol][t] for t in times])
        ax2.plot(times, (lengths - lengths[0]) / lengths[0], '-', linewidth=2, color='red',
                 label=f'Yarn {best_mol} (max contraction)')
        ax2.legend()

    plt.close(fig1)
    plt.close(fig2)

def aggregate_results(output_folder='output/simulations', merged_csv='all_results.csv'):
    pattern_dirs = [
        d for d in os.listdir(output_folder)
        if d.startswith("Pattern_") and os.path.isdir(os.path.join(output_folder, d))
    ]

    all_dataframes = []

    for pattern_dir in pattern_dirs:
        folder_path = os.path.join(output_folder, pattern_dir)
        json_files = glob.glob(os.path.join(folder_path, 'plots', '*_analysis.json'))
        
        for json_file in json_files:
            # 1) Load the single-dict JSON manually
            with open(json_file, 'r') as f:
                data_dict = json.load(f)
            
            # 2) Convert that single dictionary to a one-row DataFrame
            #    json_normalize handles any nested fields (if they exist).
            df = pd.json_normalize(data_dict)

            df['pattern_folder'] = pattern_dir
            df['source_file'] = os.path.basename(json_file)
            all_dataframes.append(df)

    if not all_dataframes:
        print("No JSON files found. Nothing to aggregate.")
        return
    
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    combined_df.to_csv(merged_csv, index=False)
    print(f"Aggregated results written to {merged_csv}")

if __name__ == "__main__":
    base_dir = os.path.join('output', 'simulations')
    # Loop over all subfolders: Pattern_*
    for pattern_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, pattern_name)
        if not os.path.isdir(folder_path) or not pattern_name.startswith("Pattern_"):
            continue

        plots_folder = os.path.join(folder_path, "plots")
        if not os.path.exists(plots_folder):
           continue

        # Look for all JSON files ending in "_holes.json"
        for filename in os.listdir(plots_folder):
            if filename.endswith("_holes.json"):
                dump_file = filename[:-11]  # remove "_holes.json"
                holes_file = os.path.join(plots_folder, filename)
                yarn_file  = os.path.join(plots_folder, dump_file + "_yarn_lengths.json")
                if os.path.exists(yarn_file):
                    # We can now run the analysis
                    analyze_files(holes_file, yarn_file, folder_path, dump_file)
                else:
                    # If the corresponding yarn file is missing, skip
                    pass
    aggregate_results()