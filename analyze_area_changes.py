import json
import os
import sys
import subprocess
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.patheffects as path_effects

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 9})

color1 = "#afd0ff"  # Dark Blue

def convert_fonts_to_outlines(input_pdf, output_pdf):
    gs_command = [
        "gs", "-o", output_pdf, "-sDEVICE=pdfwrite",
        "-dNoOutputFonts",  # convert fonts to outlines
        "-dNOPAUSE", "-dBATCH",
        "-dSubsetFonts=true", "-dEmbedAllFonts=true",
        "-dDownsampleColorImages=false",
        "-dDownsampleGrayImages=false",
        "-dDownsampleMonoImages=false",
        input_pdf
    ]
    subprocess.run(gs_command, check=True)

def analyze_area_changes(pattern, base_dir="output/simulations"):
    # Get all pattern folders (assumed to start with "Pattern_3025")
    bigA = False
    pattern_dirs = [
            d for d in os.listdir(base_dir)
            if d.startswith(pattern) and os.path.isdir(os.path.join(base_dir, d))
        ]
    overall_results = []  # For aggregated CSV output later

    for pattern_dir in pattern_dirs:
        folder_path = os.path.join(base_dir, pattern_dir)
        plots_folder = os.path.join(folder_path, "plots")
        if not os.path.exists(plots_folder):
            continue

        # Process every holes JSON file in this pattern folder.
        holes_files = glob.glob(os.path.join(plots_folder, "*_holes.json"))
        
        for holes_file in holes_files:
            with open(holes_file, 'r') as f:
                holes_data = json.load(f)
            
            # Build a dictionary mapping each label to its area at each time:
            # label_areas[label][t] = area
            label_areas = {}
            times_set = set()
            # Also build a dictionary for the unlabelled total area from "areas" field.
            total_area_unlabelled = {}
            for frame in holes_data:
                t = frame["time"]
                times_set.add(t)
                # Sum all areas (all detected holes, regardless of tracking)
                total_area_unlabelled[t] = sum(frame.get("areas", []))
                for hole in frame.get("areas_label", []):
                    label = hole["label"]
                    area = hole["area"]
                    if label not in label_areas:
                        label_areas[label] = {}
                    label_areas[label][t] = area
            times_list = sorted(times_set)
            if not times_list:
                continue

            # Select baseline labels: those that are present at t = 0 with A0 > 0.
            baseline_labels = {
                label: data[0] for label, data in label_areas.items()
                if 0 in data and data[0] > 0
            }
            # Use all baseline labels; later, if a label is missing in a frame it will be dropped for that frame.
            valid_labels = baseline_labels.copy()
              
            # Get the 12 biggest valid labels based on baseline area (A0)
            sorted_labels = sorted(valid_labels.items(), key=lambda kv: kv[1], reverse=True)
            biggest_labels = [label for label, A0 in sorted_labels[:12]]

            if not valid_labels:
                print(f"No valid labels (first & last frame) in {holes_file}; skipping analysis.")
                continue

            # For each valid label, compute the relative change (absolute difference relative to baseline):
            # r_i(t) = |A_i(t) - A_i(0)| / A_i(0)
            hole_relative_changes = {}
            for label, A0 in valid_labels.items():
                rel_changes = {}
                for t in times_list:
                    # If the label is missing in a middle frame, you could set its value to the last known area.
                    # Here, if a valid label is missing in a frame (which should be rare given our filter),
                    # we assign 0.
                    if t in label_areas[label]:
                        A_t = label_areas[label][t]
                        rel_changes[t] = abs(A_t - A0) / A0
                    else:
                        rel_changes[t] = 0.0
                hole_relative_changes[label] = rel_changes

            # Compute the overall relative area change for each frame using the baseline areas as weight:
            # R(t)= (sum_{i in valid_labels} |A_i(t)-A_i(0)|) / (sum_{i in valid_labels} A_i(0))
            # Compute the overall relative area change and loss fraction using a permanently pruned active set.
            active_set = set(valid_labels.keys())
            overall_rel = {}
            loss_fraction = {}
            for t in times_list:
                # Permanently drop labels missing at the current frame.
                active_set = {label for label in active_set if t in label_areas[label]}
                if active_set:
                    total_baseline_area = sum(valid_labels[label] for label in active_set)
                    total_abs_change = sum(valid_labels[label] * hole_relative_changes[label][t] for label in active_set)
                    overall_rel[t] = total_abs_change / total_baseline_area
                    T_tracked = sum(label_areas[label][t] for label in active_set)
                    T_all = total_area_unlabelled[t]
                    loss_fraction[t] = (T_all - T_tracked) / T_all if T_all > 0 else 0
                else:
                    overall_rel[t] = None
                    loss_fraction[t] = 0

            # ---------------------------------------
            # Read the yarn-lengths JSON
            #    to compute "max shrinking"
            # ---------------------------------------

            length_file = holes_file.replace("_holes.json", "_yarn_lengths.json")
            
            # We'll build a dict: time -> min (L(t)/L0)
            # If the file doesn't exist, store None or skip.
            time_to_min_ratio = {t: None for t in times_list}

            if os.path.isfile(length_file):
                with open(length_file, 'r') as lf:
                    length_data = json.load(lf)  # list of frames w/ 'time', 'yarns'
                
                # This structure parallels your "plot_lengths" script
                # We want (L(t)-L0)/L0 for each yarn, then take the min across yarns
                # at each time frame t.
                # First, gather data by (mol, type).
                yarn_dict = {}
                for frame in length_data:
                    t = frame['time']
                    for yd in frame['yarns']:
                        mol_val = yd['mol']
                        typ_val = yd['type']
                        length  = yd['length']
                        key = (mol_val, typ_val)
                        if key not in yarn_dict:
                            yarn_dict[key] = {
                                'times': [],
                                'lengths': []
                            }
                        yarn_dict[key]['times'].append(t)
                        yarn_dict[key]['lengths'].append(length)

                # For each (mol, type), find L0 = the length at the earliest time in that list
                # Then compute ratio = (L(t) - L0)/L0
                # We'll collect these in a dict: time -> [list of all ratio values].
                time_to_ratios = {t: [] for t in times_list}

                for (mol_val, typ_val), info in yarn_dict.items():
                    arr_times = np.array(info['times'])
                    arr_lens  = np.array(info['lengths'])
                    idx_sorted = np.argsort(arr_times)
                    t_sorted = arr_times[idx_sorted]
                    len_sorted = arr_lens[idx_sorted]

                    # Baseline length is the *first* entry in sorted order
                    L0 = len_sorted[0]
                    # Build a map from t -> ratio for these yarns
                    for i, tval in enumerate(t_sorted):
                        ratio = (len_sorted[i] - L0)/L0  # (L(t)-L0)/L0
                        # Save to time_to_ratios
                        if tval in time_to_ratios:
                            time_to_ratios[tval].append(ratio)

                # Now compute min ratio for each time
                for tval in times_list:
                    if time_to_ratios[tval]:
                        time_to_min_ratio[tval] = min(time_to_ratios[tval])
                    else:
                        time_to_min_ratio[tval] = None

            # Plot the overall relative area change (thick black line) and individual hole changes (thin lines)
            fig, ax1 = plt.subplots(figsize=(6 / 2.54, 5.6 / 2.54))
            overall_values = [overall_rel[t] for t in times_list]
            line1, = ax1.plot(times_list, overall_values, label="Weighted average", linewidth=2, color=color1, zorder=20)
            for label, rel_changes in hole_relative_changes.items():
                rel_list = [rel_changes[t] for t in times_list]
                ax1.plot(times_list, rel_list, 'k-', linewidth=0.5, alpha=0.05, zorder=2, rasterized=True)
            
            maxH = overall_values[10]
            tDelta = ax1.annotate(r'$\delta_{{10}}={:.2f}$'.format(maxH),  
                            xy=(10, maxH + 0.05),    
                            xytext=(10, maxH + 0.05),  
                            fontsize=10, color='black', ha='center', va='bottom') 
            
            tDelta.set_path_effects([path_effects.withStroke(linewidth=2, foreground="white")])

            if bigA == True:
                for label in biggest_labels:
                    rel_list = [hole_relative_changes[label][t] for t in times_list]
                    line3, = ax1.plot(times_list, rel_list, '-', linewidth=1, color="pink", alpha = 0.5, zorder = 5, label="12 Largest Areas" if label==biggest_labels[0] else "")
            ax1.set_xlabel(r'Simulation Step')
            ax1.set_ylabel("Relative Area Change")
            ax1.set_ylim((0,1))
            # Add a secondary axis for loss fraction in red.
            ax2 = ax1.twinx()
            loss_values = [loss_fraction[t] for t in times_list]
            line2, = ax2.plot(times_list, loss_values, '--', label="Loss Fraction", linewidth=1, color="red", zorder=10)
            ax2.set_ylabel("Loss Fraction")
            ax1.set_xlim((0,20))
            ax2.set_ylim((0,0.1))
    
            lines = [line1, line2]  # Collect handles
            labels = [l.get_label() for l in lines]  # Extract labels
            ax1.legend(lines, labels, loc="upper left", fontsize=7)  # Set combined legend

            ax1.xaxis.set_major_locator(ticker.MultipleLocator(10))   # Major ticks every 5 units
            ax1.xaxis.set_minor_locator(ticker.MultipleLocator(2))   # Minor ticks every 1 unit

            ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2)) # Major ticks every 0.5 units
            ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.02)) # Major ticks every 0.5 units
            ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
    

            fig.tight_layout()
            pdf_file = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_relative_area_change.pdf")
            )
            fig.savefig(pdf_file, dpi=300, bbox_inches='tight')
            plt.close(fig)

            compressed_pdf = pdf_file.replace(".pdf", "_c.pdf")
            convert_fonts_to_outlines(pdf_file, compressed_pdf)
            os.replace(compressed_pdf, pdf_file)

            # Record overall time series data for CSV aggregation.
            for t in times_list:
                overall_results.append({
                    "pattern": pattern_dir,
                    "holes_file": os.path.basename(holes_file),
                    "time": t,
                    "overall_relative_area_change": overall_rel[t],
                    "loss_fraction": loss_fraction[t],
                    # The "max shrinking" is the minimal ratio => the most negative
                    # or the biggest drop from baseline. We store it here:
                    "max_shrinking": abs(time_to_min_ratio[t])
                })
    
    # ---------------------------------------
    # 7) Write aggregated overall area change data (with max_shrinking) to CSV
    # ---------------------------------------
    if overall_results:
        df = pd.DataFrame(overall_results)
        output_csv = os.path.join(base_dir, "aggregated_area_change.csv")
        df.to_csv(output_csv, index=False)
        print(f"Aggregated results (including max shrinking) written to {output_csv}")
    else:
        print("No area change results found.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: Missing required argument.")
        print("Usage: python script.py <pattern_name>")
        sys.exit(1)  # Exit with error

    pattern_arg = sys.argv[1]
    analyze_area_changes(pattern=pattern_arg)
