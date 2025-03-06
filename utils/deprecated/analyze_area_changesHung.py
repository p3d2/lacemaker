import json
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def analyze_area_changes(base_dir="output/simulations"):
    # Get all pattern folders (assumed to start with "Pattern_3025")
    pattern_dirs = [
        d for d in os.listdir(base_dir)
        if d.startswith("Pattern_3051") and os.path.isdir(os.path.join(base_dir, d))
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
                    # -----------------------------------------
        # COPY of your code, but for HUNGARIAN labels
        label_areas_hung = {}
        times_set_hung = set()

        # total_area_unlabelled is the same, so we reuse it or build the same logic:
        # total_area_unlabelled[t] = sum(frame.get("areas", []))  # already in your code.

        for frame in holes_data:
            t = frame["time"]
            times_set_hung.add(t)
            
            # Hungarian-labeled holes appear in "areas_label_hungarian"
            for hole in frame.get("areas_label_hungarian", []):
                h_label = hole["hungarian_label"]
                area = hole["area"]
                if h_label not in label_areas_hung:
                    label_areas_hung[h_label] = {}
                label_areas_hung[h_label][t] = area

        times_list_hung = sorted(times_set_hung)
        if not times_list_hung:
            # If no Hungarian-labeled holes exist, skip
            continue

        # Now replicate your baseline label selection logic,
        # i.e., "Select baseline labels: those present at t=0 with A0 > 0"
        baseline_labels_hung = {
            label: data[0] for label, data in label_areas_hung.items()
            if 0 in data and data[0] > 0
        }

        # Then replicate your valid-labels filtering. E.g.:
        valid_labels_hung = {}
        for label, A0 in baseline_labels_hung.items():
            if all(t in label_areas_hung[label] for t in times_list_hung):
                valid_labels_hung[label] = A0

        # Possibly pick the top 12 largest baseline holes:
        sorted_labels_hung = sorted(valid_labels_hung.items(), key=lambda kv: kv[1], reverse=True)
        biggest_labels_hung = [label for label, A0 in sorted_labels_hung[:12]]

        if not valid_labels_hung:
            print(f"No valid Hungarian labels in {holes_file}; skipping Hungarian analysis.")
            # You might `continue` or just skip the plot, but let's skip for clarity:
            continue

        # Compute relative changes for each valid Hungarian-labeled hole
        hole_relative_changes_hung = {}
        hole_absolute_changes_hung = {}
        for label, A0 in valid_labels_hung.items():
            rel_changes = {}
            abs_changes = {}
            for t in times_list_hung:
                A_t = label_areas_hung[label].get(t, 0.0)  # or handle missing
                abs_changes[t] = abs(A_t - A0)
                rel_changes[t] = abs(A_t - A0) / A0 if A0 != 0 else 0
            hole_relative_changes_hung[label] = rel_changes
            hole_absolute_changes_hung[label] = abs_changes

        # Compute overall relative area change, same fix if needed
        total_baseline_hung = sum(valid_labels_hung.values())
        overall_rel_hung = {}
        for t in times_list_hung:
            total_abs_change = sum(hole_absolute_changes_hung[label][t] for label in valid_labels_hung)
            overall_rel_hung[t] = total_abs_change / total_baseline_hung if total_baseline_hung else 0

        # loss_fraction is the same since it uses total_area_unlabelled minus tracked area in Hungarian
        loss_fraction_hung = {}
        for t in times_list_hung:
            T_tracked_hung = sum(
                label_areas_hung[label][t]
                for label in valid_labels_hung
                if t in label_areas_hung[label]
            )
            T_all = total_area_unlabelled[t]
            loss_fraction_hung[t] = (T_all - T_tracked_hung) / T_all if T_all > 0 else 0

        # Build a separate dictionary for Hungarian results
        hungarian_result = {
            "pattern": pattern_dir,
            "holes_file": os.path.basename(holes_file),
            "times": times_list_hung,
            "overall_relative_area_change": overall_rel_hung,
            "hole_relative_changes": hole_relative_changes_hung,
            "valid_labels": list(valid_labels_hung.keys()),
            "loss_fraction": loss_fraction_hung
        }

        # Write it out to a second analysis JSON
        analysis_filename_hung = os.path.join(
            plots_folder,
            os.path.basename(holes_file).replace("_holes.json", "_area_change_analysis_hungarian.json")
        )
        with open(analysis_filename_hung, 'w') as f:
            json.dump(hungarian_result, f, indent=2)

        # Plot Hungarian-based figure
        fig2, ax1_2 = plt.subplots()
        overall_values_hung = [overall_rel_hung[t] for t in times_list_hung]
        ax1_2.plot(times_list_hung, overall_values_hung, label="Overall (Hungarian)", linewidth=2, color="purple", zorder=-2)
        
        for label, rel_changes in hole_relative_changes_hung.items():
            rel_list = [rel_changes[t] for t in times_list_hung]
            ax1_2.plot(times_list_hung, rel_list, 'k-', linewidth=1, alpha=0.05, zorder=2)

        # highlight largest holes
        for label in biggest_labels_hung:
            rel_list = [hole_relative_changes_hung[label][t] for t in times_list_hung]
            ax1_2.plot(times_list_hung, rel_list, '-', linewidth=1, color="orange", alpha=0.5)
        
        ax1_2.set_xlabel("Time")
        ax1_2.set_ylabel("Relative Area Change (Hungarian)")
        ax1_2.set_title(f"Hungarian Relative Area Change – {pattern_dir}")
        ax1_2.legend(loc="upper left")

        # Secondary axis for hungarian loss fraction
        ax2_2 = ax1_2.twinx()
        loss_values_hung = [loss_fraction_hung[t] for t in times_list_hung]
        ax2_2.plot(times_list_hung, loss_values_hung, label="Loss Fraction (Hungarian)", linewidth=2, color="red")
        ax2_2.set_ylabel("Loss Fraction (Hungarian)")
        ax2_2.legend(loc="upper right")

        fig2.tight_layout()
        plot_filename_hung = os.path.join(
            plots_folder,
            os.path.basename(holes_file).replace("_holes.json", "_relative_area_change_hungarian.png")
        )
        fig2.savefig(plot_filename_hung, dpi=300)
        plt.close(fig2)

        # Optionally, add your Hungarian results to overall_results if you want them in the CSV
        # Just store them with e.g. a method="hungarian"
        for t in times_list_hung:
            overall_results.append({
                "pattern": pattern_dir,
                "holes_file": os.path.basename(holes_file),
                "time": t,
                "overall_relative_area_change": overall_rel_hung[t],
                "loss_fraction": loss_fraction_hung[t],
                "method": "hungarian"  # so you can differentiate in the CSV
            })
        # -----------------------------------------
    
    # Write aggregated overall area change data to a CSV file.
    if overall_results:
        df = pd.DataFrame(overall_results)
        output_csv = os.path.join(base_dir, "aggregated_area_change.csv")
        df.to_csv(output_csv, index=False)
        print(f"Aggregated area change results written to {output_csv}")
    else:
        print("No area change results found.")

if __name__ == "__main__":
    analyze_area_changes()
