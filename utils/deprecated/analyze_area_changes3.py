#!/usr/bin/env python3
"""
analyze_area_change2.py

This script analyzes the dynamics of unlabelled hole areas (from the "areas" field)
by computing a weighted distribution at each time point. The area values are used as the weight.
For each frame, the script computes:
  1. A normalized weighted histogram (using fixed global bins).
  2. The L₁ distance (normalized) between the histogram at time 0 (baseline) and each subsequent frame.
  3. Weighted summary statistics (mean and variance) using weight = area, and then their relative changes from t = 0.
  
In addition, this version produces a second plot that shows a heatmap of the histogram evolution over time.

Results are saved to JSON files, aggregated into a CSV, and plots are generated.
"""

import json
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde, wasserstein_distance, ks_2samp

def weighted_stats(areas):
    """Compute weighted mean and weighted variance for a list of areas,
    where each sample is weighted by its area.
    Returns (weighted_mean, weighted_variance)."""
    areas = np.array(areas)
    if areas.size == 0 or np.sum(areas) == 0:
        return 0, 0
    # Weighted mean: sum(a_i^2) / sum(a_i)
    mean = np.sum(areas**2) / np.sum(areas)
    # Weighted variance: sum(a_i * (a_i - mean)^2) / sum(a_i)
    var = np.sum(areas * (areas - mean)**2) / np.sum(areas)
    return mean, var

def analyze_area_changes3(base_dir="output/simulations"):
    # Get all pattern folders (assumed to start with "Pattern_3025")
    pattern_dirs = [
        d for d in os.listdir(base_dir)
        if d.startswith("Pattern_3025") and os.path.isdir(os.path.join(base_dir, d))
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
            
            # Build a dictionary: distributions[t] = list of unlabelled areas at that frame.
            distributions = {}
            times_set = set()
            for frame in holes_data:
                t = frame["time"]
                times_set.add(t)
                distributions[t] = frame.get("areas", [])
            times_list = sorted(times_set)
            if not times_list:
                continue

            # Combine all areas (from all frames) to define global histogram bins.
            all_areas = []
            for t in times_list:
                all_areas.extend(distributions[t])
            if len(all_areas) == 0:
                continue

            global_min = min(all_areas)
            global_max = max(all_areas)
            num_bins = 50
            bins = np.linspace(global_min, global_max, num_bins+1)

            # For each frame, compute weighted histogram (weights = area) and normalize.
            histograms = {}
            for t in times_list:
                areas = np.array(distributions[t])
                if areas.size == 0:
                    hist = np.zeros(num_bins)
                else:
                    hist, _ = np.histogram(areas, bins=bins)
                    total = np.sum(hist)
                    hist = hist / total if total > 0 else np.zeros(num_bins)
                histograms[t] = hist

            # Baseline histogram at t=0.
            baseline_hist = histograms[times_list[0]]
            # Compute normalized L1 (Manhattan) distance between histograms.
            distribution_distance = {}
            for t in times_list:
                dist = np.sum(np.abs(histograms[t] - baseline_hist)) / 2.0
                distribution_distance[t] = dist

            # Use the raw unlabelled areas (the distributions dict) to compute the weighted Wasserstein distance.
            baseline_frame = times_list[0]
            baseline_areas = np.array(distributions[baseline_frame])
            # Use the areas themselves as weights.
            baseline_weights = baseline_areas  
            wasserstein_dist = {}
            for t in times_list:
                current_areas = np.array(distributions[t])
                if current_areas.size == 0:
                    wasserstein_dist[t] = 0.0
                else:
                    current_weights = current_areas
                    dist_val = wasserstein_distance(baseline_areas, current_areas,
                                                    u_weights=baseline_weights, v_weights=current_weights)
                    wasserstein_dist[t] = dist_val

            # Build the result dictionary.
            result = {
                "pattern": pattern_dir,
                "holes_file": os.path.basename(holes_file),
                "times": times_list,
                "distribution_distance": distribution_distance,  # L1 distance from baseline histogram.
                "wasserstein_distance": wasserstein_dist,
            }

            # Write the analysis result to a JSON file.
            analysis_filename = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_area_change_analysis2.json")
            )
            with open(analysis_filename, 'w') as f:
                json.dump(result, f, indent=2)

            # ----------------- Plot 2: Histogram Evolution -----------------
            # Create a 2D array where each row corresponds to the normalized weighted histogram at a given time.
            bin_centers = (bins[:-1] + bins[1:]) / 2.0
            H = np.array([histograms[t] for t in times_list])
            fig, ax = plt.subplots(figsize=(10, 6))

            # Calculate an appropriate vertical offset for each ridge
            y_offset = np.max(H) * 0.6

            # Plot each histogram as a ridge, offset by time order.
            for i, t in enumerate(times_list):
                # Add an offset to create the ridge effect
                ridge = H[i] + i * y_offset
                ax.plot(bin_centers, ridge, color='C0', lw=1)
                ax.fill_between(bin_centers, i * y_offset, ridge, color='C0', alpha=0.5)
                # Optionally annotate the ridge with its time label
                ax.text(bin_centers[0] - 0.05*(bin_centers[-1]-bin_centers[0]), i * y_offset,
                        f"{t}", va='center', ha='right', fontsize=8)

            ax.set_xlabel("Area")
            ax.set_ylabel("Histogram Ridge (Offset by Time)")
            ax.set_title(f"Evolution of Weighted Area Distribution – {pattern_dir}")
            ax.set_ylim(-y_offset, len(times_list) * y_offset)
            ax.set_yticks([])  # Hide the y-axis ticks for a cleaner look
            fig.tight_layout()

            hist_plot_filename = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_histogram_ridgeplot.png")
            )
            fig.savefig(hist_plot_filename, dpi=300)
            plt.close(fig)

            
            # --- Plot: Wasserstein Distance Evolution ---
            fig, ax = plt.subplots()
            wass_vals = [wasserstein_dist[t] for t in times_list]
            ax.plot(times_list, wass_vals, label="Wasserstein Distance", linewidth=2, color="purple")
            ax.set_xlabel("Time")
            ax.set_ylabel("Wasserstein Distance")
            ax.set_title(f"Wasserstein Distance Evolution – {pattern_dir}")
            ax.legend()
            fig.tight_layout()
            wass_plot_filename = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_wasserstein_distance.png")
            )
            fig.savefig(wass_plot_filename, dpi=300)
            plt.close(fig)
    
    # Write aggregated overall area change data to a CSV file.
    if overall_results:
        df = pd.DataFrame(overall_results)
        output_csv = os.path.join(base_dir, "aggregated_area_change3.csv")
        df.to_csv(output_csv, index=False)
        print(f"Aggregated area change2 results written to {output_csv}")
    else:
        print("No area change2 results found.")

if __name__ == "__main__":
    analyze_area_changes3()