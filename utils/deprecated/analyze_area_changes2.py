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

def analyze_area_changes2(base_dir="output/simulations"):
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
                    hist, _ = np.histogram(areas, bins=bins, weights=areas)
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

            # Compute weighted summary statistics for each frame.
            weighted_mean = {}
            weighted_variance = {}
            for t in times_list:
                mean, var = weighted_stats(distributions[t])
                weighted_mean[t] = mean
                weighted_variance[t] = var

            # Compute relative change from t = 0 for the weighted mean and variance.
            mean_rel_change = {}
            var_rel_change = {}
            baseline_mean = weighted_mean[times_list[0]]
            baseline_var = weighted_variance[times_list[0]]
            for t in times_list:
                mean_rel_change[t] = 0 if baseline_mean == 0 else (weighted_mean[t] - baseline_mean) / baseline_mean
                var_rel_change[t] = 0 if baseline_var == 0 else (weighted_variance[t] - baseline_var) / baseline_var

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

            # --- Additional Method 2: ECDF and KS Statistic ---
            # For each frame, compare the ECDF (via ks_2samp) of the current areas vs. baseline.
            ks_stat = {}
            for t in times_list:
                current_areas = distributions[t]
                if len(current_areas) == 0:
                    ks_stat[t] = 0.0
                else:
                    stat, _ = ks_2samp(baseline_areas, current_areas)
                    ks_stat[t] = stat

            # --- Additional Method 3: Quantile-Based Analysis ---
            # For each frame, compute the 25th percentile, median, and 75th percentile.
            quantiles = {}
            for t in times_list:
                areas = np.array(distributions[t])
                if areas.size == 0:
                    quantiles[t] = {"q25": 0, "median": 0, "q75": 0}
                else:
                    quantiles[t] = {
                        "q25": np.percentile(areas, 25),
                        "median": np.percentile(areas, 50),
                        "q75": np.percentile(areas, 75)
                    }
            # Compute relative changes of these quantiles from baseline.
            quantiles_rel = {}
            baseline_quant = quantiles[baseline_frame]
            for t in times_list:
                quant_rel = {}
                for key in ["q25", "median", "q75"]:
                    if baseline_quant[key] == 0:
                        quant_rel[key] = 0.0
                    else:
                        quant_rel[key] = (quantiles[t][key] - baseline_quant[key]) / baseline_quant[key]
                quantiles_rel[t] = quant_rel

            # Build the result dictionary.
            result = {
                "pattern": pattern_dir,
                "holes_file": os.path.basename(holes_file),
                "times": times_list,
                "distribution_distance": distribution_distance,  # L1 distance from baseline histogram.
                "weighted_mean": weighted_mean,
                "weighted_variance": weighted_variance,
                "mean_rel_change": mean_rel_change,
                "var_rel_change": var_rel_change,
                "wasserstein_distance": wasserstein_dist,
                "ks_statistic": ks_stat,
                "quantiles": quantiles,
                "quantiles_rel": quantiles_rel
            }

            # Write the analysis result to a JSON file.
            analysis_filename = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_area_change_analysis2.json")
            )
            with open(analysis_filename, 'w') as f:
                json.dump(result, f, indent=2)

            # ----------------- Plot 1: Summary Metrics -----------------
            # Plot distribution distance (black), and on a secondary axis plot the relative changes of weighted mean (blue)
            # and weighted variance (green).
            fig, ax1 = plt.subplots()
            ax1.plot(times_list, [distribution_distance[t] for t in times_list],
                     label="Distribution Distance", linewidth=2, color="black")
            ax1.set_xlabel("Time")
            ax1.set_ylabel("L1 Distance", color="black")
            ax1.tick_params(axis="y", labelcolor="black")
            ax1.set_title(f"Weighted Distribution Dynamics – {pattern_dir}")

            ax2 = ax1.twinx()
            ax2.plot(times_list, [mean_rel_change[t] for t in times_list],
                     label="Mean Rel Change", linewidth=2, color="blue")
            ax2.plot(times_list, [var_rel_change[t] for t in times_list],
                     label="Variance Rel Change", linewidth=2, color="green")
            ax2.set_ylabel("Relative Change", color="blue")
            ax2.tick_params(axis="y", labelcolor="blue")
            ax1.legend(loc="upper left")
            ax2.legend(loc="upper right")

            fig.tight_layout()
            plot_filename = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_area_change2.png")
            )
            fig.savefig(plot_filename, dpi=300)
            plt.close(fig)

            # ----------------- Plot 2: Histogram Evolution -----------------
            # Create a 2D array where each row corresponds to the normalized weighted histogram at a given time.
            bin_centers = (bins[:-1] + bins[1:]) / 2.0
            H = np.array([histograms[t] for t in times_list])
            fig, ax = plt.subplots(figsize=(10, 6))
            # Use imshow to display the evolution: x-axis = area (bin centers), y-axis = time.
            c = ax.imshow(H, aspect='auto', origin='lower',
                          extent=[bin_centers[0], bin_centers[-1], times_list[0], times_list[-1]],
                          interpolation='nearest', cmap='viridis')
            ax.set_xlabel("Area")
            ax.set_ylabel("Time")
            ax.set_title(f"Evolution of Weighted Area Distribution – {pattern_dir}")
            fig.colorbar(c, ax=ax, label="Normalized Weighted Frequency")
            fig.tight_layout()
            hist_plot_filename = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_histogram_evolution.png")
            )
            fig.savefig(hist_plot_filename, dpi=300)
            plt.close(fig)

            # Record overall time series data for CSV aggregation.
            for t in times_list:
                overall_results.append({
                    "pattern": pattern_dir,
                    "holes_file": os.path.basename(holes_file),
                    "time": t,
                    "distribution_distance": distribution_distance[t],
                    "mean_rel_change": mean_rel_change[t],
                    "var_rel_change": var_rel_change[t],
                    "weighted_mean": weighted_mean[t],
                    "weighted_variance": weighted_variance[t],
                })
            
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

            # --- Plot: KS Statistic Evolution ---
            fig, ax = plt.subplots()
            ks_vals = [ks_stat[t] for t in times_list]
            ax.plot(times_list, ks_vals, label="KS Statistic", linewidth=2, color="orange")
            ax.set_xlabel("Time")
            ax.set_ylabel("KS Statistic")
            ax.set_title(f"KS Statistic Evolution – {pattern_dir}")
            ax.legend()
            fig.tight_layout()
            ks_plot_filename = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_ks_statistic.png")
            )
            fig.savefig(ks_plot_filename, dpi=300)
            plt.close(fig)

            # --- Plot: Quantile Evolution (Relative Change) ---
            fig, ax = plt.subplots()
            med_vals = [quantiles_rel[t]["median"] for t in times_list]
            q25_vals = [quantiles_rel[t]["q25"] for t in times_list]
            q75_vals = [quantiles_rel[t]["q75"] for t in times_list]
            ax.plot(times_list, med_vals, label="Median Rel Change", linewidth=2, color="blue")
            ax.plot(times_list, q25_vals, label="25th Percentile Rel Change", linewidth=2, color="green")
            ax.plot(times_list, q75_vals, label="75th Percentile Rel Change", linewidth=2, color="red")
            ax.set_xlabel("Time")
            ax.set_ylabel("Relative Change")
            ax.set_title(f"Quantile Evolution – {pattern_dir}")
            ax.legend()
            fig.tight_layout()
            quant_plot_filename = os.path.join(
                plots_folder,
                os.path.basename(holes_file).replace("_holes.json", "_quantile_evolution.png")
            )
            fig.savefig(quant_plot_filename, dpi=300)
            plt.close(fig)



    
    # Write aggregated overall area change data to a CSV file.
    if overall_results:
        df = pd.DataFrame(overall_results)
        output_csv = os.path.join(base_dir, "aggregated_area_change2.csv")
        df.to_csv(output_csv, index=False)
        print(f"Aggregated area change2 results written to {output_csv}")
    else:
        print("No area change2 results found.")

if __name__ == "__main__":
    analyze_area_changes2()
