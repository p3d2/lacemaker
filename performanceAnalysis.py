import json
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



sample_map = {
    "1084": ("Square-2", "S2", "s", "FFD457"),
    "1092": ("Square-1", "S1", "s", "FFF203"),
    "2001": ("Diamond-1", "D1", "D", "BCE4E5"),
	"2002": ("Diamond-2", "D2", "D", "CFE0F3"),
	"2009": ("Diamond-3", "D3", "D", "72CDF4"),
	"2023": ("Diamond-7", "D7", "D", "002D56"),
	"2048": ("Hexa-1", "H1", "h", "A9D08E"),
	"3001": ("Rose-8", "R8", "o", "ED2024"),
	"3003": ("Rose-9", "R9", "o", "A84499"),
	"3013": ("Rose-6", "R6", "o", "B9529F"),
	"3018": ("Rose-7", "R7", "o", "EE3C96"),
	"3021": ("Rose-3", "R3", "o", "F69799"),
	"3023": ("Rose-1", "R1", "o", "EECDE2"),
	"3023v2": ("Rose-1", "R1", "o", "EECDE2"),
	"3024": ("Rose-4", "R4", "o", "E26FAB"),
	"3025": ("Rose-2", "R2", "o", "DB9FC7"),
	"3026": ("Rose-5", "R5", "o", "F26722"),
	"3034": ("Rose-10", "R10", "o", "CB1C68"),
	"3043": ("Rose-11", "R11", "o", "642265"),
	"3051": ("Rose-12", "R12", "o", "92278F"),
	"3051_1t": ("Rose-12", "R12", "o", "92278F"),
	"3051v2": ("Rose-12", "R12", "o", "92278F"),
	"3052": ("Rose-13", "R13", "o", "CB3727"),
	"3059": ("Spruce-1", "Sp1", "^", "008752"),
	"7000": ("Diamond-6", "D6", "D", "0065A4"),
	"7001": ("Diamond-5", "D5", "D", "007DC3"),
	"7002": ("Diamond-4", "D4", "D", "00B5CC"),
}

def analyze_files(area_csv_file, yarn_file, pattern_folder, base_filename, exp_str):
	# --- Read the area analysis CSV file ---
	csv_df = pd.read_csv(area_csv_file)
	# Filter rows for the current experiment based on the holes_file naming pattern.
	matching_rows = csv_df[csv_df['holes_file'].str.contains(f"contract_fix_{exp_str}_holes.json")]
	if matching_rows.empty:
		raise FileNotFoundError(
			f"No matching area analysis rows found for experiment {exp_str} in file {area_csv_file}"
		)
	# Extract rows for time 10 and time 20.
	t10_rows = matching_rows[matching_rows['time'] == 10]
	t20_rows = matching_rows[matching_rows['time'] == 20]
	if t10_rows.empty or t20_rows.empty:
		raise FileNotFoundError(
			f"Area analysis rows for time 10 or 20 missing for experiment {exp_str} in file {area_csv_file}"
		)
	area_change_t10 = t10_rows.iloc[0]['overall_relative_area_change']
	area_change_t20 = t20_rows.iloc[0]['overall_relative_area_change']
	
	# --- Process Yarn Lengths from the yarn JSON file ---
	with open(yarn_file, 'r') as f:
		yarn_data = json.load(f)
	
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
	
	# Find the yarn with maximum contraction (minimum relative change)
	best_ratio = None
	best_mol   = None
	best_time_yarn = None
	for mol, times_dict in yarn_lengths.items():
		if 0 not in times_dict:
			continue
		l0 = times_dict[0]
		if l0 == 0:
			continue
		for t, length_at_t in times_dict.items():
			ratio = (length_at_t - l0) / l0
			if l0 > 64.0:
				if (best_ratio is None) or (ratio < best_ratio):
					best_ratio = ratio
					best_mol   = mol
					best_time_yarn = t
	rel_max_contraction = best_ratio if best_ratio is not None else None
	
	filename = os.path.basename(yarn_file)
	prefix = filename.split('_')[0]  
	middle = filename.split("contract_fix_")[1].split("_yarn")[0]

	exp_id = prefix[7:]
	full_name, reduced_name, marker, color = sample_map.get(exp_id, ("UNKNOWN", "UNK", "o", "k"))

	# --- Build the analysis results dictionary ---
	analysis_results = {
		"exp": prefix,
		"sample_full": full_name,
    	"sample_reduced": reduced_name,
		"contract": middle,
		"area_change_t10": area_change_t10,
		"area_change_t20": area_change_t20,
		"rel_max_contraction": rel_max_contraction,
		"best_mol": best_mol,
		"best_time_yarn": best_time_yarn,
		"marker": marker,
		"color": color,
	}
	
	# --- Create a plot for Yarn Lengths vs. Time ---
	fig, ax = plt.subplots()
	for mol, times_dict in yarn_lengths.items():
		times_sorted = sorted(times_dict.keys())
		lengths_sorted = np.array([times_dict[t] for t in times_sorted])
		ax.plot(times_sorted, (lengths_sorted - lengths_sorted[0]) / lengths_sorted[0], '-', linewidth=1, alpha=0.2)
	ax.set_xlabel('Time')
	ax.set_ylabel('Relative Yarn Length Change')
	ax.set_title('Yarn Lengths vs. Time')
	if best_mol in yarn_lengths:
		times = sorted(yarn_lengths[best_mol].keys())
		lengths = np.array([yarn_lengths[best_mol][t] for t in times])
		ax.plot(times, (lengths - lengths[0]) / lengths[0], '-', linewidth=2, color='red',
				label=f'Yarn {best_mol} (max contraction)')
		ax.legend()
	fig.tight_layout()
	fig.savefig(os.path.join(pattern_folder, 'plots', base_filename + '_yarns_vs_time.png'), dpi=300)
	plt.close(fig)
	
	# --- Save the analysis results to a JSON file ---
	out_json = os.path.join(pattern_folder, 'plots', base_filename + '_analysis.json')
	with open(out_json, 'w') as f:
		json.dump(analysis_results, f, indent=2)
	
	# Append the generated file name to a manifest file
	manifest_path = os.path.join(pattern_folder, 'plots', 'analysis_manifest.txt')
	with open(manifest_path, 'a') as mf:
		mf.write(os.path.basename(out_json) + "\n")

def aggregate_results(output_folder='output/simulations', merged_csv='all_results.csv'):
	pattern_dirs = [
		d for d in os.listdir(output_folder)
		if d.startswith("Pattern_") and os.path.isdir(os.path.join(output_folder, d))
	]

	all_dataframes = []
	desired_keys = ["exp", "sample_full", "sample_reduced", "contract", "area_change_t10", "area_change_t20", 
					"rel_max_contraction", "best_mol", "best_time_yarn", "marker", "color"]

	for pattern_dir in pattern_dirs:
		folder_path = os.path.join(output_folder, pattern_dir)
		plots_folder = os.path.join(folder_path, 'plots')
		manifest_file = os.path.join(plots_folder, 'analysis_manifest.txt')
		
		if os.path.exists(manifest_file):
			with open(manifest_file, 'r') as mf:
				file_list = [line.strip() for line in mf if line.strip()]
			# Build full paths for the files in the manifest
			json_files = [os.path.join(plots_folder, fname) for fname in file_list]
		else:
			# Fallback: only include files matching a specific pattern (if no manifest exists)
			continue
		
		for json_file in json_files:
			with open(json_file, 'r') as f:
				data_dict = json.load(f)
			# Create one row per JSON file without exploding nested lists:
			df = pd.DataFrame([data_dict])
			df['pattern_folder'] = pattern_dir
			df['source_file'] = os.path.basename(json_file)
			# Keep only the desired keys (if present)
			df = df[[col for col in desired_keys if col in df.columns]]
			all_dataframes.append(df)

	if not all_dataframes:
		print("No JSON files found. Nothing to aggregate.")
		return

	combined_df = pd.concat(all_dataframes, ignore_index=True)
	combined_df.to_csv(merged_csv, index=False)
	print(f"Aggregated results written to {merged_csv}")

if __name__ == "__main__":
	# --- Read the configuration TSV file ---
	config_file = "assets/data/jobs.tsv"  
	if not os.path.exists(config_file):
		raise FileNotFoundError(f"Configuration TSV file not found: {config_file}")
	config_df = pd.read_csv(config_file, sep='\t')
    
	base_dir = os.path.join('output', 'simulations')
	pattern_dirs = [
		d for d in os.listdir(base_dir)
		if d.startswith("Pattern_") and os.path.isdir(os.path.join(base_dir, d))
	]

	for pattern_dir in pattern_dirs:
		plots_folder = os.path.join(base_dir, pattern_dir, "plots")
		manifest_path = os.path.join(plots_folder, 'analysis_manifest.txt')
		if os.path.exists(manifest_path):
			os.remove(manifest_path)

	for idx, row in config_df.iterrows():
		if int(row['analyse']) != 1:
			continue
		pattern_id = str(row['pattern_id'])
		# Build the pattern folder (e.g., "Pattern_1084")
		folder_path = os.path.join(base_dir, f"Pattern_{pattern_id}")
		if not os.path.isdir(folder_path):
			raise FileNotFoundError(f"Pattern folder missing: {folder_path}")
		plots_folder = os.path.join(folder_path, "plots")
		if not os.path.exists(plots_folder):
			raise FileNotFoundError(f"Plots folder missing: {plots_folder}")
		
		# Read the area analysis CSV file (e.g., "Pattern_1084_areaResults.csv")
		area_results_csv = os.path.join(plots_folder, f"Pattern_{pattern_id}_areaResults.csv")
		if not os.path.exists(area_results_csv):
			raise FileNotFoundError(f"Area results CSV missing: {area_results_csv}")
		
		# Parse the experiments from the 'exp' column (a JSON-formatted string)
		try:
			experiments = json.loads(row['exp'])
		except Exception as e:
			raise ValueError(f"Error parsing 'exp' for pattern {pattern_id}: {e}")
		
		for exp in experiments:
			# Create a string representation for the experiment, e.g. "1_2_3_4"
			exp_str = '_'.join(map(str, exp))
			# Define a base filename for output (this will be used in naming the JSON and plot files)
			base_filename = f"pattern{pattern_id}_contract_fix_{exp_str}"
			
			# Find the yarn lengths JSON file using a glob pattern.
			yarn_pattern = f"pattern{pattern_id}*contract_fix_{exp_str}_yarn_lengths.json"
			yarn_matches = glob.glob(os.path.join(plots_folder, yarn_pattern))
			if len(yarn_matches) == 0:
				raise FileNotFoundError(
					f"Yarn lengths file missing for pattern {pattern_id} experiment {exp_str} using pattern: {yarn_pattern}"
				)
			if len(yarn_matches) > 1:
				raise RuntimeError(
					f"Multiple yarn lengths files found for pattern {pattern_id} experiment {exp_str}: {yarn_matches}"
				)
			yarn_file = yarn_matches[0]
			
			# Run the analysis: use the area results CSV (instead of holes JSON), the yarn file, etc.
			analyze_files(area_results_csv, yarn_file, folder_path, base_filename, exp_str)
	aggregate_results()