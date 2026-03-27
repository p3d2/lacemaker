# Lacemaker

Software to simulate and analyse the mechanical behaviour of bobbin lace patterns under selective yarn contraction. Patterns are modelled as bead-spring networks and simulated with LAMMPS. Pore geometry is tracked over time and compared against physical experiments.

Associated article: [Title] — [Authors], [Journal], [Year], DOI: [DOI]

---

## Requirements

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate lacemaker
```

Additional requirements:
- **LAMMPS** — compiled with the custom pair style from `utils/pair_soft_exclude.cpp` (see [Custom pair potential](#custom-pair-potential))
- **Ghostscript** — used by figure scripts to convert PDF fonts to outlines (`gs` must be on PATH)

---

## Repository structure

```
lacemaker/
+-- lace_maker.py            Convert JSON pattern to LAMMPS data file
+-- run_jobs.py              Orchestrate full simulation pipeline (SLURM)
+-- sim_holes_extract.py     Extract pore geometry from simulation trajectories
+-- assets/data/jobs.tsv     Simulation parameter table (one row per pattern)
+-- all_results.csv          Aggregated simulation results (produced by performanceAnalysis.py)
+-- input/
|   +-- json_patterns/       Pattern geometry definitions (JSON)
|   \-- lammps_input/        LAMMPS input scripts (in_jobs8.lmp, in_contract_jobs8.lmp)
+-- utils/
|   +-- vid_holes_extract.py   Extract pore geometry from experimental videos
|   +-- holes_analysis.py      KDE analysis and plots of pore size distributions
|   +-- analyze_area_changes.py  Compute relative area changes over time
|   +-- performanceAnalysis.py   Aggregate results across all patterns -> all_results.csv
|   +-- lengths_analysis.py    Per-pattern yarn length plots
|   +-- fig2.py              Reproduce Figure 2 (experimental results)
|   +-- fig3_bar.py          Reproduce Figure 3 (simulation bar plots)
|   +-- genResults.py        Generate interactive 3D results plot
|   +-- pair_soft_exclude.cpp  Custom LAMMPS pair potential (source)
|   +-- pair_soft_exclude.h
|   \-- deprecated/          Unused/superseded scripts
+-- guides/
|   +-- create_json.md       Guide for defining new patterns
|   \-- diagrams.md          Pattern diagram conventions
```

---

## Pipeline

### Step 1 — Define a pattern

Pattern geometries are defined as JSON files in `input/json_patterns/`. See [guides/create_json.md](guides/create_json.md) for the full format specification.

### Step 2 — Generate LAMMPS input

Convert a JSON pattern to a LAMMPS data file:

```bash
python lace_maker.py input/json_patterns/pattern3023v2.json \
  --dist_particles=0.25 --units=1 --mass=1 --threshold=0 --ks1=30 --kb=0.1
```

Output: `output/lammps_data/pattern3023v2_0.25_1.0_0.0_30.0_0.1.data`

The parameters used for each pattern in the article are listed in `assets/data/jobs.tsv`.

### Step 3 — Run simulations

`run_jobs.py` reads `assets/data/jobs.tsv` and submits the full pipeline to a SLURM cluster:

```bash
python run_jobs.py assets/data/jobs.tsv
```

This runs for each pattern with `analyse=1` in the TSV:
1. Mesh generation via `lace_maker.py`
2. Stabilisation simulation (`input/lammps_input/in_jobs8.lmp`)
3. Contraction simulations for each yarn combination (`input/lammps_input/in_contract_jobs8.lmp`)
4. Hole extraction via `sim_holes_extract.py`

Output is written to `output/simulations/Pattern_{ID}/`.

> **Alternative:** Simulation outputs (trajectories, hole JSON, yarn lengths) can be downloaded directly from Zenodo (DOI: [Zenodo DOI]) instead of running Steps 1-4. This allows reproducing the analysis and figures without access to a cluster.

### Step 4 — Extract pore geometry from simulations

If running the analysis step independently:

```bash
python sim_holes_extract.py output/simulations/Pattern_3023v2
```

### Step 5 — Extract pore geometry from experimental videos

```bash
python utils/vid_holes_extract.py path/to/video.MOV
```

### Step 6 — Compute pore area statistics

```bash
python utils/holes_analysis.py output/simulations/Pattern_3023v2/plots
```

Produces `*_holes.json` and distribution plots inside the `plots/` subfolder.

### Step 7 — Compute area changes over time

```bash
python utils/analyze_area_changes.py output/simulations/Pattern_3023v2/plots
```

### Step 8 — Aggregate results across all patterns

```bash
python utils/performanceAnalysis.py
```

Reads all per-experiment `*_analysis.json` files and produces `all_results.csv`, which is required by the figure scripts.

### Step 9 — Yarn length plots (per pattern)

```bash
python utils/lengths_analysis.py Pattern_3023v2
```

Output: `output/simulations/Pattern_3023v2/plots/*_lengths_subplots.pdf`

---

## Reproducing the figures

### Figure 2 — experimental results

```bash
python utils/fig2.py
```

Output: `fig.pdf`. Performance values are entered manually in the script's `data` list.

### Figure 3 — simulation bar plots

Requires `all_results.csv` (produced by Step 8 above, or available in the repository).

```bash
python utils/fig3_bar.py
```

Output: `barplots_part1.pdf`, `barplots_part2.pdf`

### Interactive 3D results

Requires `all_results.csv`.

```bash
python utils/genResults.py
```

Output: `my_interactive_3d_plot.html` (open in a browser)

---

## Custom pair potential

The simulations use a custom soft-exclude pair potential defined in `utils/pair_soft_exclude.cpp`. To compile it into LAMMPS, follow the LAMMPS documentation for [adding a custom pair style](https://docs.lammps.org/Developer_write_pair.html) and copy the `.cpp` and `.h` files into the LAMMPS `src/` directory before building.

---

## Dataset

Simulation trajectories and experimental data are available on Zenodo: DOI [Zenodo DOI]

See [README_dataset.md](README_dataset.md) for a full description of the dataset structure, pattern naming, and how to visualise trajectories in OVITO.
