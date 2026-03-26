# Active textiles move in confined areas dataset

This dataset contains molecular dynamics simulation inputs/outputs and physical experiment recordings for bobbin lace patterns studied for their pore contraction behaviour under selective yarn fixing.

---

## Folder structure

```
/
+-- lammps_data/         LAMMPS input files (particle mesh, one per pattern)
+-- simulations/         MD simulation results (one subfolder per pattern)
|   +-- Pattern_{ID}/
|       +-- *.lammpstrj  Trajectories: stabilisation + each contraction step
|       +-- *.log        LAMMPS run log
|       +-- sim1.restart Restart file after stabilisation
|       \-- plots/
|           +-- *_holes.json           Pore geometry time-series
|           +-- *_yarn_lengths.json    Yarn contraction time-series
|           +-- *_analysis.json        Combined analysis output
|           +-- *_relArea.*            Relative pore area change (PDF/PNG)
|           +-- *_area_dist.*          Pore area distribution (GIF/JPG/PNG)
|           +-- *_hist.pdf             Pore area histograms
|           +-- *_lengths_subplots.pdf Yarn length evolution
|           \-- analysis_manifest.txt
\-- textile samples/     Physical experiments: video, photos, and pore analysis
    +-- {ID}-{n} {date}/ One lace sample (ID = pattern, n = replicate number)
    |   +-- *.MOV        Contraction video
    |   +-- DSC_*.JPG    Photographs (before and after contraction)
    |   \-- plots/       Same analysis outputs as simulations/plots/
    +-- knit_{date}/     Reference knitted sample (article label: K)
    +-- weave_{date}/    Reference woven sample (article label: W)
    \-- images_tools/    Photographs of experimental equipment
```

---

## Pattern naming

Two naming systems are used across the files.

The **article labels** organise patterns by structural family (K = knit, W = woven, S = square, D = diamond, H = hexagonal, R = rose, Sp = spruce). The **pattern IDs** in file and folder names follow the numbering of the reference bobbin lace pattern book: Uta Ulrich, _Grunde mit System_ (2011).

| Article label | Article full name | Pattern ID |
|---|---|---|
| K   | Knit      | knit  |
| W   | Woven     | woven |
| S1  | Square-1  | 1092  |
| S2  | Square-2  | 1084  |
| D1  | Diamond-1 | 2001  |
| D2  | Diamond-2 | 2002  |
| D3  | Diamond-3 | 2009  |
| D4  | Diamond-4 | 7002  |
| D5  | Diamond-5 | 7001  |
| D6  | Diamond-6 | 7000  |
| D7  | Diamond-7 | 2023  |
| H1  | Hexa-1    | 2048  |
| R1  | Rose-1    | 3023v2 |
| R2  | Rose-2    | 3025  |
| R3  | Rose-3    | 3021  |
| R4  | Rose-4    | 3024  |
| R5  | Rose-5    | 3026  |
| R6  | Rose-6    | 3013  |
| R7  | Rose-7    | 3018  |
| R8  | Rose-8    | 3001  |
| R9  | Rose-9    | 3003  |
| R10 | Rose-10   | 3034  |
| R11 | Rose-11   | 3043  |
| R12 | Rose-12   | 3051_1t, 3051v2 |
| R13 | Rose-13   | 3052  |
| Sp1 | Spruce-1  | 3059  |

**R12 (Rose-12) variants:** `3051_1t` uses a 4-step contraction sequence `[[1],[3],[1,2],[1,2,3,4,5,6,7,8]]`; `3051v2` uses a 2-step sequence `[[1,2],[1,2,3,4,5,6,7,8]]`. Three physical replicates are available under `textile samples/3051/`.

---

## Simulation file naming

LAMMPS data and trajectory files follow the convention:

```
pattern{ID}_{dist}_{mass}_{threshold}_{ks}_{kb}[_contract_fix_{yarns}]
```

| Field | Value | Description |
|---|---|---|
| `dist`      | 0.25 | Particle spacing (simulation units) |
| `mass`      | 1.0  | Particle mass |
| `threshold` | 0.0  | Contact force threshold |
| `ks`        | 30.0 | Stretching spring constant |
| `kb`        | 0.1  | Bending stiffness |
| `yarns`     | e.g. `1_2_3_4` | Underscore-separated IDs of yarns held fixed during contraction |

Files without the `_contract_fix_*` suffix are the stabilisation runs. Each `_contract_fix_{yarns}` file is a subsequent contraction step.

---

## Experimental samples

Each physical sample folder (`{ID}-{n} {date}/`) contains:
- A video (`*.MOV`) of the contraction experiment
- JPG photographs of the sample (before and after)
- A `plots/` subfolder with pore analysis outputs (same format as simulation analysis)

The `images_tools/` folder contains photos of the experimental tools used for the contraction tests.

---

## Additional notes

- `textile samples/2003-01/` and `textile samples/3110-10/` are physical samples not used in the article simulations.

---

## Software

Simulations were run with [LAMMPS](https://www.lammps.org/). Analysis and figure scripts are available at: [code repository URL]

---

## Visualising simulations in OVITO

The following steps use [OVITO Basic](https://www.ovito.org/) (free version).

1. **Open the input file** - File -> Open -> select the `.data` file for the pattern of interest from `lammps_data/`.

2. **Set the import format** - when the *LAMMPS Data File Import* dialog appears, choose:
   - Atom style: `bond`
   - Column mapping: `atom-ID`, `molecule-ID`, `atom-type`, `x`, `y`, `z`

3. **Adjust display sizes** - in the *Visual Elements* panel:
   - Set *Particle Radius* to `0.5`
   - Set *Bond Width* to `1.0`

4. **Load the trajectory** - in the *Pipeline* panel click *Add modification...* -> *Load Trajectory* -> select the corresponding `.lammpstrj` file from `simulations/Pattern_{ID}/`.

The stabilisation trajectory has no suffix; contraction step trajectories are named `_contract_fix_{yarns}.lammpstrj`. Load them in sequence to follow the full contraction protocol.