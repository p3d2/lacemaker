#!/bin/bash
#SBATCH --job-name=Lammps_simulations
#SBATCH --time=00:05:00
#SBATCH --mem=8GB
#SBATCH --nodes=1
#SBATCH --ntasks=8 # Ensure there are 4 tasks to match the number of CPUs requested
#SBATCH --cpus-per-task=1 # Set ths to 1 if you're using ntasks to specify processor count
#SBATCH --output=log.out

module load gcc/12.3.0 openmpi/4.1.6 openblas fftw eigen

export PATH=$PATH:$PWD/lammps/bin
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run simulation
srun lmp -var input_file "$WRKDIR/lacemaker/output/lammps_data/pattern1024_0.5_1.0_0.0_30.0_5.0.data" -var fname "1024_2408141144" \
    -var modify_bond1 1 -var modify_bond2 1 -var modify_bond3 1 -var modify_bond4 1 \
    -var modify_bond5 1 -var modify_bond6 1 -var modify_bond7 0 -var modify_bond8 0 \
    -in "$WRKDIR/lacemaker/sim_run.lmp"