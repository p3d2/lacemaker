#!/bin/bash
#SBATCH --job-name=Lammps_simulations
#SBATCH --time=00:05:00
#SBATCH --mem=8GB
#SBATCH --nodes=1
#SBATCH --ntasks=8 # Ensure there are 4 tasks to match the number of CPUs requested
#SBATCH --cpus-per-task=1 # Set ths to 1 if you're using ntasks to specify processor count
#SBATCH --output=log.out

module load openmpi/4.1.6 openblas fftw eigen
export PATH=$PATH:$PWD/lammps/bin
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run simulation
srun lmp -in "$WRKDIR/lacemaker/sim_run.lmp"
