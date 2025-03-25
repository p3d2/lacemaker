import subprocess
import os
import pandas as pd
import json
import sys

def process_jobs(tsv_filename):
    jobs = pd.read_csv(tsv_filename, sep='\t')

    for row in jobs.itertuples(index=False):
        
        if pd.isna(row.pattern_id):
            continue  # Skip this row

        pattern_id = str(row.pattern_id)
        dist_particles = str(float(row.dist_particles))
        units = str(float(row.units))
        mass = str(float(row.mass))
        threshold = str(float(row.threshold))
        ks = str(float(row.ks))
        kb = str(float(row.kb))
        status = int(row.status)
        if status == 1:
            continue # skip pattern

        filename = 'pattern' + pattern_id + '_' + dist_particles + '_' + mass + '_' + threshold + "_" + ks + "_" + kb
        input_data = os.path.join('output', 'lammps_data', filename + '.data')

        # Check if the output file already exists
        if os.path.exists(input_data):
            print(f"Lammps input file {pattern_id} already exists. Skipping mesh generation.")
        else:
            # JOB 1: Generate mesh
            command = (
    "module load scicomp-python-env/ && "
    f"srun python lace_maker.py input/json_patterns/pattern{pattern_id}.json "
    f"--dist_particles={dist_particles} "
    f"--units={units} "
    f"--mass={mass} "
    f"--threshold={threshold} "
    f"--ks1={ks} "
    f"--kb={kb}"
)
            subprocess.run(command, shell=True, executable='/bin/bash')
        
        folder_pattern = os.path.join('output', 'simulations', f'Pattern_{pattern_id}')
        if not os.path.exists(folder_pattern):
            os.makedirs(folder_pattern)

        # JOB 2: Simulation 0: stabilization       
        visc = float(row.visc)
        tstep = float(row.tstep)
        thermodump = int(row.thermodump)
        nframes = int(row.nframes)
        xmin = float(row.xmin)
        xmax = float(row.xmax)
        ymin = float(row.ymin)
        ymax = float(row.ymax)

        bond0_1 = float(row.bond0_1)
        bond0_2 = float(row.bond0_2)
        bond0_3 = float(row.bond0_3)
        bond0_4 = float(row.bond0_4)    
        bond0_5 = float(row.bond0_5)
        bond0_6 = float(row.bond0_6)
        bond0_7 = float(row.bond0_7)
        bond0_8 = float(row.bond0_8)

        b1 = int(bond0_1 != 0)
        b2 = int(bond0_2 != 0)
        b3 = int(bond0_3 != 0)
        b4 = int(bond0_4 != 0)
        b5 = int(bond0_5 != 0)
        b6 = int(bond0_6 != 0)
        b7 = int(bond0_7 != 0)
        b8 = int(bond0_8 != 0)

        sim0_data = os.path.join(folder_pattern, filename + '.lammpstrj')
        sim0_log = os.path.join(folder_pattern, filename + '.log')
        sim0_lmp_path = os.path.join('input', 'lammps_input', 'in_jobs8.lmp')
        restart_path = os.path.join(folder_pattern, 'sim1.restart')
        if os.path.exists(sim0_data):
            print(f"Stabilization simulation {pattern_id} already exists. Skipping simulation 1.")
        else:
            sbatch_script_content = f"""#!/bin/bash
#SBATCH --job-name={pattern_id}_0
#SBATCH --time=00:20:00
#SBATCH --mem=2GB
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=8
#SBATCH --output=log.out

module load gcc/12.3.0 openmpi/4.1.6 openblas fftw eigen
export PATH=$PATH:$HOME/lammps/bin
export OMP_NUM_THREADS=1

# Run simulation
srun lmp -var sim0_data {sim0_data} -var folderdata {input_data} -var sim0_log {sim0_log} -var restart {restart_path}\
    -var visc {visc} -var tstep {tstep} \
    -var thermodump {thermodump} -var nframes {nframes} -var ks {ks} \
    -var bond1 {bond0_1} -var bond2 {bond0_2} -var bond3 {bond0_3} -var bond4 {bond0_4} \
    -var bond5 {bond0_5} -var bond6 {bond0_6} -var bond7 {bond0_7} -var bond8 {bond0_8} \
    -var b1 {b1} -var b2 {b2} -var b3 {b3} -var b4 {b4} \
    -var b5 {b5} -var b6 {b6} -var b7 {b7} -var b8 {b8} \
    -var xmin {xmin} -var xmax {xmax} -var ymin {ymin} -var ymax {ymax} \
    -in {sim0_lmp_path}
"""

            script_filename = f'sbatch_script.sh'

            # Write the SBATCH script to a file
            try:
                with open(script_filename, 'w') as file:
                    file.write(sbatch_script_content)
            except IOError as e:
                print(f"Failed to write SBATCH script for job {pattern_id}: {e}")
                return False
            
            # Submit the job using sbatch
            try:
                result = subprocess.run(['sbatch', script_filename], capture_output=True, text=True, check=True)
                # Optionally, you can parse and store the job ID from result.stdout
            except subprocess.CalledProcessError as e:
                print(f"Failed to submit job {pattern_id}: {e.stderr.strip()}")
                return False
            finally:
                # Clean up the script file
                os.remove(script_filename)

        # JOB 3: Contraction of different yarns

        if not os.path.exists(sim0_data):
            print(f'Run script again after simulation 0 of pattern {pattern_id} is done')
            continue

        visc2 = float(row.visc2)
        tstep2 = float(row.tstep2)
        thermodump2 = int(row.thermodump2)
        nframes2 = int(row.nframes2)
        contract = float(row.contract)
        exp = json.loads(row.exp)
        folder_restart = os.path.join(folder_pattern, 'sim1.restart')

        for k in range(len(exp)):
            
            v_exp = [int(i) for i in exp[k]]
            v_str = [str(i) for i in v_exp]
            v_joined = "_".join(v_str)

            b1 = 1 if 1 in exp[k] else 0
            b2 = 1 if 2 in exp[k] else 0
            b3 = 1 if 3 in exp[k] else 0
            b4 = 1 if 4 in exp[k] else 0
            b5 = 1 if 5 in exp[k] else 0
            b6 = 1 if 6 in exp[k] else 0
            b7 = 1 if 7 in exp[k] else 0
            b8 = 1 if 8 in exp[k] else 0
            sim1_data = os.path.join(folder_pattern, filename + f"_contract_fix_{v_joined}.lammpstrj")
            sim1_log = os.path.join(folder_pattern, filename + f"_contract_fix_{v_joined}.log")
            sim1_lmp_path = os.path.join('input', 'lammps_input', 'in_contract_jobs8.lmp')
            if os.path.exists(sim1_data):
                print(f"Contraction simulation {pattern_id} for bonds {v_str} already exists. Skipping simulation.")
            else:
                sbatch_script_content = f"""#!/bin/bash
#SBATCH --job-name={pattern_id}_{k}
#SBATCH --time=04:00:00
#SBATCH --mem=2GB
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=8
#SBATCH --output=log.out

module load gcc/12.3.0 openmpi/4.1.6 openblas fftw eigen
export PATH=$PATH:$HOME/lammps/bin
export OMP_NUM_THREADS=1

# Run simulation
srun lmp -var sim1_data {sim1_data} -var folderrestart {folder_restart} -var sim1_log {sim1_log} -var self {sim1_lmp_path} \
    -var visc {visc2} -var tstep {tstep2} \
    -var thermodump {thermodump2} -var nframes {nframes2} -var ks {ks} -var contract {contract} \
    -var bond1 {bond0_1} -var bond2 {bond0_2} -var bond3 {bond0_3} -var bond4 {bond0_4} \
    -var bond5 {bond0_5} -var bond6 {bond0_6} -var bond7 {bond0_7} -var bond8 {bond0_8} \
    -var b1 {b1} -var b2 {b2} -var b3 {b3} -var b4 {b4} \
    -var b5 {b5} -var b6 {b6} -var b7 {b7} -var b8 {b8} \
    -var xmin {xmin} -var xmax {xmax} -var ymin {ymin} -var ymax {ymax} \
    -in {sim1_lmp_path}

"""

                script_filename = f'sbatch_script.sh'
                
                # Write the SBATCH script to a file
                try:
                    with open(script_filename, 'w') as file:
                        file.write(sbatch_script_content)
                except IOError as e:
                    print(f"Failed to write SBATCH script for job {pattern_id}: {e}")
                    return False
                
                # Submit the job using sbatch
                try:
                    result = subprocess.run(['sbatch', script_filename], capture_output=True, text=True, check=True)
                    # Optionally, you can parse and store the job ID from result.stdout
                except subprocess.CalledProcessError as e:
                    print(f"Failed to submit job {pattern_id}: {e.stderr.strip()}")
                    return False
                finally:
                    # Clean up the script file
                    os.remove(script_filename)

        # JOB 4: Analysis of *.lammpstrj
        # -------------------------------------
        lammpstrj_files = [f for f in os.listdir(folder_pattern) if f.endswith('.lammpstrj')]
        for lammpstrj_file in lammpstrj_files:
            base_name = lammpstrj_file[:-10]  # remove the ".lammpstrj" extension

            plots_folder = os.path.join(folder_pattern, 'plots')
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)

            # If the GIF already exists, skip processing
            gif_file = base_name + '_area_dist.gif'
            gif_path = os.path.join(plots_folder, gif_file)
            if os.path.exists(gif_path):
                print(f"Analysis for {lammpstrj_file} already done. Skipping.")
                continue

            # Create SBATCH script for analyzing this lammpstrj
            sbatch_script_content = f"""#!/bin/bash
#SBATCH --job-name=an_{pattern_id}
#SBATCH --time=00:10:00
#SBATCH --mem=4GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=analysis.log

module load scicomp-python-env
python sim_holes_extract.py {folder_pattern} {base_name} --roi 1.0 127.0 1.0 127.0
"""
            script_filename = 'sbatch_analyze.sh'
            try:
                with open(script_filename, 'w') as file:
                    file.write(sbatch_script_content)
            except IOError as e:
                print(f"Failed to write analysis SBATCH script for {pattern_id}: {e}")
                return False

            try:
                result = subprocess.run(['sbatch', script_filename],
                                        capture_output=True, text=True, check=True)
                # Optionally parse job ID from result.stdout
            except subprocess.CalledProcessError as e:
                print(f"Failed to submit analysis job for {pattern_id}: {e.stderr.strip()}")
                return False
            finally:
                os.remove(script_filename)


if __name__ == "__main__":
    tsv_filename = os.path.join('assets', 'data', 'jobs.tsv')
    process_jobs(tsv_filename)