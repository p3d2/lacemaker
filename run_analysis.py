import os
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process lammpstrj files.')
    parser.add_argument('folder_path', help='Input folder path')
    args = parser.parse_args()

    folder_path = args.folder_path

    # Check if 'plots' folder exists; if not, create it
    plots_folder = os.path.join(folder_path, 'plots')
    if not os.path.exists(plots_folder):
        os.makedirs(plots_folder)
        print(f"Created 'plots' folder at: {plots_folder}")

    # Iterate over all files ending with 'lammpstrj'
    for filename in os.listdir(folder_path):
        if filename.endswith('lammpstrj'):
            # Extract the base name without extension
            name = os.path.splitext(filename)[0]

            # Define the expected output file path
            output_file = os.path.join(plots_folder, f'{name}_holes.json')
            # Check if the output file already exists
            if os.path.exists(output_file):
                print(f"Output file {output_file} already exists. Skipping.")
                
            # Build the command
            cmd = ['srun','--mem=4G','python','sim_holes_extract.py', folder_path, name, '--roi', '1.0', '127.0', '1.0', '127.0']
            # Execute the command in the 'plots' directory
            print(f"Running command: {' '.join(cmd)}")
            subprocess.run(cmd)

            cmd2 = ['srun','--mem=4G','python','holes_analysis.py', folder_path, name]
            print(f"Running command: {' '.join(cmd2)}")
            subprocess.run(cmd2)
if __name__ == '__main__':
    main()