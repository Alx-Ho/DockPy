import argparse
import os
import glob
from vina import Vina
import subprocess

# Set up argument parser
parser = argparse.ArgumentParser(description='Run molecular docking using AutoDock Vina.')
parser.add_argument('--ligands', required=True, help='Ligand file or directory with ligand files (pdbqt format)')
parser.add_argument('--receptor', required=True, help='Receptor file (pdbqt format)')
parser.add_argument('--center', required=True, nargs=3, type=float, help='Center of the docking box (x y z)')
parser.add_argument('--box_size', required=True, nargs=3, type=int, help='Size of the docking box (x y z)')
parser.add_argument('--exhaustiveness', type=int, default=32, help='Exhaustiveness of the search (optional, default=32)')
parser.add_argument('--n_poses', type=int, default=20, help='Number of poses to generate (optional, default=20)')
parser.add_argument('--out_poses', type=int, default=1, help='Number of poses to write out (optional, default=1)')
parser.add_argument('--output_dir', default='.', help='Directory to write output files (optional, default is current directory)')
parser.add_argument('--overwrite', type=bool, default=False, help='Allow overwriting of poses (optional, default is False)')


# Parse arguments
args = parser.parse_args()

# Function to create a base file name from the input file
def create_base_filename(file_path):
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    return base_name

# Determine if the ligands argument is a directory or a single file
ligand_paths = []
if os.path.isdir(args.ligands):
    ligand_paths = glob.glob(os.path.join(args.ligands, '*.pdbqt'))
elif os.path.isfile(args.ligands):
    ligand_paths.append(args.ligands)
else:
    raise ValueError("The ligands argument must be either a directory or a file path.")

# Create the output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# Initialize Vina
v = Vina(sf_name='vina')
v.set_receptor(args.receptor)
v.compute_vina_maps(center=args.center, box_size=args.box_size)

# Log file to record arguments for each ligand
log_file_path = os.path.join(args.output_dir, 'docking_log.txt')
with open(log_file_path, 'w') as log_file:
    log_file.write("\nDocking Log\n===========\n")

# Process each ligand
total_ligands = len(ligand_paths)
for i, ligand in enumerate(ligand_paths, start=1):
    try: 
        ligand_base = create_base_filename(ligand)
        receptor_base = create_base_filename(args.receptor)
        output_minimized = os.path.join(args.output_dir, f'{ligand_base}_{receptor_base}_minimized.pdbqt')
        output_vina = os.path.join(args.output_dir, f'{ligand_base}_{receptor_base}_vina_out.pdbqt')

        print(f'Processing {ligand_base} ({i}/{total_ligands})...')
        with open(log_file_path, 'a') as log_file:
            log_file.write(f'\nProcessing {ligand_base} ({i}/{total_ligands})...\n')

        v.set_ligand_from_file(ligand)
        
        # Write arguments to log file
        with open(log_file_path, 'a') as log_file:
            log_file.write(f'Ligand: {ligand_base}, Receptor: {receptor_base}, Center: {args.center}, '
                        f'Box Size: {args.box_size}, Exhaustiveness: {args.exhaustiveness}, '
                        f'Number of Poses: {args.n_poses}, Output Poses: {args.out_poses}\n')

        # Score the current pose
        energy = v.score()

        print(f'Score before minimization for {ligand_base}: %.3f (kcal/mol)' % energy[0])
        with open(log_file_path, 'a') as log_file:
            log_file.write(f'\nScore before minimization for {ligand_base}: %.3f (kcal/mol)' % energy[0])

        # Minimize locally the current pose
        energy_minimized = v.optimize()

        print(f'Score after minimization for {ligand_base}: %.3f (kcal/mol)' % energy_minimized[0])
        with open(log_file_path, 'a') as log_file:
            log_file.write(f'\nScore after minimization for {ligand_base}: %.3f (kcal/mol)' % energy_minimized[0])

        v.write_pose(output_minimized, overwrite=args.overwrite)

        print(f'Wrote {output_minimized}')
        with open(log_file_path, 'a') as log_file:
            log_file.write(f'\nWrote {output_minimized}')

        # Dock the ligand
        v.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_poses)
        v.write_poses(output_vina, n_poses=args.out_poses, overwrite=args.overwrite)

        print(f'Wrote {output_vina}')
        with open(log_file_path, 'a') as log_file:
            log_file.write(f'\nWrote {output_vina}')

        print(f'Completed processing {ligand_base}. {total_ligands - i} ligand(s) remaining.')
        with open(log_file_path, 'a') as log_file:
            log_file.write(f'\nCompleted processing {ligand_base}. {total_ligands - i} ligand(s) remaining.')


        subprocess.run(['python', 'pdbqt_extract_zincid_affinity.py', '--src', args.output_dir, '--output_csv', f'{args.output_dir}/affinity_results.csv'], check=True)
        print(f'Updated {args.output_dir}/affinity_results.csv')
        with open(log_file_path, 'a') as log_file:
            log_file.write(f'\nUpdated {args.output_dir}/affinity_results.csv')

    except RuntimeError as e:
        error_message = str(e)
        print(f"Error processing {ligand_base}: {error_message}")
        with open(log_file_path, 'a') as log_file:
            log_file.write(f"\nError processing {ligand_base}: {error_message}\n")
        continue

# Indicate completion in the log file
with open(log_file_path, 'a') as log_file:
    log_file.write

subprocess.run(['python', 'pdbqt_extract_zincid_affinity.py', '--src', args.output_dir, '--output_csv', f'{args.output_dir}/affinity_results.csv'], check=True)