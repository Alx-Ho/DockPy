import argparse
import os
import glob
from vina import Vina
import subprocess
import logging
import time

# Function to create a base file name from the input file
def create_base_filename(file_path):
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    return base_name

def setup_logging(output_dir):
    log_file_path = os.path.join(output_dir, 'docking_log.log')

    # Set up logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Create file handler which logs even debug messages
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.INFO)

    # Create console handler with a higher log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    # Create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

def main(): 
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

    # Set up logging in the specified output directory
    setup_logging(args.output_dir)

    # Initialize Vina
    v = Vina(sf_name='vina')
    v.set_receptor(args.receptor)
    v.compute_vina_maps(center=args.center, box_size=args.box_size)

    # Process each ligand
    total_ligands = len(ligand_paths)
    start_time = time.time()

    # Initialize total time for processing ligands
    total_time_for_ligands = 0 

    for i, ligand in enumerate(ligand_paths, start=1):
        ligand_start_time = time.time()

        try: 
            ligand_base = create_base_filename(ligand)
            receptor_base = create_base_filename(args.receptor)
            output_minimized = os.path.join(args.output_dir, f'{ligand_base}_{receptor_base}_minimized.pdbqt')
            output_vina = os.path.join(args.output_dir, f'{ligand_base}_{receptor_base}_vina_out.pdbqt')

            logging.info(f'Processing {ligand_base} ({i}/{total_ligands})...')

            v.set_ligand_from_file(ligand)
                
            # Log the arguments
            logging.info(f'Ligand: {ligand_base}, Receptor: {receptor_base}, Center: {args.center}, '
                        f'Box Size: {args.box_size}, Exhaustiveness: {args.exhaustiveness}, '
                        f'Number of Poses: {args.n_poses}, Output Poses: {args.out_poses}')

            # Score the current pose
            energy = v.score()

            logging.info(f'Score before minimization for {ligand_base}: {energy[0]:.3f} kcal/mol')

            # Minimize locally the current pose
            energy_minimized = v.optimize()

            logging.info(f'Score after minimization for {ligand_base}: {energy_minimized[0]:.3f} kcal/mol')

            v.write_pose(output_minimized, overwrite=args.overwrite)

            logging.info(f'Wrote {output_minimized}')

            # Dock the ligand
            v.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_poses)
            v.write_poses(output_vina, n_poses=args.out_poses, overwrite=args.overwrite)
            logging.info(f'Wrote {output_vina}')

            ligand_end_time = time.time()
            time_taken_for_ligand = ligand_end_time - ligand_start_time
            total_time_for_ligands += time_taken_for_ligand

            average_time_per_ligand = total_time_for_ligands / i
            ligands_remaining = total_ligands - i
            estimated_time_left = average_time_per_ligand * ligands_remaining

            # Convert estimated_time_left to hours, minutes, and seconds
            hours, remainder = divmod(estimated_time_left, 3600)
            minutes, seconds = divmod(remainder, 60)

            logging.info(f'Completed processing {ligand_base}. {total_ligands - i} ligand(s) remaining.')
            logging.info(f'Time taken for {ligand_base}: {ligand_end_time - ligand_start_time} seconds.')
            logging.info(f'Estimated time left: {int(hours)} hours, {int(minutes)} minutes, {seconds:.2f} seconds')

            subprocess.run(['python', 'pdbqt_extract_zincid_affinity.py', '--src', args.output_dir, '--output_csv', f'{args.output_dir}/affinity_results.csv'], check=True)
            logging.info(f'Updated {args.output_dir}/affinity_results.csv')

        except RuntimeError as e:
            error_message = str(e)
            logging.error(f"Error processing {ligand_base}: {error_message}")
            continue

    # Indicate completion in the log file
    end_time = time.time()
    logging.info(f'Total docking time: {end_time - start_time} seconds.')

    subprocess.run(['python', 'pdbqt_extract_zincid_affinity.py', '--src', args.output_dir, '--output_csv', f'{args.output_dir}/affinity_results.csv'], check=True)