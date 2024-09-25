import argparse
import os
import multiprocessing
from multiprocessing import Pool, Manager
import time
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem

def process_pdbqt(pdbqt_file, dst, output_format, log_queue):
    try:
        id = os.path.splitext(os.path.basename(pdbqt_file))[0]
        output_file = os.path.join(dst, f"{id}.{output_format}")
        
        if os.path.exists(output_file):
            log_queue.put(f"Output file already exists, skipping: {output_file}")
            return

        # Read PDBQT file
        pdbqt_mol = PDBQTMolecule.from_file(pdbqt_file, skip_typing=True)
        rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)

        if not rdkit_mols:
            log_queue.put(f"No valid molecules found in PDBQT file: {pdbqt_file}")
            return

        if output_format == 'pdb':
            # Write only the first conformer of the first molecule to PDB
            with Chem.PDBWriter(output_file) as pdb_writer:
                pdb_writer.write(rdkit_mols[0])
        elif output_format == 'sdf':
            # Write all conformers of all molecules to SDF
            with Chem.SDWriter(output_file) as sdf_writer:
                for mol in rdkit_mols:
                    for i in range(mol.GetNumConformers()):
                        sdf_writer.write(mol, confId=i)

        log_queue.put(f"Successfully wrote to {output_file}")

    except Exception as e:
        log_queue.put(f"An error occurred while processing {pdbqt_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Convert PDBQT files to PDB or SDF format.')
    parser.add_argument('--pdbqt_dir', type=str, help='Input directory containing PDBQT files')
    parser.add_argument('--dst', type=str, help='Output directory for converted files', default='.')
    parser.add_argument('--num_processes', type=int, help='Number of processes to use', default=4)
    parser.add_argument('--file_limit', type=int, help='The maximum number of PDBQT files to process', default=-1)
    parser.add_argument('--output_format', type=str, choices=['pdb', 'sdf'], default='sdf',
                        help='Output format (pdb or sdf). Default is sdf.')
    args = parser.parse_args()

    dst = args.dst
    if not os.path.exists(dst):
        os.makedirs(dst)
        
    if args.num_processes == -1: 
        args.num_processes = multiprocessing.cpu_count()
        
    args.num_processes = int(args.num_processes)

    pdbqt_files = [os.path.join(args.pdbqt_dir, f) for f in os.listdir(args.pdbqt_dir) if f.endswith('.pdbqt')]
    total_files = len(pdbqt_files)
    processed_files = 0

    with Manager() as manager:
        log_queue = manager.Queue()
        with Pool(processes=args.num_processes) as pool:
            pool.starmap_async(process_pdbqt, [(pdbqt_file, dst, args.output_format, log_queue) for pdbqt_file in pdbqt_files])

            while total_files > 0:
                while not log_queue.empty():
                    log_message = log_queue.get()
                    print(log_message)
                    processed_files += 1

                    if args.file_limit != -1 and processed_files >= args.file_limit:
                        pool.terminate()
                        print(f"Processed file limit of {args.file_limit} reached. Stopping further processing.")
                        return

                print(f"Remaining PDBQT files: {total_files - processed_files}")
                if all(worker.is_alive() == False for worker in pool._pool) or (total_files - processed_files) == 0:
                    break
                time.sleep(1)

    print("Conversion complete.")

if __name__ == "__main__":
    main()
