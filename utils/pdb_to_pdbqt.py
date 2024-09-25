import argparse
from meeko import MoleculePreparation, PDBQTWriterLegacy
import rdkit
from rdkit import Chem
import os
import multiprocessing
from multiprocessing import Pool, Manager
import time

# Function to process each PDB file
def process_pdb(pdb_file, dst, log_queue):
    try:
        id = os.path.splitext(os.path.basename(pdb_file))[0]
        output_file = os.path.join(dst, f"{id}.pdbqt")
        
        if os.path.exists(output_file):
            log_queue.put(f"Output file already exists, skipping: {output_file}")
            return

        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            log_queue.put(f"Invalid PDB file: {pdb_file}")
            return

        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                with open(output_file, 'w') as f:
                    f.write(pdbqt_string)
                log_queue.put(f"Successfully wrote to {output_file}")
            else:
                log_queue.put(f"Error in writing to PDBQT for {pdb_file}: {error_msg}")

    except Exception as e:
        log_queue.put(f"An error occurred while processing {pdb_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Convert PDB files to PDBQT for molecular docking.')
    parser.add_argument('--pdb_dir', type=str, help='Input directory containing PDB files')
    parser.add_argument('--dst', type=str, help='Output directory for PDBQT files', default='.')
    parser.add_argument('--num_processes', type=int, help='Number of processes to use', default=4)
    parser.add_argument('--file_limit', type=int, help='The maximum number of PDB files to process', default=-1)
    args = parser.parse_args()

    dst = args.dst
    if not os.path.exists(dst):
        os.makedirs(dst)
        
    if args.num_processes == -1: 
        args.num_processes = multiprocessing.cpu_count()
        
    args.num_processes = int(args.num_processes)

    pdb_files = [os.path.join(args.pdb_dir, f) for f in os.listdir(args.pdb_dir) if f.endswith('.pdb')]
    total_files = len(pdb_files)
    processed_files = 0

    with Manager() as manager:
        log_queue = manager.Queue()
        with Pool(processes=args.num_processes) as pool:
            pool.starmap_async(process_pdb, [(pdb_file, dst, log_queue) for pdb_file in pdb_files])

            while total_files > 0:
                while not log_queue.empty():
                    log_message = log_queue.get()
                    print(log_message)
                    processed_files += 1

                    if args.file_limit != -1 and processed_files >= args.file_limit:
                        pool.terminate()
                        print(f"Processed file limit of {args.file_limit} reached. Stopping further processing.")
                        return

                print(f"Remaining PDB files: {total_files - processed_files}")
                if all(worker.is_alive() == False for worker in pool._pool) or (total_files - processed_files) == 0:
                    break
                time.sleep(1)

    print("Conversion complete.")

if __name__ == "__main__":
    main()
