import argparse
from meeko import MoleculePreparation, PDBQTWriterLegacy
import rdkit
import os
import multiprocessing
from multiprocessing import Pool, Manager
import time

# Function to process each SMILES line
def process_smiles(line, dst, log_queue):
    try:
        parts = line.strip().split()
        if len(parts) != 2:
            log_queue.put(f"Skipping invalid line: {line}")
            return
        smiles, id = parts

        output_file = os.path.join(dst, f"{id}.pdbqt")
        if os.path.exists(output_file):
            log_queue.put(f"Output file already exists, skipping: {output_file}")
            return

        lig = rdkit.Chem.MolFromSmiles(smiles)
        if lig is None:
            log_queue.put(f"Invalid SMILES string for ID {id}: {smiles}")
            return

        protonated_lig = rdkit.Chem.AddHs(lig)
        rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(protonated_lig)

        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                with open(output_file, 'w') as f:
                    f.write(pdbqt_string)
                log_queue.put(f"Successfully wrote to {output_file}")
            else:
                log_queue.put(f"Error in writing to PDBQT for ID {id}: {error_msg}")

    except Exception as e:
        log_queue.put(f"An error occurred while processing ID {id}: {e}")


def main():
    parser = argparse.ArgumentParser(description='Convert SMILES to PDBQT for molecular docking.')
    parser.add_argument('--smiles_file', type=str, help='Input file with SMILES strings and IDs')
    parser.add_argument('--dst', type=str, help='Output directory for PDBQT files', default='.')
    parser.add_argument('--num_processes', type=int, help='Number of processes to use', default=multiprocessing.cpu_count())
    args = parser.parse_args()

    dst = args.dst
    if not os.path.exists(dst):
        os.makedirs(dst)

    with open(args.smiles_file, 'r') as file:
        lines = file.readlines()

    total_smiles = len(lines)

    with Manager() as manager:
        log_queue = manager.Queue()
        with Pool(processes=args.num_processes) as pool:
            pool.starmap_async(process_smiles, [(line, dst, log_queue) for line in lines])

            # Continuously check the queue for new log messages
            while total_smiles > 0:
                while not log_queue.empty():
                    log_message = log_queue.get()
                    print(log_message)
                    total_smiles -= 1  # Decrement the counter
                print(f"Remaining SMILES strings: {total_smiles}")
                if all(worker.is_alive() == False for worker in pool._pool):
                    break
                time.sleep(1)  # Wait a bit before checking the queue again

    print("Conversion complete.")


if __name__ == "__main__":
    main()