import argparse
import os
import pandas as pd
from multiprocessing import Pool

def extract_zinc_id_and_affinity(file_path):
    zinc_id = None
    affinity = None

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('REMARK  Name ='):
                zinc_id = line.split('=')[-1].strip()
            elif line.startswith('REMARK VINA RESULT:'):
                affinity = line.split()[3]

    return zinc_id, affinity, file_path

def process_file(file_path):
    zinc_id, affinity, path = extract_zinc_id_and_affinity(file_path)
    return {'zinc_id': zinc_id, 'affinity': affinity, 'path': path}

def process_directory(directory_path):
    files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('vina_out.pdbqt')]
    
    with Pool() as pool:
        results = pool.map(process_file, files)

    return results

def main():
    parser = argparse.ArgumentParser(description='Process PDBQT files and extract ZINC ID and affinity.')
    parser.add_argument('--src', type=str, help='Directory containing PDBQT files')
    parser.add_argument('--output_csv', type=str, help='Output CSV file name')

    args = parser.parse_args()

    data = process_directory(args.src)
    df = pd.DataFrame(data).sort_values(by='affinity')
    df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()
