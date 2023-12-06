import os
import subprocess
import argparse
import glob

def main(curl_file_path, destination, smiles_limit):

    # Ensure the destination directory exists
    os.makedirs(destination, exist_ok=True)

    # Create subdirectories in the destination
    downloaded_dir = os.path.join(destination, 'downloaded_zinc')
    extracted_dir = os.path.join(destination, 'extracted_pdbqts')
    os.makedirs(downloaded_dir, exist_ok=True)
    os.makedirs(extracted_dir, exist_ok=True)

    # Download
    subprocess.run(['python', 'utils/download_zinc.py', '--curl_file', curl_file_path, '--dst', downloaded_dir], check=True)

    # Uncompress
    subprocess.run(['python', 'utils/uncompress_zinc.py', '--src', downloaded_dir], check=True)

    # Format SMILES
    for smiles_file in glob.glob(f'{downloaded_dir}/**/*.smi', recursive=True): 
        print(f'Converting {smiles_file} to PDBQT')
        subprocess.run(['python', 'utils/smiles_to_pdbqt.py', '--smiles_file', smiles_file, '--dst', os.path.join(downloaded_dir, 'converted_smiles'), '--smiles_limit', str(smiles_limit)], check=True)

    # Extract
    subprocess.run(['python', 'utils/extract_zinc.py', downloaded_dir, extracted_dir], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated script to process ZINC curl files.")
    parser.add_argument("--curl_file", type=str, help="Path to the .curl file")
    parser.add_argument("--dst", type=str, help="Path to the destination directory")
    parser.add_argument("--smiles_limit", type=int, help="The maximum number of SMILES strings to process. NOTE: This limit is approximate since conversions are done with multiprocessing.", default=None)

    args = parser.parse_args()
    main(args.curl_file, args.dst, args.smiles_limit)
