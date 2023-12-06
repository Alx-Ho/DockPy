import os
import subprocess
import argparse

def main(curl_file_path, destination):

    # Ensure the destination directory exists
    os.makedirs(destination, exist_ok=True)

    # Create subdirectories in the destination
    downloaded_dir = os.path.join(destination, 'downloaded_zinc')
    extracted_dir = os.path.join(destination, 'extracted_pdbqts')
    os.makedirs(downloaded_dir, exist_ok=True)
    os.makedirs(extracted_dir, exist_ok=True)

    # Step 1: Download
    subprocess.run(['python', 'download_zinc.py', '--curl_file', curl_file_path, '--dst', downloaded_dir], check=True)

    # Step 2: Uncompress
    subprocess.run(['python', 'uncompress_zinc.py', '--src', downloaded_dir], check=True)

    # Step 3: Extract
    subprocess.run(['python', 'extract_zinc.py', downloaded_dir, extracted_dir], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated script to process ZINC 3D files.")
    parser.add_argument("--curl_file", type=str, help="Path to the .curl file")
    parser.add_argument("--dst", type=str, help="Path to the destination directory")

    args = parser.parse_args()
    main(args.curl_file, args.dst)
