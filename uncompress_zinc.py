import os
import tarfile
import gzip
import shutil
import argparse
import multiprocessing
import glob

def extract_file(file_path, output_directory):
    if file_path.endswith('.tgz'):
        with tarfile.open(file_path, 'r:gz') as tar:
            tar.extractall(output_directory)
            print(f"Extracted {file_path} to {output_directory}")
    elif file_path.endswith('.gz'):
        uncompressed_file_path = os.path.join(output_directory, os.path.basename(file_path)[:-3])
        with gzip.open(file_path, 'rb') as f_in:
            with open(uncompressed_file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Extracted {file_path} to {uncompressed_file_path}")

def main(src, num_processes):
    # List all .tgz and .gz files in the directory
    tgz_files = glob.glob(f"{src}/**/*.tgz", recursive=True)
    gz_files = glob.glob(f"{src}/**/*.gz", recursive=True)

    # Combine the lists, ensuring only .tgz and .gz files are included
    files_to_process = tgz_files + gz_files

    # Create a multiprocessing pool with the specified number of processes
    pool = multiprocessing.Pool(processes=num_processes)

    # Define a function that specifies the arguments for extract_file
    def process_file(file):
        return file, src

    # Use the multiprocessing pool to extract files concurrently
    pool.starmap(extract_file, map(process_file, files_to_process))

    print("All .tgz and .gz files have been uncompressed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Uncompress .tgz and .gz files in a directory.")
    parser.add_argument("--src", type=str, help="Path to the directory with .tgz and .gz files")
    parser.add_argument("--num_processes", type=int, default=multiprocessing.cpu_count(), 
                        help="Number of processes to use for parallel processing (default is all available cores)")

    args = parser.parse_args()
    main(args.src, args.num_processes)
