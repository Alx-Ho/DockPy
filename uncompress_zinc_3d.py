import os
import tarfile
import argparse
import multiprocessing
import glob

def extract_tgz_file(file_path, output_directory):
    with tarfile.open(file_path, 'r:gz') as tar:
        tar.extractall(output_directory)
        print(f"Extracted {file_path} to {output_directory}")

def main(src, num_processes):
    # List all files in the directory
    files = glob.glob(f"{src}/**/*.tgz")

    # Create a multiprocessing pool with the specified number of processes
    pool = multiprocessing.Pool(processes=num_processes)

    # Define a function that specifies the arguments for extract_tgz_file
    def process_file(file):
        return file, src

    # Use the multiprocessing pool to extract .tgz files concurrently
    pool.starmap(extract_tgz_file, map(process_file, [file for file in files if file.endswith('.tgz')]))

    print("All .tgz files have been uncompressed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Uncompress .tgz files in a directory.")
    parser.add_argument("--src", type=str, help="Path to the directory with .tgz files")
    parser.add_argument("--num_processes", type=int, default=multiprocessing.cpu_count(), 
                        help="Number of processes to use for parallel processing (default is all available cores)")

    args = parser.parse_args()
    main(args.src, args.num_processes)
