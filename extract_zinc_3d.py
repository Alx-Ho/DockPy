import os
import shutil
import argparse
import glob

def copy_pdbqt_files(src_dir, output_dir):
    # Using glob.glob to find all .pdbqt files in src_dir and its subdirectories
    for src_path in glob.glob(os.path.join(src_dir, '**', '*.pdbqt'), recursive=True):
        file_name = os.path.basename(src_path)
        dst_path = os.path.join(output_dir, file_name)
        shutil.copy(src_path, dst_path)
        print(f"Copying {src_path} to {dst_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy .pdbqt files from a directory tree to a flat output directory.")
    parser.add_argument("src", help="Source directory where .pdbqt files are located.")
    parser.add_argument("dst", help="Destination directory where .pdbqt files will be copied.")
    args = parser.parse_args()

    if not os.path.exists(args.dst):
        os.makedirs(args.dst)

    copy_pdbqt_files(args.src, args.dst)
