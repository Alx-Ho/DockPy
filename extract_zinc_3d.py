import os
import shutil
import argparse

def copy_pdbqt_files(src_dir, output_dir):
    for root, _, files in os.walk(src_dir):
        for file in files:
            if file.endswith(".pdbqt"):
                src_path = os.path.join(root, file)
                dst_path = os.path.join(output_dir, file)
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
