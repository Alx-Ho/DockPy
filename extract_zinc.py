import os
import shutil
import argparse
import glob

def copy_files(src_dir, output_dir, extensions):
    for ext in extensions:
        # Using glob.glob to find all files with the given extension in src_dir and its subdirectories
        for src_path in glob.glob(os.path.join(src_dir, '**', f'*.{ext}'), recursive=True):
            file_name = os.path.basename(src_path)
            dst_path = os.path.join(output_dir, file_name)
            shutil.copy(src_path, dst_path)
            print(f"Copying {src_path} to {dst_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy files with specific extensions from a directory tree to a flat output directory.")
    parser.add_argument("src", help="Source directory where files are located.")
    parser.add_argument("dst", help="Destination directory where files will be copied.")
    args = parser.parse_args()

    if not os.path.exists(args.dst):
        os.makedirs(args.dst)

    # List of file extensions to copy
    extensions = ['pdbqt', 'smi']
    copy_files(args.src, args.dst, extensions)
