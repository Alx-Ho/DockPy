import pandas as pd
import os
import shutil
import argparse
import sys
from pathlib import Path

def copy_top_molecules(input_dir, output_dir, n):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Path to the CSV file
    csv_path = os.path.join(input_dir, "affinity_results.csv")
    
    # Check if CSV file exists
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} does not exist.")
        sys.exit(1)
    
    # Read CSV file
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)
    
    # Sort by affinity (lowest values = best)
    df_sorted = df.sort_values(by='affinity', ascending=True)
    
    # Select top N molecules
    df_top_n = df_sorted.head(n)
    
    # Copy the CSV file to the output directory
    shutil.copy2(csv_path, os.path.join(output_dir, "affinity_results.csv"))
    print(f"Copied affinity_results.csv to {output_dir}")
    
    # Find all PDBQT files in the input directory and its subdirectories
    all_pdbqt_files = {}
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".pdbqt"):
                all_pdbqt_files[file] = os.path.join(root, file)
    
    # Calculate the required number of digits for padding based on N
    padding_digits = len(str(n))
    
    # Process each selected molecule
    for i, (_, row) in enumerate(df_top_n.iterrows(), 1):
        # Extract the original file name from the path
        original_path = row['path']
        original_filename = os.path.basename(original_path)
        
        # Create new file name with dynamic padding based on N
        new_filename = f"rank_{i:0{padding_digits}d}_{original_filename}"
        
        # Find the source file
        if original_filename in all_pdbqt_files:
            source_path = all_pdbqt_files[original_filename]
            dest_path = os.path.join(output_dir, new_filename)
            
            try:
                shutil.copy2(source_path, dest_path)
                print(f"Copied {original_filename} to {new_filename}")
            except Exception as e:
                print(f"Error copying {original_filename}: {e}")
        else:
            print(f"Warning: File {original_filename} not found. Skipping.")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Copy top N molecules with best affinities to an output directory.')
    parser.add_argument('input_dir', type=str, help='Input directory containing affinity_results.csv')
    parser.add_argument('output_dir', type=str, help='Output directory to copy files to')
    parser.add_argument('n', type=int, help='Number of top molecules to select')
    
    args = parser.parse_args()
    
    # Call the function
    copy_top_molecules(args.input_dir, args.output_dir, args.n)
    
    print(f"Successfully copied top {args.n} molecules to {args.output_dir}")

if __name__ == "__main__":
    main()