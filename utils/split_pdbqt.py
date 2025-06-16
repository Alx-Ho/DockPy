import os
import sys
import argparse
import csv
import re

def split_pdbqt_models(input_file, output_dir=None, create_csv=True):
    """
    Splits a PDBQT file containing multiple models into individual PDBQT files
    and extracts the docking scores into a CSV file.
    
    Args:
        input_file (str): Path to the input PDBQT file
        output_dir (str, optional): Directory to save the output files. 
                                   If None, uses the same directory as the input file.
        create_csv (bool): Whether to create a CSV file with docking scores
    
    Returns:
        list: List of paths to the created output files
    """
    # Set output directory
    if output_dir is None:
        output_dir = os.path.dirname(input_file)
    
    # Create output directory if it doesn't exist
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Get base filename without extension
    base_filename = os.path.basename(input_file)
    base_name, ext = os.path.splitext(base_filename)
    
    # Read the input file
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Split the content based on ENDMDL
    models = content.split('ENDMDL')
    
    # The last element will be empty or contain trailing content, remove it if empty
    if models[-1].strip() == '':
        models = models[:-1]
    
    output_files = []
    scores = []
    
    # Regular expression to extract the docking score
    score_pattern = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)")
    
    # Process each model
    for i, model in enumerate(models):
        # Add ENDMDL back since it was removed by the split
        if not model.strip().endswith('ENDMDL'):
            model += '\nENDMDL'
        
        # Create output filename with zero-padded model number
        output_filename = f"{base_name}_model_{i+1:02d}{ext}"
        output_path = os.path.join(output_dir, output_filename)
        
        # Extract the docking score
        score_match = score_pattern.search(model)
        score = float(score_match.group(1)) if score_match else None
        scores.append((i+1, score))
        
        # Write the model to a new file
        with open(output_path, 'w') as f:
            f.write(model.strip() + '\n')
        
        output_files.append(output_path)
        print(f"Created: {output_path} (Score: {score})")
    
    # Create CSV file with docking scores
    if create_csv and scores:
        csv_filename = os.path.join(output_dir, f"{base_name}_scores.csv")
        with open(csv_filename, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['Model', 'Docking Score'])
            for model_num, score in scores:
                csv_writer.writerow([f"Model {model_num}", score])
        print(f"Created CSV file with docking scores: {csv_filename}")
    
    return output_files

def main():
    parser = argparse.ArgumentParser(description='Split PDBQT file into individual model files and extract docking scores')
    parser.add_argument('input_file', help='Input PDBQT file containing multiple models')
    parser.add_argument('-o', '--output-dir', help='Directory to save the output files')
    parser.add_argument('--no-csv', action='store_true', help='Do not create a CSV file with docking scores')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' does not exist")
        sys.exit(1)
    
    if not args.input_file.lower().endswith('.pdbqt'):
        print(f"Warning: Input file '{args.input_file}' does not have a .pdbqt extension")
    
    split_pdbqt_models(args.input_file, args.output_dir, not args.no_csv)
    print(f"Successfully split {args.input_file} into individual model files")

if __name__ == "__main__":
    main()