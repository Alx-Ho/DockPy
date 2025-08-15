#!/usr/bin/env python3
"""
Convert CSV file with SMILES column to .smi format for DockPy processing.
"""

import pandas as pd
import argparse
import os

def csv_to_smi(csv_file, smiles_column, output_file, id_column=None, max_rows=None):
    """
    Convert CSV with SMILES to .smi format.
    
    Args:
        csv_file (str): Path to input CSV file
        smiles_column (str): Name of column containing SMILES strings
        output_file (str): Path to output .smi file
        id_column (str, optional): Name of column to use as identifier
        max_rows (int, optional): Maximum number of rows to process
    """
    
    # Read CSV file
    print(f"Reading CSV file: {csv_file}")
    df = pd.read_csv(csv_file)
    
    # Check if SMILES column exists
    if smiles_column not in df.columns:
        raise ValueError(f"Column '{smiles_column}' not found in CSV. Available columns: {list(df.columns)}")
    
    # Limit rows if specified
    if max_rows:
        df = df.head(max_rows)
        print(f"Processing first {len(df)} rows")
    
    # Remove rows with missing SMILES
    initial_count = len(df)
    df = df.dropna(subset=[smiles_column])
    final_count = len(df)
    
    if initial_count != final_count:
        print(f"Removed {initial_count - final_count} rows with missing SMILES")
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
    
    # Write .smi file
    print(f"Writing .smi file: {output_file}")
    with open(output_file, 'w') as f:
        for idx, row in df.iterrows():
            smiles = str(row[smiles_column]).strip()
            
            # Skip empty SMILES
            if not smiles or smiles.lower() == 'nan':
                continue
                
            if id_column and id_column in df.columns and pd.notna(row[id_column]):
                # Use provided ID column
                compound_id = str(row[id_column]).strip()
                f.write(f"{smiles} {compound_id}\n")
            else:
                # Use row index as ID
                compound_id = f"compound_{idx}"
                f.write(f"{smiles} {compound_id}\n")
    
    print(f"Successfully converted {final_count} SMILES to .smi format")

def main():
    parser = argparse.ArgumentParser(description='Convert CSV with SMILES to .smi format for DockPy')
    parser.add_argument('--csv_file', required=True, help='Path to input CSV file')
    parser.add_argument('--smiles_column', required=True, help='Name of column containing SMILES strings')
    parser.add_argument('--output_file', required=True, help='Path to output .smi file')
    parser.add_argument('--id_column', help='Name of column to use as compound identifier (optional)')
    parser.add_argument('--max_rows', type=int, help='Maximum number of rows to process (optional)')
    
    args = parser.parse_args()
    
    csv_to_smi(
        csv_file=args.csv_file,
        smiles_column=args.smiles_column,
        output_file=args.output_file,
        id_column=args.id_column,
        max_rows=args.max_rows
    )

if __name__ == "__main__":
    main()
