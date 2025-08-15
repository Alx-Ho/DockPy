# DockPy Tools

This repository hosts a suite of tools designed to facilitate molecular docking and processing using the AutoDock Vina software, specifically tailored for handling ZINC database files.

# Table of Contents
1. [DockPy Tools](#dockpy-tools)
2. [Requirements](#requirements)
3. [Conda Environment Setup](#conda-environment-setup)
    - [Creating the Environment](#creating-the-environment)
    - [Activating the Environment](#activating-the-environment)
4. [Overview for preprocess.py](#overview-for-preprocesspy)
    - [Required Arguments](#required-arguments)
    - [Optional Arguments](#optional-arguments)
    - [Example Command](#example-command)
    - [Workflow](#workflow)
    - [Notes](#notes)
5. [Overview for dock.py](#overview-for-dockpy)
    - [Required Arguments](#required-arguments-1)
    - [Optional Arguments](#optional-arguments-1)
    - [Example Command](#example-command-1)
    - [Workflow](#workflow-1)
    - [Notes](#notes-1)
6. [Preparing Custom Ligands for Docking](#preparing-custom-ligands-for-docking)
7. [Post-Processing Docking Results](#post-processing-docking-results)


## Requirements

- Conda environment with dependencies specified in `vina_env.yml`.
- AutoDock Vina and associated dependencies (`dock.py` only).

## Conda Environment Setup

Before running the script, set up the Conda environment using the provided `vina_env.yml` file. This file contains all the necessary dependencies.

### Creating the Environment

Run the following command to create the Conda environment:

```bash
conda env create -f vina_env.yml
```

### Activating the Environment

Activate the Conda environment before running the script:

```bash
conda activate vina
```

## Overview for preprocess.py

`preprocess.py` is an automated script designed for processing ZINC database files. It downloads, extracts, and converts ZINC database files, handling the workflow from downloading ZINC files, uncompressing them, formatting SMILES strings, and extracting the relevant data into a specified directory.

## Required Arguments

- `--curl_file`: Path to the .curl file containing ZINC download links.
- `--dst`: Path to the destination directory where all processed files will be saved.

## Optional Arguments
- `--smiles_limit`: The maximum number of SMILES strings to process. This limit is approximate due to the use of multiprocessing in conversions.

### Example Command

```bash
python preprocess.py --curl_file path/to/curl_file.curl --dst path/to/destination --smiles_limit 1000
```

### Workflow

1. **Directory Setup**: Creates necessary subdirectories in the destination directory.
2. **Download**: Downloads ZINC files using the specified .curl file.
3. **Uncompress**: Uncompresses the downloaded ZINC files.
4. **Format SMILES**: Converts SMILES strings to PDBQT format. Processes files recursively in the download directory.
5. **Extract**: Extracts the required data from the downloaded and processed files.

### Notes

- Ensure that all the external scripts (`download_zinc.py`, `uncompress_zinc.py`, etc.) are in the same directory as `preprocess.py` or in the system path.
- The `--smiles_limit` argument is an approximation due to the parallel processing nature of the script.

## Overview for dock.py

`dock.py` is an automated script designed for molecular docking using AutoDock Vina to automate the docking process of ligands to a given receptor and log the results. It is capable of processing multiple ligands, scoring, minimizing, and docking them to a specified receptor. The script also logs detailed information about the docking process and extracts affinity results.

### Required Arguments

- `--ligands`: Path to a ligand file or a directory containing multiple ligand files (in pdbqt format).
- `--receptor`: Path to the receptor file (in pdbqt format).
- `--center`: Center coordinates of the docking box (x y z).
- `--box_size`: Size of the docking box (x y z).

### Optional Arguments
- `--exhaustiveness`: Exhaustiveness of the search, default is 32.
- `--n_poses`: Number of docking poses to generate, default is 20.
- `--out_poses`: Number of top docking poses to output for each ligand, default is 1 (only the best pose).
- `--output_dir`: Directory to write output files, default is the current directory.
- `--overwrite`: Boolean flag to allow overwriting of output files, default is False.
- `--keep_minimized`: Boolean flag to write ligand pose after local minimization, default is False.

### Example Command

```bash
python dock.py --ligands path/to/ligands --receptor path/to/receptor.pdbqt --center 0 0 0 --box_size 20 20 20
```

### Workflow

1. **Argument Parsing**: Parses command-line arguments for ligand and receptor files, docking box specifications, and other options.
2. **Logging Setup**: Initializes logging in the specified output directory, capturing detailed information throughout the process.
3. **Ligand Processing**: Determines if the ligand path is a directory or a single file and processes each ligand file accordingly.
4. **Vina Initialization**: Sets up AutoDock Vina with the specified receptor and docking box parameters.
5. **Docking Process**: For each ligand, the script scores the pose, optionally minimizes it, docks it with the specified exhaustiveness, and writes the resulting poses.
6. **Result Extraction**: After docking, the script runs `pdbqt_extract_zincid_affinity.py` to extract and log affinity results in a CSV file.
7. **Logging and Time Tracking**: Logs detailed information for each ligand processed, including time taken and estimated time remaining for the batch.

### Notes

- Ensure the presence of AutoDock Vina in your PATH environment variable (check with `echo $PATH` in terminal).
- Remember to remove any ligands in the target box of the receptor that may interfere with docking results.
- The script provides detailed logging, useful for tracking the progress and diagnosing issues.


# Preparing Custom Ligands for Docking

This section covers how to prepare your own set of ligands (not from ZINC database) for molecular docking using the DockPy utilities.

## Overview

If you have a CSV file containing SMILES strings of compounds you want to dock, you can use the `utils/csv_to_smi.py` script to convert them to the .smi format, then process them through the standard DockPy workflow.

## Step-by-Step Process

### 1. Convert CSV to .smi Format

First, convert your CSV file containing SMILES strings to the .smi format required by DockPy:

```bash
python utils/csv_to_smi.py --csv_file your_compounds.csv --smiles_column "SMILES" --output_file compounds.smi
```

#### Required Arguments
- `--csv_file`: Path to your input CSV file
- `--smiles_column`: Name of the column containing SMILES strings
- `--output_file`: Path where the .smi file will be created

#### Optional Arguments
- `--id_column`: Name of column to use as compound identifier (if not provided, compounds will be numbered automatically)
- `--max_rows`: Maximum number of rows to process from the CSV

#### Example Commands

**Basic conversion:**
```bash
python utils/csv_to_smi.py --csv_file molecules.csv --smiles_column "smiles" --output_file ligands.smi
```

**With custom ID column:**
```bash
python utils/csv_to_smi.py --csv_file molecules.csv --smiles_column "canonical_smiles" --id_column "compound_id" --output_file ligands.smi
```

**Process only first 500 compounds:**
```bash
python utils/csv_to_smi.py --csv_file molecules.csv --smiles_column "smiles" --output_file ligands.smi --max_rows 500
```

### 2. Convert SMILES to PDBQT Format

Once you have the .smi file, convert the SMILES strings to PDBQT format for docking:

```bash
python utils/smiles_to_pdbqt.py --smiles_file ligands.smi --dst ligand_pdbqts/ --num_processes 8
```

#### Arguments
- `--smiles_file`: Path to the .smi file created in step 1
- `--dst`: Output directory for PDBQT files
- `--num_processes`: Number of CPU cores to use (optional, default=4)
- `--smiles_limit`: Maximum number of SMILES to process (optional, -1 for all)

### 3. Run Molecular Docking

Now you can use the generated PDBQT files with the main docking script:

```bash
python dock.py --ligands ligand_pdbqts/ --receptor receptor.pdbqt --center 0 0 0 --box_size 20 20 20 --output_dir docking_results/
```

## Complete Workflow Example

Here's a complete example starting from a CSV file:

```
# 1. Convert CSV to .smi format
python utils/csv_to_smi.py \
    --csv_file my_compounds.csv \
    --smiles_column "SMILES" \
    --id_column "compound_name" \
    --output_file my_ligands.smi \
    --max_rows 1000

# 2. Convert SMILES to PDBQT
python utils/smiles_to_pdbqt.py \
    --smiles_file my_ligands.smi \
    --dst prepared_ligands/ \
    --num_processes 8

# 3. Run docking
python dock.py \
    --ligands prepared_ligands/ \
    --receptor protein.pdbqt \
    --center 10.5 -5.2 8.1 \
    --box_size 25 25 25 \
    --output_dir docking_results/ \
    --exhaustiveness 32
```

## ChEMBL Small Molecule Drugs Example

For the included ChEMBL small molecule drugs dataset:

```
# 1. Convert ChEMBL CSV to .smi format
python utils/csv_to_smi.py \
    --csv_file chembl_sm_drugs_08152025.csv \
    --smiles_column "Smiles" \
    --id_column "Parent Molecule" \
    --output_file chembl_sm_drugs.smi

# 2. Convert SMILES to PDBQT (process first 1000 for testing)
python utils/smiles_to_pdbqt.py \
    --smiles_file chembl_sm_drugs.smi \
    --dst chembl_ligands/ \
    --num_processes 8 \
    --smiles_limit 1000

# 3. Run docking against your target
python dock.py \
    --ligands chembl_ligands/ \
    --receptor your_receptor.pdbqt \
    --center X Y Z \
    --box_size 25 25 25 \
    --output_dir chembl_docking_results/
```

## CSV File Requirements

Your CSV file should:
- Have a header row with column names
- Contain a column with valid SMILES strings
- Optionally include a column with unique compound identifiers
- Use standard CSV formatting (comma-separated values)

**Example CSV structure:**
```csv
compound_id,SMILES,molecular_weight
COMP001,CCO,46.07
COMP002,CC(=O)O,60.05
COMP003,c1ccccc1,78.11
```

## Notes

- The script will automatically skip rows with missing or invalid SMILES strings
- If no ID column is specified, compounds will be automatically numbered as `compound_0`, `compound_1`, etc.
- Large datasets can be processed in chunks using the `--max_rows` parameter
- The conversion process may take considerable time for large compound libraries due to 3D structure generation
- For the ChEMBL dataset, consider starting with a subset using `--max_rows` to test your workflow before processing the entire dataset

# Post-Processing Docking Results

After running molecular docking with `dock.py`, you may want to analyze the results and convert the output files to different formats. DockPy provides utilities for selecting the best-performing ligands and converting PDBQT files to other molecular formats.

## Selecting Top-Performing Ligands

### Overview for get_top_n_affinities.py

`get_top_n_affinities.py` is a utility script that identifies and copies the top N molecules with the best (lowest) binding affinities from your docking results. This is useful for focusing on the most promising compounds for further analysis.

#### Required Arguments

- `input_dir`: Directory containing the docking results with `affinity_results.csv`
- `output_dir`: Directory where the top molecules will be copied
- `n`: Number of top molecules to select (integer)

#### Example Commands

**Select top 10 molecules:**
```bash
python utils/get_top_n_affinities.py docking_results/ top_10_molecules/ 10
```

**Select top 50 molecules from a large screening:**
```bash
python utils/get_top_n_affinities.py large_screen_results/ best_50_hits/ 50
```

#### Workflow

1. **CSV Analysis**: Reads the `affinity_results.csv` file from the input directory
2. **Sorting**: Sorts molecules by binding affinity (lowest/best values first)
3. **Selection**: Selects the top N molecules based on affinity scores
4. **File Copying**: Copies the corresponding PDBQT files to the output directory with ranked filenames
5. **Ranking**: Renames files with rank prefixes (e.g., `rank_01_`, `rank_02_`, etc.)

#### Output Structure

The script creates an output directory containing:
- `affinity_results.csv`: Copy of the original results file
- `rank_01_[original_filename].pdbqt`: Best-scoring molecule
- `rank_02_[original_filename].pdbqt`: Second-best molecule
- And so on...

#### Notes

- The script automatically determines the number of digits needed for ranking based on N
- Files are ranked from best (lowest affinity) to worst within the selected top N
- The original `affinity_results.csv` is copied to the output directory for reference
- If a PDBQT file cannot be found, a warning is displayed and the file is skipped

## Converting PDBQT Files to Other Formats

### Overview for pdbqt_converter.py

`pdbqt_converter.py` converts PDBQT files (the output format from AutoDock Vina) to more commonly used molecular formats like PDB or SDF. This is useful for visualization in molecular viewers or further analysis with other software tools.

#### Required Arguments

- `--pdbqt_dir`: Directory containing PDBQT files to convert

#### Optional Arguments

- `--dst`: Output directory for converted files (default: current directory)
- `--num_processes`: Number of CPU cores to use for parallel processing (default: 4, use -1 for all cores)
- `--file_limit`: Maximum number of files to process (default: -1 for all files)
- `--output_format`: Output format - either 'pdb' or 'sdf' (default: 'pdb')

#### Example Commands

**Convert all PDBQT files to PDB format:**
```bash
python utils/pdbqt_converter.py --pdbqt_dir docking_results/ --dst converted_pdb/ --output_format pdb
```

**Convert top molecules to SDF format:**
```bash
python utils/pdbqt_converter.py --pdbqt_dir top_10_molecules/ --dst top_molecules_sdf/ --output_format sdf --num_processes 8
```

**Convert only first 100 files for testing:**
```bash
python utils/pdbqt_converter.py --pdbqt_dir large_results/ --dst test_conversion/ --file_limit 100
```

#### Format Differences

**PDB Format:**
- Contains only the first conformer of the first molecule
- Widely supported by molecular visualization software
- Good for structural analysis and visualization

**SDF Format:**
- Contains all conformers of all molecules
- Preserves multiple binding poses from docking
- Better for comprehensive analysis of docking results
- Suitable for cheminformatics workflows

#### Workflow

1. **File Discovery**: Identifies all PDBQT files in the specified directory
2. **Parallel Processing**: Uses multiprocessing to convert files efficiently
3. **Format Conversion**: Converts PDBQT to the specified output format using RDKit and Meeko
4. **Output Generation**: Creates converted files with the same base name but different extension

#### Notes

- The script preserves the original filename but changes the extension (e.g., `molecule.pdbqt` → `molecule.pdb`)
- Uses multiprocessing for faster conversion of large datasets
- Automatically skips files that already exist in the output directory
- Provides detailed logging of the conversion process
- Requires RDKit and Meeko libraries (included in the conda environment)

## Complete Post-Processing Workflow Example

Here's a typical workflow for analyzing docking results:

```bash
# 1. Run docking (assuming this was already completed)
python dock.py --ligands ligands/ --receptor protein.pdbqt --center 0 0 0 --box_size 25 25 25 --output_dir docking_results/

# 2. Select top 20 molecules with best affinities
python utils/get_top_n_affinities.py docking_results/ top_20_hits/ 20

# 3. Convert the top hits to PDB format for visualization
python utils/pdbqt_converter.py --pdbqt_dir top_20_hits/ --dst top_20_pdb/ --output_format pdb

# 4. Optionally, convert to SDF format for cheminformatics analysis
python utils/pdbqt_converter.py --pdbqt_dir top_20_hits/ --dst top_20_sdf/ --output_format sdf
```

This workflow allows you to efficiently identify promising compounds from large virtual screening campaigns and prepare them for further analysis or experimental validation.
