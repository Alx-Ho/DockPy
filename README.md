# AutoDock ZINC Tools

This repository hosts a suite of tools designed to facilitate molecular docking and processing using the AutoDock Vina software, specifically tailored for handling ZINC database files.

# Table of Contents
1. [AutoDock ZINC Tools](#autodock-zinc-tools)
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
- `--out_poses`: Number of docking poses to output, default is 1.
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
- The script provides detailed logging, useful for tracking the progress and diagnosing issues.
