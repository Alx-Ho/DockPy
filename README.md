# AutoDock ZINC Tools

## Requirements

- Conda environment with dependencies specified in `vina_env.yml`.

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

Execute `preprocess.py` from the command line with the following arguments:

- `--curl_file`: Path to the .curl file containing ZINC download links.
- `--dst`: Path to the destination directory where all processed files will be saved.
- `--smiles_limit`: (Optional) The maximum number of SMILES strings to process. This limit is approximate due to the use of multiprocessing in conversions.

### Example Command

```bash
python preprocess.py --curl_file path/to/curl_file.curl --dst path/to/destination --smiles_limit 1000
```

## Workflow

1. **Directory Setup**: Creates necessary subdirectories in the destination directory.
2. **Download**: Downloads ZINC files using the specified .curl file.
3. **Uncompress**: Uncompresses the downloaded ZINC files.
4. **Format SMILES**: Converts SMILES strings to PDBQT format. Processes files recursively in the download directory.
5. **Extract**: Extracts the required data from the downloaded and processed files.

## Notes

- Ensure that all the external scripts (`download_zinc.py`, `uncompress_zinc.py`, etc.) are in the same directory as `preprocess.py` or in the system path.
- The `--smiles_limit` argument is an approximation due to the parallel processing nature of the script.
