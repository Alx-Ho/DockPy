import argparse
from meeko import MoleculePreparation, PDBQTWriterLegacy
import rdkit
import os

# Create an argument parser
parser = argparse.ArgumentParser(description='Convert PDB to PDBQT for molecular docking.')
parser.add_argument('pdb_file', type=str, help='Input PDB file')
parser.add_argument('-o', '--output', type=str, help='Optional output PDBQT file location', default=None)

# Parse arguments
args = parser.parse_args()

ligand_pdb_file = args.pdb_file
output_file = args.output

# Default output file name if not specified
if not output_file:
    output_file = os.path.splitext(ligand_pdb_file)[0] + '.pdbqt'
    print(f"Output file not specified. Using default: {output_file}")

print("Reading the molecule from the PDB file...")
# Read the molecule
lig = rdkit.Chem.MolFromPDBFile(ligand_pdb_file, removeHs=False)

print("Protonating the molecule...")
protonated_lig = rdkit.Chem.AddHs(lig)

print("Embedding the molecule...")
rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)

print("Preparing the molecule...")
# Prepare the molecule
preparator = MoleculePreparation()
mol_setups = preparator.prepare(protonated_lig)

print("Writing the molecule to PDBQT format...")
# Write to PDBQT
for setup in mol_setups:
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
    if is_ok:
        with open(output_file, 'w') as f:
            f.write(pdbqt_string)
        print(f"Successfully wrote to {output_file}")
    else:
        print(f"Error in writing to PDBQT: {error_msg}")

print("Conversion complete.")
