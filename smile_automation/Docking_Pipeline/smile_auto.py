"""
Module generates 2D and 3D structures using SMILES
"""


from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import pandas as pd
import os
import argparse

class StructureGenerator:
    """
    Class to generate 2D and 3D molecular structures from SMILES in an Excel file.

    Attributes:
        input_filename (str): Path to the Excel file containing SMILES and MMV ID columns.
    """

    def __init__(self, input_filename):
        """
        Initializes the StructureGenerator object with the provided input filename.

        Args:
            input_filename (str): Path to the Excel file containing SMILES and MMV ID columns.
        """
        self.input_filename = input_filename

    def generate_structures(self):
        """
        Generates 2D and 3D molecular structures from SMILES in the Excel file.
        """
        # Load the Excel file with SMILES
        try:
            excel_file = pd.read_excel(self.input_filename)
        except Exception as e:
            print(f'Error loading Excel file: {e}')
            return
        
        # Check if 'SMILES' column exists in the Excel file
        if 'SMILES' not in excel_file.columns:
            print("Error: 'SMILES' column not found in the Excel file.")
            return
        
        smiles_column = excel_file['SMILES']
        mmv_id_column = excel_file.get('MMV ID', None)

        # Create folders to save structures
        os.makedirs('2D_structures', exist_ok=True)
        os.makedirs('3D_structures', exist_ok=True)

        # Loop through SMILES and generate structures with progress bar
        for i, (smiles, mmv_id) in enumerate(tqdm(zip(smiles_column, mmv_id_column), desc='Processing Structures')):
            # Check for invalid or missing SMILES values
            if pd.isna(smiles) or not isinstance(smiles, str):
                print(f'Skipping invalid SMILES at index {i+1}: {smiles}')
                continue

            mol = Chem.MolFromSmiles(smiles)

            # Handle invalid SMILES
            if mol is None:
                print(f'Invalid SMILES at index {i+1}: {smiles}')
                continue

            # Add explicit hydrogens
            mol = Chem.AddHs(mol)

            # Generate 2D structure
            AllChem.Compute2DCoords(mol)
            mol_name_2d = f'o{i+1}_{mmv_id}_2D.pdb' if mmv_id else f'o{i+1}_2D.pdb'
            mol_path_2d = os.path.join('2D_structures', mol_name_2d)
            if not os.path.exists(mol_path_2d):
                Chem.MolToPDBFile(mol, mol_path_2d)

            # Generate 3D structure
            try:
                AllChem.EmbedMolecule(mol, randomSeed=42)
                mol_name_3d = f'o{i+1}_{mmv_id}_3D.pdb' if mmv_id else f'o{i+1}_3D.pdb'
                mol_path_3d = os.path.join('3D_structures', mol_name_3d)
                if not os.path.exists(mol_path_3d):
                    Chem.MolToPDBFile(mol, mol_path_3d)
            except Exception as e:
                print(f'Error generating 3D structure for SMILES {smiles}: {e}')

if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Generate 2D and 3D structures from SMILES in an Excel file.')
    parser.add_argument('filename', type=str, help='Excel file containing SMILES and MMV ID columns.')

    # Parse command line arguments
    args = parser.parse_args()

    # Check if the filename argument is provided
    if args.filename:
        structure_generator = StructureGenerator(args.filename)
        structure_generator.generate_structures()
    else:
        print('Error: Please provide the Excel filename as a command line argument.')
