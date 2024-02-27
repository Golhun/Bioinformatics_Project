"""
Module converts .pdb to .png files
"""


import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

class StructureConverter:
    """
    Class to convert PDB files to SDF and then to PNG images.

    Attributes:
        input_folder (str): Path to the folder containing PDB files.
        output_folder (str): Path to the folder where SDF and PNG files will be saved.
    """

    def __init__(self, input_folder, output_folder):
        """
        Initializes the StructureConverter object with input and output folders.

        Args:
            input_folder (str): Path to the folder containing PDB files.
            output_folder (str): Path to the folder where SDF and PNG files will be saved.
        """
        self.input_folder = input_folder
        self.output_folder = output_folder

    def convert_pdb_to_sdf(self, input_pdb, output_sdf):
        """
        Converts a PDB file to SDF format using Open Babel.

        Args:
            input_pdb (str): Path to the input PDB file.
            output_sdf (str): Path to save the output SDF file.
        """
        try:
            # Use Open Babel to convert PDB to SDF
            subprocess.run(["obabel", input_pdb, "-O", output_sdf, "-h"], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error converting {input_pdb} to SDF: {e}")

    def sdf_to_png(self, input_sdf, output_png):
        """
        Converts an SDF file to a PNG image.

        Args:
            input_sdf (str): Path to the input SDF file.
            output_png (str): Path to save the output PNG image.
        """
        try:
            # Load the molecule from the SDF file
            mol_supplier = Chem.SDMolSupplier(input_sdf)
            
            # Extract the first molecule from the supplier
            mol = next(mol_supplier, None)
            
            # Check if a molecule is obtained
            if mol is not None:
                # Generate 2D coordinates for the molecule
                AllChem.Compute2DCoords(mol)

                # Generate a PNG image from the 2D coordinates and save
                img = Draw.MolToImage(mol, size=(300, 300))
                img.save(output_png)
            else:
                print(f"Unable to process molecule from {input_sdf}")
        except Exception as e:
            print(f"Error converting {input_sdf} to PNG: {e}")

    def convert_all(self):
        """
        Converts all PDB files in the input folder to SDF and PNG files.
        """
        # Create the output folder if it doesn't exist
        os.makedirs(self.output_folder, exist_ok=True)

        # Loop through all files in the input folder
        for file_name in os.listdir(self.input_folder):
            file_path = os.path.join(self.input_folder, file_name)

            # Check if it's a file and has the .pdb extension
            if os.path.isfile(file_path) and file_name.endswith(".pdb"):
                try:
                    # Generate the corresponding SDF file path
                    sdf_file = os.path.splitext(file_name)[0] + ".sdf"
                    sdf_path = os.path.join(self.output_folder, sdf_file)

                    # Convert PDB to SDF
                    self.convert_pdb_to_sdf(file_path, sdf_path)

                    # Convert the SDF to PNG and save with the same name
                    png_file = os.path.splitext(file_name)[0] + ".png"
                    png_path = os.path.join(self.output_folder, png_file)
                    self.sdf_to_png(sdf_path, png_path)
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")

        print("Conversion completed.")

if __name__ == "__main__":
    input_folder = "2D_structures"
    output_folder = "2D_PNG2_structures"  # Change this to the desired output folder

    structure_converter = StructureConverter(input_folder, output_folder)
    structure_converter.convert_all()
