"""
The module converts .pdb(3D format) files to .pdbqt files for docking
"""


import os
import subprocess
from tqdm import tqdm

class LigandConverter:
    """
    Class to convert ligand files to PDBQT format using Open Babel.

    Attributes:
        input_folder (str): Path to the folder containing ligand files.
        output_folder (str): Path to the folder where PDBQT files will be saved.
    """

    def __init__(self, input_folder, output_folder):
        """
        Initializes the LigandConverter object with input and output folders.

        Args:
            input_folder (str): Path to the folder containing ligand files.
            output_folder (str): Path to the folder where PDBQT files will be saved.
        """
        self.input_folder = input_folder
        self.output_folder = output_folder

    def convert_to_pdbqt(self, input_file, output_file):
        """
        Converts a ligand file to PDBQT format using Open Babel.

        Args:
            input_file (str): Path to the input ligand file.
            output_file (str): Path to save the output PDBQT file.
        """
        try:
            command = f'obabel -i{input_file.split(".")[-1]} {input_file} -opdbqt -O {output_file}'
            subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error converting {input_file} to PDBQT: {e}")

    def convert_all(self):
        """
        Converts all ligand files in the input folder to PDBQT format.
        """
        # Create output folder if it doesn't exist
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Get a list of ligand files in the input folder
        ligand_files = [file for file in os.listdir(self.input_folder) if file.endswith((".mol2", ".pdb", ".sdf"))]

        # Iterate over ligand files and convert each one
        for file_name in tqdm(ligand_files, desc="Converting", unit="file"):
            input_file_path = os.path.join(self.input_folder, file_name)
            output_file_path = os.path.join(self.output_folder, file_name.replace(".mol2", ".pdbqt").replace(".pdb", ".pdbqt").replace(".sdf", ".pdbqt"))
            self.convert_to_pdbqt(input_file_path, output_file_path)

if __name__ == "__main__":
    input_folder = "./3D_structures"
    output_folder = "Docking_Ligands"

    ligand_converter = LigandConverter(input_folder, output_folder)
    ligand_converter.convert_all()
