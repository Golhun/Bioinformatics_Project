"""
Module performs molecular docking
"""


import os
import subprocess
import csv
from tqdm import tqdm

class VinaRunner:
    """
    Class to run Autodock Vina for molecular docking.

    Attributes:
        protein_folder (str): Path to the folder containing receptor files.
        ligand_folder (str): Path to the folder containing ligand files.
        output_folder (str): Path to the folder where output files will be saved.
        conf_file (str): Path to the configuration file containing docking parameters.
    """

    def __init__(self, protein_folder, ligand_folder, output_folder, conf_file):
        """
        Initializes the VinaRunner object with input and output folders and configuration file.

        Args:
            protein_folder (str): Path to the folder containing receptor files.
            ligand_folder (str): Path to the folder containing ligand files.
            output_folder (str): Path to the folder where output files will be saved.
            conf_file (str): Path to the configuration file containing docking parameters.
        """
        self.protein_folder = protein_folder
        self.ligand_folder = ligand_folder
        self.output_folder = output_folder
        self.conf_file = conf_file

    def run_vina(self, protein_file, ligand_file, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness):
        """
        Runs Autodock Vina for molecular docking if the output files do not exist.

        Args:
            protein_file (str): Path to the receptor file.
            ligand_file (str): Path to the ligand file.
            center_x (float): X-coordinate of the center of the search space.
            center_y (float): Y-coordinate of the center of the search space.
            center_z (float): Z-coordinate of the center of the search space.
            size_x (float): Size of the search space along the X-axis.
            size_y (float): Size of the search space along the Y-axis.
            size_z (float): Size of the search space along the Z-axis.
            exhaustiveness (int): Exhaustiveness of the search algorithm.
        """
        try:
            output_base = os.path.basename(ligand_file).split('.')[0]
            output_pdbqt = os.path.join(self.output_folder, f"{output_base}_complex.pdbqt")
            output_csv = os.path.join(self.output_folder, f"{output_base}_docking_results.csv")
            log_file = os.path.join(self.output_folder, f"{output_base}_docking_log.txt")

            # Check if output files already exist
            if os.path.exists(output_pdbqt) and os.path.exists(output_csv):
                print(f"Docking results for {ligand_file} already exist. Skipping...")
                return

            # Run Vina with the calculated parameters
            subprocess.run(["vina", "--receptor", protein_file, "--ligand", ligand_file, "--out", output_pdbqt, "--log", log_file,
                            "--center_x", str(center_x), "--center_y", str(center_y), "--center_z", str(center_z),
                            "--size_x", str(size_x), "--size_y", str(size_y), "--size_z", str(size_z), "--exhaustiveness", str(exhaustiveness)],
                            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Parse the log file and extract docking results
            with open(log_file, 'r') as f:
                lines = f.readlines()
                docking_results = [line.strip().split() for line in lines if line.startswith('   ')]

            # Write docking results to a CSV file
            with open(output_csv, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Mode', 'Affinity (kcal/mol)', 'RMSD lower bound', 'RMSD upper bound'])
                writer.writerows(docking_results)
                
            # Remove the second line from the CSV file
            with open(output_csv, 'r') as f:
                lines = f.readlines()
            with open(output_csv, 'w', newline='') as f:
                f.writelines(lines[:1] + lines[2:])
        except subprocess.CalledProcessError as e:
            print(f"Error running Vina: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")
        finally:
            # Clean up the log file
            if os.path.exists(log_file):
                os.remove(log_file)

    def read_conf_file(self):
        """
        Reads parameters from the configuration file.

        Returns:
            tuple: A tuple containing receptor file, exhaustiveness, center coordinates, and search space dimensions.
        """
        receptor_file = None
        exhaustiveness = None
        center_x = None
        center_y = None
        center_z = None
        size_x = None
        size_y = None
        size_z = None

        with open(self.conf_file, 'r') as f:
            conf_data = f.readlines()

        for line in conf_data:
            if line.startswith('receptor'):
                receptor_file = line.split('=')[1].strip()
            elif line.startswith('exhaustiveness'):
                exhaustiveness = int(line.split('=')[1].strip())
            elif line.startswith('center_x'):
                center_x = float(line.split('=')[1].strip())
            elif line.startswith('center_y'):
                center_y = float(line.split('=')[1].strip())
            elif line.startswith('center_z'):
                center_z = float(line.split('=')[1].strip())
            elif line.startswith('size_x'):
                size_x = float(line.split('=')[1].strip())  # No need to convert to Angstroms
            elif line.startswith('size_y'):
                size_y = float(line.split('=')[1].strip())  # No need to convert to Angstroms
            elif line.startswith('size_z'):
                size_z = float(line.split('=')[1].strip())  # No need to convert to Angstroms

        return receptor_file, exhaustiveness, center_x, center_y, center_z, size_x, size_y, size_z

    def main(self):
        """
        Main function to run molecular docking using Vina.
        """
        try:
            if not os.path.exists(self.output_folder):
                os.makedirs(self.output_folder)

            # Read parameters from the conf.txt file
            receptor_file, exhaustiveness, center_x, center_y, center_z, size_x, size_y, size_z = self.read_conf_file()

            # Reduce the search space dimensions if necessary
            max_search_space = 27000  # Maximum search space volume
            current_volume = size_x * size_y * size_z
            if current_volume > max_search_space:
                scaling_factor = (max_search_space / current_volume) ** (1/3)
                size_x *= scaling_factor
                size_y *= scaling_factor
                size_z *= scaling_factor

            # Get list of ligand files
            ligand_files = [os.path.join(self.ligand_folder, f) for f in os.listdir(self.ligand_folder) 
                            if os.path.isfile(os.path.join(self.ligand_folder, f))]

            # Run docking for each ligand
            for ligand_file in tqdm(ligand_files, desc="Docking", unit="file"):
                protein_file = os.path.join(self.protein_folder, receptor_file)

                # Run Vina
                self.run_vina(protein_file, ligand_file, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness)
        except Exception as e:
            print(f"An error occurred: {e}")

if __name__ == "__main__":
    protein_folder = "./Docking_Proteins"
    ligand_folder = "./Docking_Ligands"
    output_folder = "results"
    conf_file = "conf.txt"

    vina_runner = VinaRunner(protein_folder, ligand_folder, output_folder, conf_file)
    vina_runner.main()
