"""
Module pperforms a molecular docking pipeline
"""

import subprocess
import os
import sys

class ScriptRunner:
    """
    Class to run various scripts for molecular docking pipeline.
    """

    def run_script(self, script, args=[]):
        """
        Runs the specified script with the given arguments.

        Args:
            script (str): Name of the script to run.
            args (list): List of arguments to pass to the script.
        """
        command = ['python', script] + args
        subprocess.run(command)

    def smile_auto(self, xlsx_file):
        """
        Runs the smile_auto.py script.

        Args:
            xlsx_file (str): Path to the input XLSX file.
        """
        self.run_script('smile_auto.py', [xlsx_file])

    def pdb_to_png(self):
        """Runs the pdb_to_png.py script."""
        self.run_script('pdb_to_png.py')

    def ligand_prep(self):
        """Runs the ligand_prep.py script."""
        self.run_script('ligand_prep.py')

    def auto_dock(self, protein_folder, config_file):
        """
        Runs the auto_dock.py script.

        Args:
            protein_folder (str): Path to the folder containing protein files.
            config_file (str): Path to the configuration file.
        """
        self.run_script('auto_dock.py', [protein_folder, config_file])

    def lead_compound_identification(self):
        """Runs the lead_compound_identification.py script."""
        self.run_script('lead_compound_identification.py')

class MolecularDockingPipeline:
    """
    Class to define the steps of the molecular docking pipeline.
    """

    def __init__(self, runner):
        """
        Initializes the MolecularDockingPipeline object.

        Args:
            runner (ScriptRunner): Instance of the ScriptRunner class.
        """
        self.runner = runner

    def run_pipeline(self, protein_folder, config_file, xlsx_file):
        """
        Runs the entire molecular docking pipeline.

        Args:
            protein_folder (str): Path to the folder containing protein files.
            config_file (str): Path to the configuration file.
            xlsx_file (str): Path to the input XLSX file.
        """
        try:
            # Step 1: Generate 2D and 3D structures
            self.runner.smile_auto(xlsx_file)

            # Step 2: Convert 2D structures to PNG
            self.runner.pdb_to_png()

            # Step 3: Prepare ligands
            self.runner.ligand_prep()

            # Step 4: Perform docking
            self.runner.auto_dock(protein_folder, config_file)

            # Step 5: Identify lead compounds
            self.runner.lead_compound_identification()
        except Exception as e:
            print(f"An error occurred during pipeline execution: {e}")
            sys.exit(1)

if __name__ == "__main__":
    try:
        # Parse command line arguments
        args = sys.argv[1:]
        parsed_args = {}
        for arg in args:
            key, value = arg.split("=")
            parsed_args[key] = value

        # Check for missing arguments
        required_args = ['--protein_folder', '--config_file', '--xlsx_file']
        missing_args = [arg for arg in required_args if arg not in parsed_args]
        if missing_args:
            raise ValueError(f"Missing required argument(s): {', '.join(missing_args)}")

        # Check if files/folders exist
        for arg_key, arg_value in parsed_args.items():
            if arg_key.endswith('_folder') and not os.path.isdir(arg_value):
                raise FileNotFoundError(f"{arg_key}: {arg_value} does not exist.")
            elif arg_key.endswith('_file') and not os.path.isfile(arg_value):
                raise FileNotFoundError(f"{arg_key}: {arg_value} does not exist.")

        # Initialize the ScriptRunner and MolecularDockingPipeline objects
        script_runner = ScriptRunner()
        pipeline = MolecularDockingPipeline(script_runner)

        # Run the molecular docking pipeline
        pipeline.run_pipeline(parsed_args['--protein_folder'], parsed_args['--config_file'], parsed_args['--xlsx_file'])
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
