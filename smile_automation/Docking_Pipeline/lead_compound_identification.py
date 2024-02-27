"""
Module identifies lead compounds based on highest affinity
"""


import os
import glob
import pandas as pd

class DockingResultProcessor:
    """
    Class to process docking results and identify top molecules.
    """

    def __init__(self, folder_path, top_molecules):
        """
        Initializes the DockingResultProcessor object.

        Args:
            folder_path (str): Path to the folder containing docking result CSV files.
            top_molecules (int): Number of top molecules to select.
        """
        self.folder_path = folder_path
        self.top_molecules = top_molecules

    def process_docking_results(self):
        """
        Processes docking results to identify top compounds.

        Returns:
            list: List of dictionaries containing file names and their highest affinity scores.
        """
        # Initialize an empty list to store file names and their corresponding highest affinity scores
        file_scores = []

        # Iterate through CSV files in the folder
        for file_path in glob.glob(os.path.join(self.folder_path, "*.csv")):
            print(f"Processing file: {file_path}")
            try:
                # Load CSV file
                df = pd.read_csv(file_path)
                print(f"Columns in DataFrame: {df.columns}")

                # Check if 'Affinity (kcal/mol)' column exists
                if 'Affinity (kcal/mol)' not in df.columns:
                    print(f"Warning: 'Affinity (kcal/mol)' column not found in {file_path}")
                    continue

                # Extract numeric values from the 'Affinity (kcal/mol)' column
                affinities = df['Affinity (kcal/mol)'].apply(pd.to_numeric, errors='coerce').dropna().values

                if len(affinities) == 0:
                    print(f"Warning: No valid numeric values found in 'Affinity (kcal/mol)' column of {file_path}")
                    continue

                # Get the highest affinity score
                max_affinity = max(affinities)

                # Store the file name and its highest affinity score
                file_name = os.path.basename(file_path)
                file_scores.append({'File_Name': file_name, 'Highest_Affinity': max_affinity})
            except Exception as e:
                print(f"Error processing file {file_path}: {e}")

        # Sort files by their highest affinity scores in ascending order
        sorted_files = sorted(file_scores, key=lambda x: x['Highest_Affinity'])

        # Select top molecules based on highest affinity scores
        top_files = sorted_files[:self.top_molecules]

        return top_files

if __name__ == "__main__":
    try:
        # Define the folder containing docking results
        docking_folder = "./results"

        # Define the number of top molecules to select
        top_molecules = 20

        # Initialize the DockingResultProcessor object
        processor = DockingResultProcessor(docking_folder, top_molecules)

        # Process docking results to identify top compounds
        selected_files = processor.process_docking_results()
        print("Selected Files:", selected_files)

        # Create a folder for lead compounds if it doesn't exist
        lead_folder = "Lead_compounds"
        if not os.path.exists(lead_folder):
            os.makedirs(lead_folder)

        # Save selected files to a CSV file in the Lead_compounds folder
        csv_file_path = os.path.join(lead_folder, "lead_compound_files.csv")
        df = pd.DataFrame(selected_files)
        df.to_csv(csv_file_path, index=False)
    except Exception as e:
        print(f"An error occurred: {e}")
