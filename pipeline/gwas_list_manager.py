import os
import sys
import pandas as pd

class GWASListManager:
    """
    A class to manage the GWAS list and profiles data folder.
    Attributes:
        metadata_path (str): Path to the metadata CSV file.
        gwas_list_path (str): Path to the output GWAS list text file.
        profiles_data_folder (str): Path to the folder containing profile data files.
    Methods:
        load_metadata():
            Load the metadata CSV file.
        prepare_gwas_list(df):
            Prepare the GWAS list for the output text file.
        check_profiles_data_folder():
            Check if the profiles data folder exists.
        get_profile_files():
            Return the list of files in the profiles data folder.
        has_partial_match(name, profile_files):
            Check if the name partially matches any file in the profiles data folder.
        filter_gwas_list(filtered_df, profile_files):
            Filter out rows where the 'name' has a partial match with any file in the profiles data folder.
        update_gwas_list(ready_df):
            Update the gwas_list.txt file.
        gwas_list_setup():
            Main method to setup the GWAS list.
    """
    def __init__(self, metadata_path, gwas_list_path, profiles_data_folder):
        self.metadata_path = metadata_path
        self.gwas_list_path = gwas_list_path
        self.profiles_data_folder = profiles_data_folder

    def load_metadata(self):
        """Load the metadata CSV file."""
        return pd.read_csv(self.metadata_path, na_filter=False)

    def prepare_gwas_list(self, df):
        """Prepare the GWAS list for the output text file."""
        txt_columns = ['name', 'path', 'population', 'n', 'sampling', 'prevalence', 'mean', 'sd', 'label']
        filtered_df = df[txt_columns]

        # Fill empty cells with 'NA'
        return filtered_df.replace('', 'NA').fillna('NA')

    def check_profiles_data_folder(self):
        """Check if the profiles data folder exists."""
        if not os.path.isdir(self.profiles_data_folder):
            raise FileNotFoundError(f"Profiles data folder does not exist: {self.profiles_data_folder}")

    def get_profile_files(self):
        """Return the list of files in the profiles data folder."""
        return set(os.listdir(self.profiles_data_folder))

    def has_partial_match(self, name, profile_files):
        """Check if the name partially matches any file in the profiles data folder."""
        return any(name in file_name for file_name in profile_files)

    def filter_gwas_list(self, filtered_df, profile_files):
        """Filter out rows where the 'name' has a partial match with any file in the profiles data folder."""
        return filtered_df[~filtered_df['name'].apply(lambda name: self.has_partial_match(name, profile_files))]

    def update_gwas_list(self, ready_df):
        """Update the gwas_list.txt file."""
        ready_df.to_csv(self.gwas_list_path, sep=' ', header=True, index=False)

    def gwas_list_setup(self):
        """
        Main method to setup the GWAS list.
        This method performs the following steps:
        1. Loads the metadata CSV file.
        2. Prepares the GWAS list by filtering the metadata.
        3. Checks if the profiles data folder exists.
        4. Retrieves the list of profile files.
        5. Filters the GWAS list based on the profile files.
        6. Checks if there are any new GWAS files to process.
        7. Updates the gwas_list.txt file if there are new GWAS files.
        Returns:
            int: Returns 1 if there are no new GWAS files to process, otherwise returns 0.
        """
        """Main method to setup the GWAS list."""
        # Load the metadata CSV
        df = self.load_metadata()

        # Prepare the GWAS list
        filtered_df = self.prepare_gwas_list(df)

        # Check if the profiles data folder exists
        self.check_profiles_data_folder()

        # Get the list of profile files
        profile_files = self.get_profile_files()
        
        # Filter the GWAS list
        ready_df = self.filter_gwas_list(filtered_df, profile_files)
        
        # Check if there are any new GWAS files to process
        if ready_df.empty:
            print("No new GWAS in the metadata file.")
            return 1

        # Update the gwas_list.txt file
        self.update_gwas_list(ready_df)
        return 0

if __name__ == "__main__":
    # Set up the paths
    genopred_dir='/imppc/labs/dnalab/share/software/GenoPred'
    metadata_path = "input_data/gwas_metadata.csv"
    gwas_list_path = f"{genopred_dir}/pipeline/input_gcat/gwas_list.txt"
    profiles_data_folder = "profiles_data"

    # Instantiate and run the GWAS list setup
    manager = GWASListManager(metadata_path, gwas_list_path, profiles_data_folder)
    if manager.gwas_list_setup() == 1:
        sys.exit(1)
