import configparser
import os

class DataPaths:
    """
    A class to manage file paths for various data files used in the pipeline.

    Attributes
    ----------
    paths : dict
        A dictionary to store file paths.

    Methods
    -------
    load_from_config(config_path: str):
        Load file paths from the given configuration file.
    """

    def __init__(self, config_file):
        """
        Initialize DataPaths with default values and load from the config file.

        Parameters
        ----------
        config_file : str
            Path to the configuration file.
        """
        self.database_file = None
        self.schema_file = None
        self.metadata_file = None
        self.metab_file = None
        self.icd_file = None
        self.phecodes_file = None
        self.indiv_file = None
        self.phenotypes_file = None
        self.cohorts_file = None
        self.population_file = None
        self.pc_file = None
        self.metab_ehr_file = None
        self.quest_data = None
        self.quest_meta = None
        
        # Load values from the config file
        self.load_from_config(config_file)

    def load_from_config(self, config_path):
        """
        Load file paths from the given configuration file.

        Parameters
        ----------
        config_path : str
            Path to the configuration file.

        Raises
        ------
        FileNotFoundError
            If the configuration file does not exist.
        """
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file '{config_path}' does not exist.")
        
        config = configparser.ConfigParser()
        config.read(config_path)
        
        # Update each attribute from the config file
        self.database_file = config.get('Paths', 'database_file', fallback=self.database_file)
        self.schema_file = config.get('Paths', 'schema_file', fallback=self.schema_file)
        self.metadata_file = config.get('Paths', 'metadata_file', fallback=self.metadata_file)
        self.metab_file = config.get('Paths', 'metab_file', fallback=self.metab_file)
        self.icd_file = config.get('Paths', 'icd_file', fallback=self.icd_file)
        self.phecodes_file = config.get('Paths', 'phecodes_file', fallback=self.phecodes_file)
        self.indiv_file = config.get('Paths', 'indiv_file', fallback=self.indiv_file)
        self.phenotypes_file = config.get('Paths', 'phenotypes_file', fallback=self.phenotypes_file)
        self.cohorts_file = config.get('Paths', 'cohorts_file', fallback=self.cohorts_file)
        self.population_file = config.get('Paths', 'population_file', fallback=self.population_file)
        self.pc_file = config.get('Paths', 'pc_file', fallback=self.pc_file)
        self.metab_ehr_file = config.get('Paths', 'metab_ehr_file', fallback=self.metab_ehr_file)
        self.quest_data = config.get('Paths','questionaire_data_file', fallback=self.quest_data)
        self.quest_meta = config.get('Paths', 'questionaire_metadata_file', fallback=self.quest_meta)
    
    def validate_paths(self):
        invalid_paths = []

        # Check if file paths exist
        if not os.path.isfile(self.database_file):
            invalid_paths.append(f"Invalid file path: {self.database_file}")
        if not os.path.isfile(self.schema_file):
            invalid_paths.append(f"Invalid file path: {self.schema_file}")
        if not os.path.isfile(self.metadata_file):
            invalid_paths.append(f"Invalid file path: {self.metadata_file}")
        if not os.path.isfile(self.metab_file):
            invalid_paths.append(f"Invalid file path: {self.metab_file}")
        if not os.path.isfile(self.icd_file):
            invalid_paths.append(f"Invalid file path: {self.icd_file}")
        if not os.path.isfile(self.phecodes_file):
            invalid_paths.append(f"Invalid file path: {self.phecodes_file}")
        if not os.path.isfile(self.indiv_file):
            invalid_paths.append(f"Invalid file path: {self.indiv_file}")
        if not os.path.isfile(self.phenotypes_file):
            invalid_paths.append(f"Invalid file path: {self.phenotypes_file}")
        if not os.path.isfile(self.cohorts_file):
            invalid_paths.append(f"Invalid file path: {self.cohorts_file}")
        if not os.path.isfile(self.population_file):
            invalid_paths.append(f"Invalid file path: {self.population_file}")
        if not os.path.isfile(self.pc_file):
            invalid_paths.append(f"Invalid file path: {self.pc_file}")
        if not os.path.isfile(self.metab_ehr_file):
            invalid_paths.append(f"Invalid file path: {self.metab_ehr_file}")
        if not os.path.isfile(self.quest_meta):
            invalid_paths.append(f"Invalid file path: {self.quest_meta}")
        if not os.path.isfile(self.quest_data):
            invalid_paths.append(f"Invalid file path: {self.quest_data}")
        
        # Raise an error if any invalid paths are found
        if invalid_paths:
            raise FileNotFoundError(f"The following paths are invalid:\n" + "\n".join(invalid_paths))
        
