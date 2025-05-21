"""
Module: genopred_adapter
------------------------
This module defines classes for processing GenoPred output files and preparing the data for correlation calculations.

Classes:
--------
- GenopredAdapter: Base class for processing GenoPred output files and returning a DataFrame ready to be used by the CorrelationCalculator.
- MetaboliteAdapter: Adapter class for handling metabolite data in the GEPETO PRS Visualization Tool.
- DiseaseCodesAdapter: Adapter class for handling disease codes data in the GEPETO PRS Visualization Tool.
- QuestionaireLogAdapter: Adapter class for handling logistic questionnaire data data in the context of genetic prediction.
- QuestionaireLinAdapter: Adapter class for handling linear questionnaire data in the context of genetic prediction.
- GenopredAdapterFactory: Factory class for creating and returning specific GenopredAdapter instances based on the target type.

Usage:
------
The GenopredAdapter class processes the GenoPred output files, loads necessary data, and prepares it for correlation calculations.
The MetaboliteAdapter class extends GenopredAdapter to provide functionality specific to metabolite data.
The DiseaseCodesAdapter class extends GenopredAdapter to provide functionality specific to disease codes data.
The QuestionaireLogAdapter class extends GenopredAdapter to process and adapt logistic questionnaire data.
The QuestionaireLinAdapter class extends GenopredAdapter to process and adapt linear questionnaire data.
The GenopredAdapterFactory class provides methods to create and return specific GenopredAdapter instances based on the target type.
"""

import numpy as np
import pandas as pd

class GenopredAdapter:
    """
    GenopredAdapter: this class processes the GenoPred output files and returns a DataFrame ready to be used by the CorrelationCalculator.

    Attributes
    ----------
    prs : pandas.DataFrame
        DataFrame containing the processed PRS data.
    db_handler : object
        Database handler for fetching individual data.
    gwas_name : str
        Name of the GWAS study.
    dp : DataPaths
        DataPaths object containing file paths.

    Methods
    -------
    get_prs():
        Returns the processed PRS DataFrame.
    _load_pc_data():
        Loads and processes the PC data (principal components).
    _load_individual_data():
        Fetches individual data and returns a mapping (iid -> entity_id).
    _map_entity_ids(prs, pc, ehr, gen2ent):
        Maps entity IDs to PRS, PC, and EHR dataframes.
    _clean_profile(profiles_file):
        Cleans the profiles file and returns a DataFrame.
    _transform_to_z_scores(series):
        Transforms a series to Z-scores.
    """

    def __init__(self, profiles_file, db_handler, gwas_name, dp):
        """
        Initialize the GenopredAdapter with the given parameters.

        Parameters
        ----------
        profiles_file : str
            Path to the profiles file.
        db_handler : object
            Database handler for fetching individual data.
        gwas_name : str
            Name of the GWAS study.
        dp : DataPaths
            DataPaths object containing file paths.
        """

        if profiles_file is not None:
            self.prs = self._clean_profile(profiles_file)
            self.prs['prs'] = self._transform_to_z_scores(self.prs.iloc[:, 2])
            self.prs = self.prs.drop(self.prs.columns[2], axis=1)
        self.db_handler = db_handler
        self.gwas_name = gwas_name
        self.dp = dp  # Data paths container

    def get_prs(self):
        """
        Returns the processed PRS DataFrame.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the processed PRS data.
        """
        return self.prs

    def _load_pc_data(self):
        """
        Loads and processes the PC data (principal components).

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the PC data.
        """        
        return pd.read_csv(self.dp.pc_file, sep='\t')

    def _load_individual_data(self):
        """
        Fetches individual data and returns a mapping (iid -> entity_id).

        Returns
        -------
        dict
            Dictionary mapping individual IDs (iid) to entity IDs (entity_id).
        """
        consulta = self.db_handler.get_all_indiv()
        return dict(zip(consulta['iid'], consulta['entity_id']))

    def _map_entity_ids(self, prs, pc, ehr, gen2ent):
        """
        Maps entity IDs to PRS, PC, and EHR dataframes.

        Parameters
        ----------
        prs : pandas.DataFrame
            DataFrame containing the PRS data.
        pc : pandas.DataFrame
            DataFrame containing the PC data.
        ehr : pandas.DataFrame
            DataFrame containing the EHR data.
        gen2ent : dict
            Dictionary mapping individual IDs (iid) to entity IDs (entity_id).

        Returns
        -------
        tuple
            Tuple containing the mapped PRS, PC, and EHR DataFrames.
        """
        # Map entity IDs to PRS
        prs['entity_id'] = prs['IID'].map(gen2ent)
        pc['entity_id'] = pc['IID'].map(gen2ent)
        ehr['entity_id'] = ehr['indiv_id'].map(gen2ent)
        disdf = pd.DataFrame(gen2ent.values(), columns=['entity_id'])

        return prs, pc, ehr, disdf

    def _prepare_prs(self):
        """
        Formats PRS (Polygenic Risk Scores) data into quartiles and deciles.

        This method processes the PRS data associated with the instance and 
        divides it into quartiles and deciles for further analysis or visualization.

        Returns:
            dict: A dictionary containing the PRS data divided into quartiles and deciles.
        """
        return self._divide_prs(self.prs, self.gwas_name)

    def _divide_prs(self, prs, pheno):
        """
        Divides Polygenic Risk Scores (PRS) into quartiles and deciles, and assigns categorical labels.
        Parameters:
        prs (pd.DataFrame): DataFrame containing the PRS values.
        pheno (str): Phenotype name to be used as a prefix for the new columns.
        Returns:
        pd.DataFrame: DataFrame with additional columns for quartiles, deciles, and their categorical labels.
            - '{pheno}_quartile': Quartile assignment (0-3).
            - '{pheno}_quartile_low': Categorical label for quartiles ('low', 'intermediate', 'high').
            - '{pheno}_quartile_rest': Categorical label for quartiles ('rest', 'high').
            - '{pheno}_decile': Decile assignment (0-9).
            - '{pheno}_decile_low': Categorical label for deciles ('low', 'intermediate', 'high').
            - '{pheno}_decile_rest': Categorical label for deciles ('rest', 'high').
        """
        prs[f'{pheno}_quartile'] = pd.qcut(prs['prs'], q=4, labels=False)
        prs[f'{pheno}_quartile_low'] = prs[f'{pheno}_quartile'].map({0:'low',1:'intermediate',2:'intermediate', 3:'high'})
        prs[f'{pheno}_quartile_rest'] = prs[f'{pheno}_quartile'].map({0:'rest',1:'rest',2:'rest', 3:'high'})
        prs[f'{pheno}_decile'] = pd.qcut(prs['prs'], 10, labels=False)
        prs[f'{pheno}_decile_low'] = prs[f'{pheno}_decile'].map({0:'low', 9:'high'}).fillna('intermediate')
        prs[f'{pheno}_decile_rest'] = prs[f'{pheno}_decile'].map({0:'rest',1:'rest',2:'rest', 3:'rest', 4:'rest',
                                                    5:'rest', 6:'rest', 7:'rest', 8:'rest',9:'high'})
        
        
        return prs

    def _clean_profile(self, profiles_file):
        """
        Cleans the profiles file and returns a DataFrame.

        Parameters
        ----------
        profiles_file : str
            Path to the profiles file.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the cleaned profiles data.
        """

        data = pd.read_csv(profiles_file, delimiter=' ')
        first_two_cols = data[['FID', 'IID']]
        asterisk_col = data.filter(like='*')
        asterisk_col = asterisk_col.iloc[:, 0]
        data = pd.concat([first_two_cols, asterisk_col], axis=1)

        return data

    def _transform_to_z_scores(self, values):
        """
        Transforms a series to Z-scores.

        Parameters
        ----------
        series : pandas.Series
            Series to be transformed to Z-scores.

        Returns
        -------
        pandas.Series
            Series transformed to Z-scores.
        """
        mean = np.mean(values)
        std_dev = np.std(values)
        return (values - mean) / std_dev

    def _prepare_logistic_reg(self, ehr, prs, disdf, pc):
        """
        Prepares the data for logistic regression by processing electronic health records (EHR), polygenic risk scores (PRS),
        disease data frame (disdf), and principal components (pc).

        Parameters
        ----------
        ehr : pd.DataFrame
            DataFrame containing electronic health records with columns 'type' and 'target_id'.
        prs : pd.DataFrame
            DataFrame containing polygenic risk scores with columns 'entity_id' and PRS values.
        disdf : pd.DataFrame
            DataFrame containing disease information with columns 'entity_id' and disease-related data.
        pc : pd.DataFrame
            DataFrame containing principal components with columns 'entity_id', 'GENDER', and 'CURRENT_AGE'.
        
        Returns
        ----------
        tuple: A tuple containing:
            - df (pd.DataFrame): The merged DataFrame ready for logistic regression.
            - target_count (pd.DataFrame): DataFrame containing target counts filtered by count > 10.
        """
        

        target_count = []

        target_count = self.db_handler.get_target_counts(ehr['type'].iloc[0])
        target_count = target_count[target_count['phenotype_count'] > 10]
        target_count = target_count.rename(columns={"phenotype_count": "N"})

        new_cols = {}
        if not target_count.empty:
            for disease in target_count[target_count['N'] > 10]['target_id'].unique():
                new_cols[disease] = np.where(disdf['entity_id'].isin(ehr[ehr['target_id'] == disease]['indiv_id'].tolist()), 1, 0)

        new_cols_df = pd.DataFrame(new_cols, index=disdf.index)
        disdf = pd.concat([disdf, new_cols_df], axis=1)
        
        pheno = pd.merge(disdf, pc[['entity_id', 'GENDER', 'CURRENT_AGE']], on='entity_id')
        pheno = pheno.drop('GENDER', axis=1)
        pheno = pheno.drop('CURRENT_AGE', axis=1)
        pheno_prs = pd.merge(pheno, prs, on='entity_id')
        df = pd.merge(pheno_prs, pc, on='entity_id')
        
        return df, target_count

class MetaboliteAdapter(GenopredAdapter):
    """
    Adapter class for handling metabolite data in the GEPETO PRS Visualization Tool.

    This class extends the GenopredAdapter to provide functionality specific to metabolite data,
    including loading and processing PRS, PC, and EHR data.

    Attributes:
        prs (DataFrame): Polygenic Risk Scores (PRS) data.
        db_handler (object): Database handler for accessing data.
        gwas_name (str): Name of the GWAS study.
        dp (object): Data provider object containing file paths and other configurations.

    Methods:
        __init__(prs, db_handler, gwas_name, dp):
            Initializes the MetaboliteAdapter with the given PRS data, database handler, GWAS name, and data provider.

        get_adapted_data(target_type=None):
            Loads and processes the PC, EHR, and individual mapping data, formats the PRS data,
            and merges the EHR and PC data on 'entity_id'.
            Returns the merged EHR data and formatted PRS data.
    """

    def __init__(self, prs, db_handler, gwas_name, dp):
        """
        Initialize the GenopredAdapter.

        Args:
            prs (object): The PRS (Polygenic Risk Score) data.
            db_handler (object): The database handler object.
            gwas_name (str): The name of the GWAS (Genome-Wide Association Study).
            dp (object): The data processor object.
        """
        # Call the parent constructor with None for profiles_file, as it's not needed here
        super().__init__(None, db_handler, gwas_name, dp)  
        self.prs = prs  # Set the PRS data

    def get_adapted_data(self, target_type=None):
        """
        Loads and processes various datasets to produce adapted data for further analysis.

        Args:
            target_type (str, optional): The target type for adaptation. Defaults to None.

        Returns:
            tuple: A tuple containing:
                - ehr (pd.DataFrame): The merged EHR and PC data.
                - prs (pd.DataFrame): The formatted PRS data.
        """
        # Load PC data (handled by parent)
        pc = self._load_pc_data()

        # Load Metabolite EHR data
        ehr = pd.read_csv(self.dp.metab_ehr_file, sep='\t')
        ehr = ehr.rename(columns={'IID': 'indiv_id'})

        # Load individual mappings (handled by parent)
        gen2ent = self._load_individual_data()

        # Map entity IDs to PRS, PC, and EHR (handled by parent)
        prs, pc, ehr, disdf = self._map_entity_ids(self.prs, pc, ehr, gen2ent)

        # Format PRS into quartiles and deciles (handled by parent)
        prs = self._prepare_prs()

        # Merge EHR and PC data on 'entity_id'
        ehr = pd.merge(ehr, pc, on='entity_id')

        return ehr, prs
    
class DiseaseCodesAdapter(GenopredAdapter):

    def __init__(self, prs, db_handler, gwas_name, dp):
        # Call the parent constructor with None for profiles_file, as it's not needed here
        super().__init__(None, db_handler, gwas_name, dp)  
        self.prs = prs  # Set the PRS data

    def get_adapted_data(self, target_type):
        # Load PC data (handled by parent)
        pc = self._load_pc_data()

        # Load phenotype data (specific to DiseaseCodesAdapter)
        ehr = self.db_handler.get_phenotypes(target_type)

        # Load individual mappings (handled by parent)
        gen2ent = self._load_individual_data()

        # Map entity IDs to PRS, PC, and EHR (handled by parent)
        prs, pc, ehr, disdf = self._map_entity_ids(self.prs, pc, ehr, gen2ent)

        # Format PRS into quartiles and deciles (handled by parent)
        prs = self._prepare_prs()

        # Prepare logistic regression data (specific to DiseaseCodesAdapter)
        df, target_counts = self._prepare_logistic_reg(ehr, prs, disdf, pc)

        return df, target_counts

class QuestionaireLogAdapter(GenopredAdapter):
    """
    Adapter class for handling questionnaire log data in the context of genetic prediction.
    This class extends the GenopredAdapter to process and adapt data from questionnaire logs,
    polygenic risk scores (PRS), and principal components (PC) for further analysis.
    Attributes:
        prs (pd.DataFrame): DataFrame containing polygenic risk scores.
        db_handler (DatabaseHandler): Handler for database operations.
        gwas_name (str): Name of the genome-wide association study.
        dp (DataProvider): Provider for data paths and configurations.
    Methods:
        __init__(prs, db_handler, gwas_name, dp):
            Initializes the adapter with PRS data, database handler, GWAS name, and data provider.
        get_adapted_data(regression_type):
            Loads and processes the necessary data, merges it, and returns the adapted phenotypes
            and target counts for the specified regression type.
    """
    
    # Call the parent constructor with None for profiles_file, as it's not needed here
    def __init__(self, prs, db_handler, gwas_name, dp):
        super().__init__(None, db_handler, gwas_name, dp)  
        self.prs = prs  # Set the PRS data
    
    def get_adapted_data(self, regression_type):
        """
        Retrieves and processes data for regression analysis.
        This method performs the following steps:
        1. Loads principal component (PC) data.
        2. Loads phenotype data from a specified file.
        3. Filters phenotype data to keep only relevant columns.
        4. Loads individual mappings to map entity IDs.
        5. Maps entity IDs to Polygenic Risk Scores (PRS), PC, and Electronic Health Records (EHR).
        6. Formats PRS into quartiles and deciles.
        7. Merges EHR and PC data on 'entity_id'.
        8. Calculates target counts for binary targets with more than 10 occurrences.
        Args:
            regression_type (str): The type of regression to be performed.
        Returns:
            tuple: A tuple containing:
                - phenotypes (pd.DataFrame): The merged phenotype, PRS, and PC data.
                - target_count (pd.DataFrame): The counts of binary targets with more than 10 occurrences.
        """
        # Load PC data (handled by parent)
        pc = self._load_pc_data()
        
        # Load phenotypes
        phenotypes = pd.read_csv(self.dp.quest_data, sep='\t')
        targets = self.db_handler.get_target_codes('binary')
        columns_to_keep = ['entity_id'] + targets
        phenotypes = phenotypes[columns_to_keep]
        phenotypes['entity_id'] = phenotypes['entity_id'].str.replace('=', '', regex=False)

        # Load individual mappings (handled by parent)
        gen2ent = self._load_individual_data()

        # Map entity IDs to PRS, PC, and EHR (handled by parent)
        prs = self.prs
        prs['entity_id'] = prs['IID'].map(gen2ent)
        pc['entity_id'] = pc['IID'].map(gen2ent)
        
        # Format PRS into quartiles and deciles (handled by parent)
        prs = self._prepare_prs()

        # Merge EHR and PC data on 'entity_id'
        phenotypes = pd.merge(phenotypes, prs, on='entity_id')
        phenotypes = pd.merge(phenotypes, pc, on='entity_id')
 
        # Calculate target counts
        target_count = self.db_handler.get_target_counts('binary')
        #target_count = target_count[target_count['phenotype_count'] > 10]
        
        return phenotypes, target_count

class QuestionaireLinAdapter(GenopredAdapter):
    """
    Adapter class for handling questionnaire linear data in the GEPETO PRS Visualization Tool.
    This class extends the GenopredAdapter to process and adapt data from a questionnaire
    for polygenic risk score (PRS) analysis.
    Attributes:
        prs (DataFrame): DataFrame containing PRS data.
        db_handler (DatabaseHandler): Handler for database operations.
        gwas_name (str): Name of the GWAS study.
        dp (DataProvider): Data provider for accessing various data sources.
    Methods:
        __init__(prs, db_handler, gwas_name, dp):
            Initializes the adapter with PRS data, database handler, GWAS name, and data provider.
        get_adapted_data(target_type):
            Loads and adapts the data for the specified target type.
            Returns a tuple containing the adapted phenotypes and PRS data.
    """
    
    def __init__(self, prs, db_handler, gwas_name, dp):
        # Call the parent constructor with None for profiles_file, as it's not needed here
        super().__init__(None, db_handler, gwas_name, dp)  
        self.prs = prs  # Set the PRS data

    def get_adapted_data(self, target_type):
        """
        Load and adapt data for analysis.

        This method performs the following steps:
        1. Loads principal component (PC) data.
        2. Loads phenotype data from a specified file.
        3. Filters the phenotype data to keep only relevant columns.
        4. Cleans the 'entity_id' column in the phenotype data.
        5. Loads individual mappings.
        6. Maps entity IDs to Polygenic Risk Scores (PRS), PC, and Electronic Health Records (EHR).
        7. Formats PRS into quartiles and deciles.
        8. Merges EHR and PC data on 'entity_id'.

        Args:
            target_type (str): The type of target codes to retrieve (e.g., 'continuous').

        Returns:
            tuple: A tuple containing:
            - phenotypes (pd.DataFrame): The merged phenotype and PC data.
            - prs (pd.DataFrame): The formatted PRS data.
        """
        # Load PC data (handled by parent)
        pc = self._load_pc_data()

         # Load phenotypes
        phenotypes = pd.read_csv(self.dp.quest_data, sep='\t')
        targets = self.db_handler.get_target_codes('continuous')
        columns_to_keep = ['entity_id'] + targets
        phenotypes = phenotypes[columns_to_keep]
        phenotypes['entity_id'] = phenotypes['entity_id'].str.replace('=', '', regex=False)

        # Load individual mappings (handled by parent)
        gen2ent = self._load_individual_data()

        # Map entity IDs to PRS, PC, and EHR (handled by parent)
        prs = self.prs
        prs['entity_id'] = prs['IID'].map(gen2ent)
        pc['entity_id'] = pc['IID'].map(gen2ent)

        # Format PRS into quartiles and deciles (handled by parent)
        prs = self._prepare_prs()

        # Merge EHR and PC data on 'entity_id'
        phenotypes = pd.merge(phenotypes, pc, on='entity_id')

        return phenotypes, prs


class AdapterFactory:
    """
    A factory class for creating various types of GenopredAdapters.
    Methods
    -------
    create_general_adapter(profiles_file, db_handler, gwas_name, data_paths)
        Creates and returns a general GenopredAdapter.
    get_specific_adapter(general_adapter, target_type) -> GenopredAdapter
        Returns a specific adapter based on the target type.
        Supported target types:
            - 'Metabolite'
            - 'ICD code'
            - 'Phecode'
            - 'Questionaire_log'
            - 'Questionaire_lin'
        Raises:
            ValueError: If the target type is unknown.
    """
     
    @staticmethod
    def create_general_adapter(profiles_file, db_handler, gwas_name, data_paths):
        """
        Creates and returns a general GenopredAdapter.

        Args:
            profiles_file (str): The path to the profiles file.
            db_handler (object): The database handler object.
            gwas_name (str): The name of the GWAS.
            data_paths (dict): A dictionary containing data paths.

        Returns:
            GenopredAdapter: An instance of the GenopredAdapter class.
        """
        return GenopredAdapter(profiles_file, db_handler, gwas_name, data_paths)

    @staticmethod
    def get_specific_adapter(general_adapter, target_type) -> GenopredAdapter:
        """
        Returns a specific adapter instance based on the target type.

        Parameters:
        general_adapter (GenopredAdapter): The general adapter instance containing common attributes.
        target_type (str): The type of target for which a specific adapter is required. 
                           Possible values are 'Metabolite', 'ICD code', 'Phecode', 'Questionaire_log', and 'Questionaire_lin'.

        Returns:
        GenopredAdapter: An instance of a specific adapter class corresponding to the target type.

        Raises:
        ValueError: If the target type is unknown.
        """
        if target_type == 'Metabolite':
            return MetaboliteAdapter(general_adapter.prs, general_adapter.db_handler, general_adapter.gwas_name, general_adapter.dp)
        elif target_type == 'ICD code' or target_type == 'Phecode':
            return DiseaseCodesAdapter(general_adapter.prs, general_adapter.db_handler, general_adapter.gwas_name, general_adapter.dp)
        elif target_type == 'Questionaire_log':
            return QuestionaireLogAdapter(general_adapter.prs, general_adapter.db_handler, general_adapter.gwas_name, general_adapter.dp)
        elif target_type == 'Questionaire_lin':
            return QuestionaireLinAdapter(general_adapter.prs, general_adapter.db_handler, general_adapter.gwas_name, general_adapter.dp)
        else:
            raise ValueError(f"Unknown target type: {target_type}")
