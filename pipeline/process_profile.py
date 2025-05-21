import logging
import pandas as pd
import sys
import os
from pipeline.correlation_calculator import CorrelationCalculatorFactory
from pipeline.genopred_adapter import AdapterFactory
from pipeline.prevalence_calculator import PrevalenceCalculator

from sqlitedb.db_handler import DBHandler
from sqlitedb.db_loader import DBLoader
from pipeline.data_paths import DataPaths

# Logger for db_loader messages
db_logger = logging.getLogger('db_logger')
db_logger.setLevel(logging.INFO)
db_handler = logging.FileHandler('logs/db_loader.log')
db_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
db_logger.addHandler(db_handler)

def load_prs(data, gwas_code, db_handler):
    """Function to load PRS data into the database."""
    data = data.iloc[:, -2:]  # Assuming last two columns are 'iid' and 'prs_score'
    data.columns = ['iid', 'prs_score']
    data['gwas_id'] = gwas_code

    gender_table = db_handler.get_iid_gender()
    data = data.merge(gender_table, on='iid', how='left')

    data_sorted = data.sort_values(by='prs_score')

    # Compute percentiles
    for sex in ['all', 'Male', 'Female']:
        if sex in ['Male', 'Female']:
            tmp = data_sorted[data_sorted['gender'] == sex].copy()
        else:
            tmp = data_sorted.copy()

        # Normalize the prs_score (Z-score normalization)
        tmp['prs_score'] = (tmp['prs_score'] - tmp['prs_score'].mean()) / tmp['prs_score'].std()
        
        # Compute percentiles based on normalized values
        tmp['prs_percentile'] = pd.qcut(tmp['prs_score'], q=100, labels=False)
        percentile_mapping = tmp.set_index('iid')['prs_percentile']

        # Map the percentile data to data_sorted
        if sex != 'all':
            sex = 'male' if sex == 'Male' else 'female'
        data_sorted[f'prs_percentile_{sex}'] = data_sorted['iid'].map(percentile_mapping)

    data_sorted = data_sorted.rename(columns={'iid': 'indiv_id'})

    # Insert data into the 'prs' table
    db_handler.set_prs(data_sorted)

def process_file(file_path, gwas_name, db_handler, paths):
    db_logger.info(f"Processing file: {file_path}")

    adapter = AdapterFactory.create_general_adapter(file_path, db_handler, gwas_name, paths)
    prs_data = adapter.get_prs()
    load_prs(prs_data, gwas_name, db_handler)
    db_logger.info(f"Stored new prs: {gwas_name}")

    # Mapping of target types to regression types
    target_to_regression_map = {
        'Phecode': 'logistic',                    
        'ICD code': 'logistic',
    }

    # Compute correlations
    for target_type, regression_type in target_to_regression_map.items():
        # Adapt data
        specific_adapter = AdapterFactory.get_specific_adapter(adapter, target_type)
        df1, df2 = specific_adapter.get_adapted_data(target_type)
        # Get specific calculator and compute correlations
        if not df1.empty and not df2.empty:
            specific_calculator = CorrelationCalculatorFactory.get_calculator(regression_type)
            result = specific_calculator.compute_correlations(df1, df2, gwas_name, target_type, db_handler)

            # Store the result in the database
            db_handler.set_correlations(result)
            db_logger.info(f"Stored {target_type} correlations for {gwas_name}")

    # Compute prevalences
    prevalence_calculator = PrevalenceCalculator()
    prevalences = prevalence_calculator.calculate_prevalences(db_handler, gwas_name)
    db_logger.info("Calculated new prevalences")
    db_handler.set_prevalences(prevalences)
    db_logger.info(f"Stored new prevalences.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python process_profile.py <file_path> <gwas_name>")
        sys.exit(1)

    file_path = sys.argv[1]
    gwas_name = os.path.basename(file_path).split('-')[1]  # Assuming the gwas_name is derived from the filename

    # Initialize db_handler
    #Load config
    config_path = '/imppc/labs/dnalab/share/GEPETO/PolyGenie/input_data/config.ini'
    data_paths = DataPaths(config_path)

    db_handler = DBHandler(data_paths.database_file)
    
    # Define paths or pass them if needed
    paths = {}  # Replace with actual paths or configuration

    process_file(file_path, gwas_name, db_handler, paths)


