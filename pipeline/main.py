import os
import sys

# Get the root directory
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the root directory to sys.path
sys.path.append(root_dir)

from sqlitedb.db_handler import DBHandler
from sqlitedb.db_loader import DBLoader
from pipeline.data_paths import DataPaths

def main():

    param_ehr = int(sys.argv[1])
    param_ind = int(sys.argv[2])
    param_pipeline = int(sys.argv[3])

    #Load config
    config_path = 'input_data/config.ini'
    data_paths = DataPaths(config_path)

    db_handler = DBHandler(data_paths.database_file)
    db_loader = DBLoader(db_handler, data_paths)

    if not os.path.exists(data_paths.database_file):
        db_handler.create_db(data_paths.schema_file, data_paths.database_file)
        db_loader.load_populations_data()
    data_paths.validate_paths()
    if param_ehr == 1:
        db_loader.load_targets_data()
    elif param_ind == 1:
        db_loader.load_individuals_data()
        db_loader.load_phenotypes()
    elif param_pipeline == 1:
        db_loader.load_pipeline_data(int(sys.argv[4]))
        #db_loader.load_pipeline_data_sge(int(sys.argv[4]))


if __name__ == '__main__':
    main()
