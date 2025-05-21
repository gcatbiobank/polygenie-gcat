import pytest
from unittest import mock

from pipeline.data_paths import DataPaths


# Sample content of the config file
CONFIG_CONTENT = """
[Paths]
database_file = sqlitedb/polygenie.db
schema_file = sqlitedb/schema.sql
metadata_file = input_data/gwas_metadata.csv
data_folder = app/data/
metab_file = input_data/metabolites.tsv
icd_file = input_data/icdcodes.tsv
phecodes_file = input_data/phecodes.tsv
indiv_file = input_data/individuals.tsv
phenotypes_file = input_data/phenotypes.tsv
cohorts_file = input_data/cohort.tsv
population_file = input_data/populations.tsv
"""

@pytest.fixture
def config_file(tmp_path):
    """Creates a temporary config file for testing."""
    config_path = tmp_path / "test_config.ini"
    config_path.write_text(CONFIG_CONTENT)
    return str(config_path)


def test_load_from_config_file_not_found():
    # Mock os.path.exists to return False (simulate file not found)
    with mock.patch('os.path.exists', return_value=False):
        with pytest.raises(FileNotFoundError, match="Config file 'nonexistent_config.ini' does not exist."):
            DataPaths('nonexistent_config.ini')


def test_load_from_config_success(config_file):
    # Mock os.path.exists to return True (simulate file exists)
    with mock.patch('os.path.exists', return_value=True):
        # Create a DataPaths instance using the config file
        dp = DataPaths(config_file)

        # Assert that the attributes have the correct values from the config file
        assert dp.database_file == 'sqlitedb/polygenie.db'
        assert dp.schema_file == 'sqlitedb/schema.sql'
        assert dp.metadata_file == 'input_data/gwas_metadata.csv'
        assert dp.metab_file == 'input_data/metabolites.tsv'
        assert dp.icd_file == 'input_data/icdcodes.tsv'
        assert dp.phecodes_file == 'input_data/phecodes.tsv'
        assert dp.indiv_file == 'input_data/individuals.tsv'
        assert dp.phenotypes_file == 'input_data/phenotypes.tsv'
        assert dp.cohorts_file == 'input_data/cohort.tsv'
        assert dp.population_file == 'input_data/populations.tsv'


def test_fallback_values(config_file, tmp_path):
    # Modify the config content to remove some keys
    incomplete_config = """
    [Paths]
    database_file = sqlitedb/polygenie.db
    schema_file = sqlitedb/schema.sql
    """
    incomplete_config_path = tmp_path / "incomplete_config.ini"
    incomplete_config_path.write_text(incomplete_config)

    # Mock os.path.exists to return True
    with mock.patch('os.path.exists', return_value=True):
        # Create a DataPaths instance using the incomplete config file
        dp = DataPaths(str(incomplete_config_path))

        # Assert that the missing values fall back to None
        assert dp.database_file == 'sqlitedb/polygenie.db'
        assert dp.schema_file == 'sqlitedb/schema.sql'
        assert dp.metadata_file is None  # This was missing in the config
        assert dp.metab_file is None     # This was missing in the config
        assert dp.icd_file is None       # This was missing in the config
        assert dp.phecodes_file is None  # This was missing in the config
        assert dp.indiv_file is None     # This was missing in the config
        assert dp.phenotypes_file is None  # This was missing in the config
        assert dp.cohorts_file is None     # This was missing in the config
        assert dp.population_file is None  # This was missing in the config


def test_validate_paths_all_valid(config_file):
    # Simulate all files and directory are valid
    with mock.patch('os.path.isfile', return_value=True), \
         mock.patch('os.path.isdir', return_value=True):
        
        data_paths = DataPaths(config_file)
        
        # Should not raise any errors
        data_paths.validate_paths()
