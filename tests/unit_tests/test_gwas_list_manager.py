import pytest
import pandas as pd
import os
from unittest.mock import patch

from pipeline.gwas_list_manager import GWASListManager

# 1. Test Initialization
def test_gwas_list_manager_init():
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    assert manager.metadata_path == "metadata.csv"
    assert manager.gwas_list_path == "gwas_list.txt"
    assert manager.profiles_data_folder == "profiles/"

# 2. Test Loading Metadata
@patch("pandas.read_csv")
def test_load_metadata(mock_read_csv):
    mock_read_csv.return_value = pd.DataFrame({"name": ["GWAS1"], "path": ["path1"]})
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    result = manager.load_metadata()
    assert isinstance(result, pd.DataFrame)
    mock_read_csv.assert_called_once_with("metadata.csv", na_filter=False)

# 3. Test Preparing GWAS List
def test_prepare_gwas_list():
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    # Input DataFrame
    df = pd.DataFrame({
        'name': ['GWAS1', 'GWAS2'],
        'path': ['path1', 'path2'],
        'population': ['pop1', 'pop2'],
        'n': [1000, 2000],
        'sampling': ['random', 'targeted'],
        'prevalence': [0.1, 0.2],
        'mean': [0.5, 0.6],
        'sd': [0.1, 0.2],
        'label': ['label1', 'label2']
    })
    
    # Call the method
    prepared_df = manager.prepare_gwas_list(df)
    
    # Verify that the DataFrame has the expected columns and values
    assert prepared_df.shape == df.shape
    assert (prepared_df.columns == ['name', 'path', 'population', 'n', 'sampling', 'prevalence', 'mean', 'sd', 'label']).all()
    assert prepared_df.isin(['NA']).sum().sum() == 0  # Ensure no empty cells

# 4. Test Checking Profiles Data Folder
def test_check_profiles_data_folder_exists():
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    with patch("os.path.isdir", return_value=True):
        try:
            manager.check_profiles_data_folder()
        except FileNotFoundError:
            pytest.fail("FileNotFoundError raised unexpectedly!")

def test_check_profiles_data_folder_missing():
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    with patch("os.path.isdir", return_value=False):
        with pytest.raises(FileNotFoundError):
            manager.check_profiles_data_folder()

# 5. Test Getting Profile Files
@patch("os.listdir")
def test_get_profile_files(mock_listdir):
    mock_listdir.return_value = ["file1.txt", "file2.txt"]
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    profile_files = manager.get_profile_files()
    assert profile_files == {"file1.txt", "file2.txt"}

# 6. Test Partial Match Check
def test_has_partial_match():
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    profile_files = {"file1.txt", "gwas_match.txt"}
    assert manager.has_partial_match("gwas", profile_files) is True
    assert manager.has_partial_match("GWAS2", profile_files) is False

# 7. Test Filtering GWAS List
@patch.object(GWASListManager, 'has_partial_match', side_effect=[True, False])
def test_filter_gwas_list(mock_has_partial_match):
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    df = pd.DataFrame({
        'name': ['GWAS1', 'GWAS2'],
        'path': ['path1', 'path2'],
        'population': ['pop1', 'pop2'],
        'n': [1000, 2000],
        'sampling': ['random', 'targeted'],
        'prevalence': [0.1, 0.2],
        'mean': [0.5, 0.6],
        'sd': [0.1, 0.2],
        'label': ['label1', 'label2']
    })
    
    profile_files = {"gwas1.txt"}
    filtered_df = manager.filter_gwas_list(df, profile_files)
    
    # Assert that GWAS1 was filtered out and GWAS2 remains
    assert filtered_df.shape[0] == 1
    assert filtered_df.iloc[0]['name'] == 'GWAS2'

# 8. Test Updating GWAS List
@patch("pandas.DataFrame.to_csv")
def test_update_gwas_list(mock_to_csv):
    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    df = pd.DataFrame({
        'name': ['GWAS1'],
        'path': ['path1'],
        'population': ['pop1'],
        'n': [1000],
        'sampling': ['random'],
        'prevalence': [0.1],
        'mean': [0.5],
        'sd': [0.1],
        'label': ['label1']
    })
    
    manager.update_gwas_list(df)
    mock_to_csv.assert_called_once_with("gwas_list.txt", sep=' ', header=True, index=False)

# 9. Test the Main Method: gwas_list_setup
@patch.object(GWASListManager, "load_metadata")
@patch.object(GWASListManager, "prepare_gwas_list")
@patch.object(GWASListManager, "check_profiles_data_folder")
@patch.object(GWASListManager, "get_profile_files")
@patch.object(GWASListManager, "filter_gwas_list")
@patch.object(GWASListManager, "update_gwas_list")
def test_gwas_list_setup(
    mock_update_gwas_list, mock_filter_gwas_list, mock_get_profile_files, mock_check_profiles_data_folder,
    mock_prepare_gwas_list, mock_load_metadata
):
    # Mock the DataFrame and method returns
    df = pd.DataFrame({'name': ['GWAS1'], 'path': ['path1']})
    filtered_df = pd.DataFrame({'name': ['GWAS1'], 'path': ['path1']})
    
    mock_load_metadata.return_value = df
    mock_prepare_gwas_list.return_value = filtered_df
    mock_filter_gwas_list.return_value = filtered_df
    mock_get_profile_files.return_value = {"profile1.txt"}

    manager = GWASListManager("metadata.csv", "gwas_list.txt", "profiles/")
    
    # Test the return code when GWAS list is not empty
    result = manager.gwas_list_setup()
    assert result == 0
    mock_update_gwas_list.assert_called_once_with(filtered_df)
    
    # Test the return code when filtered GWAS list is empty
    mock_filter_gwas_list.return_value = pd.DataFrame()  # Return an empty DataFrame
    result = manager.gwas_list_setup()
    assert result == 1
    mock_update_gwas_list.assert_called_once_with(filtered_df)  # Only called once
