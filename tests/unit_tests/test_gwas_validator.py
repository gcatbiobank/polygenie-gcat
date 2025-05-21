import pytest
import pandas as pd
from unittest.mock import patch, mock_open

from pipeline.gwas_validator import GWASValidator

# Test data
test_columns = ["SNP", "A1", "A2", "BETA", "P"]
test_df = pd.DataFrame(columns=test_columns)


# 1. Test Initialization and Reverse Dictionary Creation
def test_gwas_validator_init():
    validator = GWASValidator()
    assert validator.log_file == "pipeline/logs/header_check.log"
    assert isinstance(validator.reverse_dict, dict)
    assert validator.reverse_dict["SNP"] == "SNP"
    assert validator.reverse_dict["ALLELE1"] == "A1"


# 2. Test Mapping Columns
def test_map_columns():
    validator = GWASValidator()
    input_columns = ["SNP", "ALLELE1", "P_VALUE"]
    mapped_columns = validator.map_columns(input_columns)
    
    expected_output = ["SNP", "A1", "P"]
    assert mapped_columns == expected_output


# 3. Test Checking Required Columns
def test_check_required_columns_pass():
    validator = GWASValidator()
    test_data = pd.DataFrame(columns=["SNP", "A1", "A2", "BETA", "P"])
    
    with patch("builtins.open", mock_open()) as mocked_file:
        result = validator.check_required_columns(test_data, "test_file.csv")
        assert result is True
        mocked_file.assert_called_once_with("pipeline/logs/header_check.log", 'a')


def test_check_required_columns_fail():
    validator = GWASValidator()
    test_data = pd.DataFrame(columns=["SNP", "A1", "A2", "BETA"])  # Missing P column
    
    with patch("builtins.open", mock_open()) as mocked_file:
        result = validator.check_required_columns(test_data, "test_file.csv")
        assert result is False
        mocked_file.assert_called_once_with("pipeline/logs/header_check.log", 'a')


# 4. Test Reading Files
@patch("pandas.read_csv")
def test_read_file(mock_read_csv):
    validator = GWASValidator()
    mock_read_csv.return_value = pd.DataFrame(columns=["SNP", "A1", "A2", "BETA", "P"])

    result = validator.read_file("test_file.csv")
    assert isinstance(result, pd.DataFrame)
    mock_read_csv.assert_called_once_with("test_file.csv", delimiter=',', nrows=1)


# 5. Test Validate Files
@patch.object(GWASValidator, "read_file")
@patch.object(GWASValidator, "check_required_columns")
def test_validate_files(mock_check_required_columns, mock_read_file):
    mock_check_required_columns.return_value = False
    mock_read_file.return_value = test_df

    validator = GWASValidator()
    test_files_df = pd.DataFrame({"path": ["file1.csv", "file2.csv"]})

    with patch("os.path.exists", return_value=True), patch("os.remove"):
        with pytest.raises(SystemExit):  # Expecting SystemExit due to failed validation
            validator.validate_files(test_files_df)

    mock_read_file.assert_called()
    mock_check_required_columns.assert_called()
