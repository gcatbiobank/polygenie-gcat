import pytest
import pandas as pd
import pickle
from unittest.mock import MagicMock, patch

from pipeline.correlation_calculator import CorrelationCalculatorFactory, LinearCorrelationCalculator, LogisticCorrelationCalculator


@pytest.fixture
def mock_db_handler():
    # Mock the database handler with basic configurations.
    mock_db = MagicMock()
    mock_db.get_target_scope.return_value = None  # Can be adjusted as needed
    return mock_db

def load_pickled_data(file_path):
    """Helper function to load pickled data files."""
    with open(file_path, "rb") as file:
        data = pickle.load(file)
    return data[0][0], data[0][1]  # Returns df1, df2 from the nested list format

### Tests for CorrelationCalculatorFactory ###
def test_correlation_calculator_factory():
    assert isinstance(CorrelationCalculatorFactory.get_calculator('linear'), LinearCorrelationCalculator)
    assert isinstance(CorrelationCalculatorFactory.get_calculator('logistic'), LogisticCorrelationCalculator)
    with pytest.raises(ValueError):
        CorrelationCalculatorFactory.get_calculator('unsupported_type')
