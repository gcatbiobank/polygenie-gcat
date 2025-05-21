import pytest
import pandas as pd
from unittest.mock import MagicMock
from pipeline.prevalence_calculator import PrevalenceCalculator

@pytest.fixture
def mock_db_handler():
    """Fixture to mock the database handler."""
    db_handler = MagicMock()
    db_handler.get_gwas_codes.return_value = ['GWAS1', 'GWAS2']
    db_handler.get_target_codes.side_effect = lambda x: ['T1', 'T2'] if x in ['Phecode', 'ICD code'] else ['T3']
    db_handler.get_all_risk_scores.return_value = pd.DataFrame({
        'gwas_id': ['GWAS1', 'GWAS1', 'GWAS2'],
        'indiv_id': [1, 2, 3],
        'prs_percentile_all': [90, 80, 70],
        'prs_percentile_male': [90, None, 70],
        'prs_percentile_female': [None, 80, None]
    })
    db_handler.get_phenotypes.side_effect = lambda x: pd.DataFrame({
        'target_id': ['T1', 'T1', 'T2'],
        'indiv_id': [1, 4, 2]
    })
    return db_handler

@pytest.fixture
def prevalence_calculator():
    """Fixture for the PrevalenceCalculator instance."""
    return PrevalenceCalculator()

def test_process_batch(prevalence_calculator, mock_db_handler):
    """Test the _process_batch method."""
    batch = [('GWAS1', 'T1'), ('GWAS2', 'T2')]
    results = prevalence_calculator._process_batch(batch, mock_db_handler)

    # Verify the result structure
    assert isinstance(results, list)
    assert all(isinstance(item, tuple) for item in results)

