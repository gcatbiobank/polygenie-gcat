import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock

from pipeline.genopred_adapter import AdapterFactory, DiseaseCodesAdapter, GenopredAdapter, MetaboliteAdapter, QuestionaireLinAdapter, QuestionaireLogAdapter

# Sample test data
sample_profile_data = pd.DataFrame({
    'FID': ['F1', 'F2'],
    'IID': ['I1', 'I2'],
    'ad_ldak_Model1': [1.5, 2.5],
    'ad_ldak_Model131*': [0.5, 1.5]
})

sample_cleaned_profile = pd.DataFrame({
    'FID': ['F1', 'F2'],
    'IID': ['I1', 'I2'],
    'ad_ldak_Model131*': [0.5, 1.5]
})

class TestGenopredAdapter:
    
    def test_clean_profile(self):
        adapter = GenopredAdapter.__new__(GenopredAdapter)  # Bypass __init__

        # Mocking read_csv inside _clean_profile to simulate the file reading
        with patch('pandas.read_csv', return_value=sample_profile_data):
            result = adapter._clean_profile('mock_profiles_file')

        # Expected result
        expected_result = sample_cleaned_profile

        # Assert the cleaned DataFrame matches expected output
        pd.testing.assert_frame_equal(result, expected_result)

    def test_transform_to_z_scores(self):
        adapter = GenopredAdapter.__new__(GenopredAdapter)  # Bypass __init__

        # Sample values to transform
        values = pd.Series([1.0, 2.0, 3.0])
        
        # Expected Z-score transformation
        expected_z_scores = (values - np.mean(values)) / np.std(values)

        # Call the method and assert the result
        result = adapter._transform_to_z_scores(values)
        pd.testing.assert_series_equal(result, expected_z_scores, check_names=False)


    @patch('pipeline.genopred_adapter.GenopredAdapter._clean_profile')
    @patch('pipeline.genopred_adapter.GenopredAdapter._transform_to_z_scores')
    def test_init(self, mock_transform_to_z_scores, mock_read_csv):
        # Setup mock return values
        mock_read_csv.return_value = sample_cleaned_profile  # Mock the cleaned profile data from pandas.read_csv
        mock_transform_to_z_scores.return_value = pd.Series([0, 1])

        # Initialize db_handler and gwas_name mock objects
        mock_db_handler = MagicMock()
        mock_dp = MagicMock()
        gwas_name = 'test_gwas'

        # Initialize the GenopredAdapter class
        adapter = GenopredAdapter(profiles_file='mock_profiles_file', db_handler=mock_db_handler, dp=mock_dp, gwas_name=gwas_name)

        # Check if prs['prs'] column is created correctly
        pd.testing.assert_series_equal(adapter.prs['prs'], pd.Series([0, 1]), check_names=False)

        # Check if the third column ('ad_ldak_Model131*') is dropped from prs
        assert 'ad_ldak_Model131*' not in adapter.prs.columns

    @patch('pipeline.genopred_adapter.GenopredAdapter._clean_profile')
    def test_map_entity_ids(self, mock_read_csv):
        # Sample data
        prs = pd.DataFrame({'IID': ['I1', 'I2']})
        pc = pd.DataFrame({'IID': ['I1', 'I2']})
        ehr = pd.DataFrame({'indiv_id': ['I1', 'I2']})
        gen2ent = {'I1': 101, 'I2': 102}

        # Initialize the GenopredAdapter class
        adapter = GenopredAdapter(profiles_file='mock_profiles_file', db_handler=MagicMock(), dp=MagicMock(), gwas_name='test_gwas')

        # Call _map_entity_ids and verify result
        mapped_prs, mapped_pc, mapped_ehr, mapped_disdf = adapter._map_entity_ids(prs, pc, ehr, gen2ent)

        # Expected data
        expected_prs = prs.copy()
        expected_prs['entity_id'] = [101, 102]
        expected_pc = pc.copy()
        expected_pc['entity_id'] = [101, 102]
        expected_ehr = ehr.copy()
        expected_ehr['entity_id'] = [101, 102]

        pd.testing.assert_frame_equal(mapped_prs, expected_prs)
        pd.testing.assert_frame_equal(mapped_pc, expected_pc)
        pd.testing.assert_frame_equal(mapped_ehr, expected_ehr)

    @patch('pipeline.genopred_adapter.GenopredAdapter._clean_profile')
    def test_prepare_prs(self, mock_read_csv):
        # Sample PRS data
        sample_prs = pd.DataFrame({
            'prs': [0.1, 0.5, 0.9, 1.2, 1.8]
        })

        # Initialize the GenopredAdapter class
        adapter = GenopredAdapter(profiles_file='mock_profiles_file', db_handler=MagicMock(), dp=MagicMock(), gwas_name='test_gwas')
        adapter.prs = sample_prs

        # Call _prepare_prs and verify result
        result = adapter._prepare_prs()

        # Verify quartile and decile columns are added
        assert f'test_gwas_quartile' in result.columns
        assert f'test_gwas_decile' in result.columns

    def test_prepare_logistic_reg_with_mock_db(self):
        # Sample input data for the test
        ehr = pd.DataFrame({
            'target_id': ['D1', 'D1', 'D2', 'D3'],
            'entity_id': [101, 102, 103, 104],
            'type': ['test_type', 'test_type', 'test_type', 'test_type']
        })

        prs = pd.DataFrame({
            'entity_id': [101, 102, 103, 104],
            'prs': [0.1, 0.5, 0.9, 1.2]
        })

        disdf = pd.DataFrame({
            'entity_id': [101, 102, 103, 104]
        })

        pc = pd.DataFrame({
            'entity_id': [101, 102, 103, 104],
            'PC1': [0.1, 0.2, 0.3, 0.4],
            'GENDER': ['Male', 'Male', 'Female', 'Male'],
            'CURRENT_AGE': [20, 30, 40, 50]
        })

        # Mocked target count data from the database
        mocked_target_count = pd.DataFrame({
            'target_id': ['D1', 'D2', 'D3'],
            'phenotype_count': [25, 15, 5]
        })

        # Expected final merged DataFrame (df)
        expected_df = pd.DataFrame({
            'entity_id': [101, 102, 103, 104],
            'D1': [1, 1, 0, 0],
            'prs': [0.1, 0.5, 0.9, 1.2],
            'PC1': [0.1, 0.2, 0.3, 0.4],
            'GENDER': ['Male', 'Male', 'Female', 'Male'],
            'CURRENT_AGE': [20, 30, 40, 50]
        })

        # Initialize the GenopredAdapter class
        adapter = GenopredAdapter.__new__(GenopredAdapter)  # Bypass __init__
        adapter.db_handler = MagicMock()  # Mock the db_handler

        # Mock the `get_target_counts` method to return the mocked target count
        adapter.db_handler.get_target_counts.return_value = mocked_target_count

        # Call the method with the mocked data
        df, target_count = adapter._prepare_logistic_reg(ehr, prs, disdf, pc)

        # Assert the final merged DataFrame (df) matches the expected output
        pd.testing.assert_frame_equal(df, expected_df)


class TestAdapterFactory:
    @patch('pipeline.genopred_adapter.GenopredAdapter._clean_profile')
    def test_adapter_factory_metabolite(self, mock_read_csv):
        general_adapter = AdapterFactory.create_general_adapter( 'mock_profiles_file', MagicMock(), MagicMock(), 'test_gwas')
        adapter = AdapterFactory.get_specific_adapter(general_adapter, 'Metabolite')
        assert isinstance(adapter, MetaboliteAdapter)

    @patch('pipeline.genopred_adapter.GenopredAdapter._clean_profile')
    def test_adapter_factory_disease(self, mock_read_csv):
        general_adapter = AdapterFactory.create_general_adapter( 'mock_profiles_file', MagicMock(), MagicMock(), 'test_gwas')
        adapter = AdapterFactory.get_specific_adapter(general_adapter, 'ICD code')
        assert isinstance(adapter, DiseaseCodesAdapter)

    @patch('pipeline.genopred_adapter.GenopredAdapter._clean_profile')
    def test_adapter_factory_quest_log(self, mock_read_csv):
        general_adapter = AdapterFactory.create_general_adapter( 'mock_profiles_file', MagicMock(), MagicMock(), 'test_gwas')
        adapter = AdapterFactory.get_specific_adapter(general_adapter, 'Questionaire_log')
        assert isinstance(adapter, QuestionaireLogAdapter)

    @patch('pipeline.genopred_adapter.GenopredAdapter._clean_profile')
    def test_adapter_factory_quest_lin(self, mock_read_csv):
        general_adapter = AdapterFactory.create_general_adapter( 'mock_profiles_file', MagicMock(), MagicMock(), 'test_gwas')
        adapter = AdapterFactory.get_specific_adapter(general_adapter, 'Questionaire_lin')
        assert isinstance(adapter, QuestionaireLinAdapter)

    @patch('pipeline.genopred_adapter.GenopredAdapter._clean_profile')
    def test_adapter_factory_invalid_type(self, mock_read_csv):
        with pytest.raises(ValueError):
            general_adapter = AdapterFactory.create_general_adapter('mock_profiles_file', MagicMock(), MagicMock(), 'test_gwas')
            adapter = AdapterFactory.get_specific_adapter(general_adapter,'Unknown')
