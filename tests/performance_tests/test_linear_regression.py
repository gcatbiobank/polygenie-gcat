from unittest.mock import Mock
import pytest
import time
import pickle
from pathlib import Path
from statistics import mean, stdev

from tests.performance_tests.correlation_calculator_tests import LinearCorrelationCalculatorCloud, LinearCorrelationCalculatorLocal

# Path to the pickled test data file
TEST_DATA_PATH = Path("tests/test_data/Metabolite_adapted_data.pkl")
RESULTS_FILE = Path("tests/test_results/performance_tests_results.txt")
OUTPUT_FILE = Path("tests/test_results/full_run_output.pkl")

@pytest.fixture
def test_data():
    """Fixture to load the pickled test data."""
    with open(TEST_DATA_PATH, 'rb') as f:
        data = pickle.load(f)
    return data[0][0], data[0][1]  # Returns df1, df2 from the nested list format

def test_compute_correlations_local_performance(test_data):
    """Test the performance of the compute_correlations method and log the results."""
    df1, df2 = test_data
    calculator = LinearCorrelationCalculatorLocal()

    gwas_name = "psoriasis"
    target_type = "Metabolite"
    
    # Create a mock for db_handler with get_target_scope method always returning "both"
    db_handler_mock = Mock()
    db_handler_mock.get_target_scope.return_value = 'both'


    # Measure the time taken to compute correlations
    for n in [1, 10, 20, 30, 40, 50]:
        durations = []
        for rep in range(5):    #For each, get average of 5 runs
            start_time = time.time()
            calculator.compute_correlations(df1, df2, gwas_name, target_type, db_handler_mock, n)
            durations.append(time.time() - start_time)

        # Save the performance results to a file
        with open(RESULTS_FILE, 'a') as file:
            file.write(f"Performance test - compute_correlations (local): {mean(durations):.4f} +/- {stdev(durations):.04f} seconds for {n} regressions\n")

def test_compute_correlations_cloud_performance(test_data):
    """Test the performance of the compute_correlations method and log the results."""
    df1, df2 = test_data
    calculator = LinearCorrelationCalculatorCloud()

    gwas_name = "psoriasis"
    target_type = "Metabolite"
    
    # Create a mock for db_handler with get_target_scope method always returning "both"
    db_handler_mock = Mock()
    db_handler_mock.get_target_scope.return_value = 'both'


    # Measure the time taken to compute correlations
    for n in [1, 10, 20, 30, 40, 50]:
        durations = []
        for rep in range(5):    #For each, get average of 5 runs
            start_time = time.time()
            calculator.compute_correlations(df1, df2, gwas_name, target_type, db_handler_mock, n, rep)
            durations.append(time.time() - start_time)

        # Save the performance results to a file
        with open(RESULTS_FILE, 'a') as file:
            file.write(f"Performance test - compute_correlations (local - lithops): {mean(durations):.4f} +/- {stdev(durations):.04f} seconds for {n} regressions\n")


def test_complete_gwas_run_cloud_performance(test_data):
    """Test the performance of the compute_correlations method and log the results."""
    df1, df2 = test_data
    calculator = LinearCorrelationCalculatorCloud()

    gwas_name = "psoriasis"
    target_type = "Metabolite"
    
    # Create a mock for db_handler with get_target_scope method always returning "both"
    db_handler_mock = Mock()
    db_handler_mock.get_target_scope.return_value = 'both'

    n = -1
    durations = []
    
    start_time = time.time()
    results = calculator.compute_correlations(df1, df2, gwas_name, target_type, db_handler_mock, n, 1)
    durations.append(time.time() - start_time)

    results.to_pickle(OUTPUT_FILE)

    # Save the performance results to a file
    n_display = "all" if n == -1 else n
    with open(RESULTS_FILE, 'a') as file:
        file.write(f"Performance test - compute_correlations (local - lithops) all metabolites: {mean(durations):.4f} +/- {stdev(durations):.04f} seconds for {n_display} regressions\n")
