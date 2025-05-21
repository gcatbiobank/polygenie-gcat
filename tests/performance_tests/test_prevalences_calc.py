import cProfile
import logging
import time
from sqlitedb.db_handler import DBHandler
from pipeline.prevalence_calculator import PrevalenceCalculator

# Logger for db_loader messages
logger = logging.getLogger('prevalences_logger')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('logs/prevalences_performance.log')
handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

def test_prevalences_profile():
    profiler = cProfile.Profile()
    profiler.enable()

    start = time.time()
    # Run the prevalence calculation
    calculator = PrevalenceCalculator()
    db_handler = DBHandler('sqlitedb/polygenie.db')
    prevalences = calculator.calculate_prevalences(db_handler, 'psoriasis')
    logger.info(f"Time elapsed: {time.time()-start}")

    profiler.disable()
    profiler.dump_stats('tests/test_results/prevalences_profiling_results_improved.prof')