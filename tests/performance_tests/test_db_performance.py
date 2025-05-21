import sqlite3
import time

from sqlitedb.db_handler import DBHandler


# Function to benchmark a query
def benchmark_correlations(iterations=5):
    db_handler = DBHandler('sqlitedb/polygenie.db')
    gwas_codes = db_handler.get_gwas_codes()
    references = ['low']
    divisions = ['quartile']
    target_types = ['Metabolite', 'ICD code', 'Phecode']

    total_time = 0
    for gwas in gwas_codes:
        for reference in references:
            for division in divisions:
                for target_type in target_types:
                    for _ in range(iterations):
                        start_time = time.time()
                        db_handler.get_correlations(gwas, reference, division, target_type)
                        total_time += time.time() - start_time
    return total_time / (iterations * len(gwas_codes) * len(references) * len(divisions) * len(target_types))

def benchmark_targets(iterations=5):
    db_handler = DBHandler('sqlitedb/polygenie.db')
    
    target_types = ['Metabolite', 'ICD code', 'Phecode']

    total_time = 0
    for target_type in target_types:
        for _ in range(iterations):
            start_time = time.time()
            db_handler.get_target_codes(target_type)
            total_time += time.time() - start_time
                        
    return total_time / (iterations * len(target_types))


def test_benchmark_index_correlations():
    # Connect to the SQLite database
    conn = sqlite3.connect('sqlitedb/polygenie.db')
    cursor = conn.cursor()

    # Step 1: Benchmark without the index
    cursor.execute("DROP INDEX IF EXISTS idx_correlations_gwas_ref_div")
    conn.commit()

    no_index_time = benchmark_correlations()
    with open('tests/test_results/index_performance_results.txt', 'w') as f:
        f.write(f"Average execution time without idx_correlations_gwas_ref_div index: {no_index_time:.6f} seconds\n")

    # Step 2: Add the index
    cursor.execute("CREATE INDEX idx_correlations_gwas_ref_div ON correlations (gwas, reference, division)")
    conn.commit()

    # Step 3: Benchmark with the index
    index_time = benchmark_correlations()
    with open('tests/test_results/index_performance_results.txt', 'a') as f:
        f.write(f"Average execution time with idx_correlations_gwas_ref_div index: {index_time:.6f} seconds\n")

    # Step 4: Close the connection
    conn.close()


def test_benchmark_index_targets():
    # Connect to the SQLite database
    conn = sqlite3.connect('sqlitedb/polygenie.db')
    cursor = conn.cursor()

    # Step 1: Benchmark without the index
    cursor.execute("DROP INDEX IF EXISTS idx_targets_type")
    conn.commit()

    no_index_time = benchmark_targets()
    with open('tests/test_results/index_performance_results.txt', 'a') as f:
        f.write(f"Average execution time without idx_targets_type index: {no_index_time:.6f} seconds\n")

    # Step 2: Add the index
    cursor.execute("CREATE INDEX idx_targets_type ON targets (type)")
    conn.commit()

    # Step 3: Benchmark with the index
    index_time = benchmark_correlations()
    with open('tests/test_results/index_performance_results.txt', 'a') as f:
        f.write(f"Average execution time with idx_targets_type index: {index_time:.6f} seconds\n")

    # Step 4: Close the connection
    conn.close()
