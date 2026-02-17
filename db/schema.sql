-- Detailed PheWAS / regression results (store regression metadata and statistics)
-- Expanded to store coefficients, SE, CIs, p-values, ORs and analysis metadata

-- Analysis run registry (one row per pipeline run)
CREATE TABLE IF NOT EXISTS analysis_run (
    run_id TEXT PRIMARY KEY,
    prs_name TEXT,
    label TEXT,
    n_groups INTEGER,
    include_intermediates BOOLEAN,
    normalize BOOLEAN,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_analysis_prs ON analysis_run(prs_name);

CREATE TABLE IF NOT EXISTS phewas_result (
    run_id TEXT,
    prs_name TEXT,
    target_code TEXT,

    -- main statistics
    odds_ratio REAL,
    ci_low REAL,
    ci_high REAL,
    p_value REAL,
    beta REAL,

    -- additional regression metadata
    SE REAL,
    CI_lower REAL,
    CI_upper REAL,
    n INTEGER,
    n_groups INTEGER,
    include_intermediates BOOLEAN,
    covariates TEXT,
    formula TEXT,
    sex_filter TEXT,

    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (run_id, target_code)
);

-- PRS percentiles / prevalence outputs (from compute_percentiles.py)
CREATE TABLE IF NOT EXISTS percentile_result (
    prs_name TEXT,
    prs_column TEXT,
    target_code TEXT,
    percentile INTEGER,
    sex TEXT,
    value REAL,
    n INTEGER,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (prs_name, prs_column, target_code, percentile, sex)
);

-- PRS manifest (from preprocessing/prs_present.csv)
CREATE TABLE IF NOT EXISTS prs_manifest (
    prs_name TEXT PRIMARY KEY,
    path TEXT,
    label TEXT,
    sex TEXT,
    full_path TEXT,
    discovered_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- GWAS metadata (ingested from data/gwas_metadata.csv)
CREATE TABLE IF NOT EXISTS gwas_metadata (
    name TEXT PRIMARY KEY,
    path TEXT,
    label TEXT,
    n_cases INTEGER,
    n_controls INTEGER,
    n INTEGER,
    population TEXT,
    sex TEXT,
    sampling TEXT,
    prevalence REAL,
    mean REAL,
    sd REAL,
    source TEXT,
    sumstats_source TEXT,
    prevalence_mean_source TEXT,
    comments TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_gwas_label ON gwas_metadata(label);

-- Legacy prevalence table (populated from percentile aggregates for compatibility)
CREATE TABLE IF NOT EXISTS prevalence (
    run_id TEXT,
    prs_name TEXT,
    prs_column TEXT,
    target_code TEXT,
    percentile INTEGER,
    sex TEXT,
    prevalence REAL,
    n INTEGER,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (run_id, prs_name, prs_column, target_code, percentile)
);

-- Phenotype / Target metadata populated during preprocessing
CREATE TABLE IF NOT EXISTS target (
    target_code TEXT PRIMARY KEY,
    description TEXT,
    target_class TEXT,
    class_file TEXT,
    domain TEXT,
    target_type TEXT,
    sex TEXT,
    covariates TEXT,
    full_path TEXT,
    file_exists BOOLEAN,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS disease_prevalence (
    target_code TEXT NOT NULL,
    target_class TEXT NOT NULL,
    sex TEXT NOT NULL CHECK (sex IN ('male', 'female', 'both')),
    prevalence REAL NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (target_code, sex)
);

CREATE TABLE IF NOT EXISTS cohort_distribution (
    target_code TEXT NOT NULL,
    variable TEXT NOT NULL,        -- age | bmi | hs
    category TEXT NOT NULL,        -- bin label
    sex TEXT NOT NULL,
    count INTEGER NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (target_code, variable, category, sex)
);

-- Indexes for fast lookup
CREATE INDEX IF NOT EXISTS idx_phewas_run ON phewas_result(run_id);
CREATE INDEX IF NOT EXISTS idx_phewas_prs ON phewas_result(prs_name);
CREATE INDEX IF NOT EXISTS idx_percentile_prs ON percentile_result(prs_name);
CREATE INDEX IF NOT EXISTS idx_prev_run ON prevalence(run_id);
-- Composite index to speed up lookups by PRS name and target and to help ORDER BY percentile
CREATE INDEX IF NOT EXISTS idx_prevalence_prs_target_percentile ON prevalence(prs_name, target_code, percentile);
CREATE INDEX IF NOT EXISTS idx_target_class ON target(target_class);
CREATE INDEX IF NOT EXISTS idx_prev_target ON disease_prevalence (target_code, sex);
CREATE INDEX IF NOT EXISTS idx_dist_target ON cohort_distribution (target_code, variable, sex);
