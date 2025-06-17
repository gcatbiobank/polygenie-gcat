# PolyGenie

**PolyGenie** is a scalable pipeline for phenome-wide association studies (PheWAS) leveraging polygenic risk scores (PRS). It processes GWAS summary statistics, individual-level genotype data, and diverse phenotypic inputs to uncover PRS-based associations across health records, metabolite levels, and survey data. Results are stored in a database and made explorable via an interactive web dashboard.

---

## Features

* Automated ingestion of GWAS summary statistics
* PRS computation using GenoPred and MegaPRS
* Phenome-wide association analysis (ICD, Phecodes, metabolites, traits)
* Logistic and linear regression with covariates
* False discovery rate (FDR) correction
* SQLite database for results storage
* Dash web app for interactive exploration

---

## Requirements

* Python 3.8+
* SQLite (for local development)
* Dependencies: `dash`, `plotly`, `pandas`, `numpy`, `sqlalchemy`, `scikit-learn`, etc.

Install via:

Environment requirements listed in environment.yml

---

## Input Configuration

All input paths are defined in `input_data/config.ini`. Core files include:

| Parameter                                              | Description                     |
| ------------------------------------------------------ | ------------------------------- |
| `metadata_file`                                        | GWAS metadata file (CSV)        |
| `metab_file`                                           | Metabolite values (TSV)         |
| `icd_file`                                             | Diagnoses using ICD codes (TSV) |
| `phecodes_file`                                        | Phecode mappings (TSV)          |
| `indiv_file`                                           | Individual information (TSV)    |
| `phenotypes_file`                                      | General phenotype file (TSV)    |
| `cohorts_file`, `population_file`, `pc_file`           | Cohort, ancestry, PCA data      |
| `questionaire_data_file`, `questionaire_metadata_file` | Survey data and codebook        |

For GWAS summary stats, include a `gwas_metadata.csv` file with:

| Column                          | Description                |
| ------------------------------- | -------------------------- |
| `name`                          | Trait ID                   |
| `path`                          | File path to summary stats |
| `label`                         | Display name               |
| `n`, `population`, `prevalence` | Metadata fields (optional) |

---

##  Pipeline Execution

Launch the main script:

```bash
bash polygenie.sh
```

This shows an interactive menu:

```
1) Launch Pipeline            # Run PRS computation and associations
2) Manage Database           # Update EHR, individuals, or pipeline data
3) Update PolyGenie Website  # Refresh the web interface
4) Exit
```

Each option walks through specific update or analysis actions. For example:

* `Launch Pipeline`: Runs full PRS → PheWAS pipeline
* `Manage Database`: Submenu to update individual, clinical, or PRS results

---

## Web Application

To serve the Dash web interface:

```bash
cd polygenie_app
gunicorn -w 4 -b 0.0.0.0:8050 app.app:server
```

Use a reverse proxy (e.g., nginx) and `systemd` for production deployments.

**Note on devices**: While functional on modern browsers (Chrome, Firefox, Edge), layout and performance may be limited on small screens or older browsers.

---

## Use Cases

* Explore disease prevalence by PRS percentile
* Visualize PRS-disease-metabolite links

---

## Project Structure

```
polygenie/
├── input_data/             # Input TSV and config files
├── app/          			# Dash web app
├── sqlitedb/               # Local DB and schema
├── scripts/                # Scripts that execute parts of the pipeline
├── pipeline/               # Pipeline code
├── profiles_data/          # PRS data
├── logs/               	# Pipeline logs
├── polygenie.sh            # Main control script
├── environment.yml
└── README.md
```

---

## Limitations

* Mobile responsiveness is limited in the Dash UI
* Associations depend on GWAS quality and cohort completeness

---

## License

MIT License. See `LICENSE` file.

---

## Contributing

Pull requests are welcome. Open an issue to propose changes or features.

---

## Contact

For questions or collaborations: \[[gcatbiobank@igtp.cat](gcatbiobank@igtp.cat)]


## Credits

- **Rafael de Cid Ibeas** - *Project Manager*
- **Xavier Farré Ramon** - *Project Supervisor*
- **Mireia Gasco Agorreta** - *Developer*
