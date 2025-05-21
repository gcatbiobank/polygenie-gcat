# PolyGenie - Technical Documentation

This documentation contains a detailed specification of the PolyGenie pipeline and database structures, providing insights to the 
most relevant design decisions and the formats expected for the input files required for each script.

## Table of Contents
1. [Overview](#overview)
2. [Pipeline Description](#pipeline-description)
    1. [Configuration to add a new GWAS](#configuration-to-add-a-new-gwas)
    2. [Metadata format](#metadata-format)
    3. [Genotypes format](#genotype-data-format)
    4. [Execution Instructions](#execution-instructions)
3. [Database Description](#database-description)
    1. [EER Model](#eer-model)
    2. [Logical Model](#logical-model)
    3. [Physical Model - SQLite](#physical-model---sqlite)
    4. [How to Use the db_manager](#how-to-use-the-db_manager)
4. [Web deployment](#web-deployment)
5. [Credits](#credits)

## Overview
The following graph represents the directory structure of the project and also the main files of the GenoPred pipeline structure.

```bash
Project
 |
 |---- launch_pipeline.py
 |---- deploy_app.py
 |
 |---- pipeline
 |         |
 |         |---- load_data.py
 |         |---- data_processing.py
 |         |---- gwas_list_setup.py
 |         |---- header_check.py
 |         |---- logs
 |         |---- data
 |
 |------- app
 |         |         |        
 |         |---- metadata
 |         |       |
 |         |       |---- gwas_metadata.csv
 |         |     
 |         |---- assets
 |         |       |
 |         |       |---- style.css
 |         |
 |         |---- app.py
 |
 |---- sqlitedb
 |         |
 |         |---- db_upload_phe_ind_targ.py
 |         |---- db_upload_pipeline_data.py
 |         |---- db_handler.py
 |         |---- db_manager.sh
 |         |---- schema.sql
 |         |---- polygenie.db
 |
 |---- tests
 |------- environment.yml

GenoPred
 |
 |---- pipeline
 |
 |---- input_gcat
 |        |
 |        |---- config.yaml
 |        |---- gwas_list.yaml
 |        |---- target_list.yaml
 |
 |---- results_gcat
 |        |
 |        |---- output
 |        |---- gwas_data

```
This structure contains the three main modules of this software: the pipeline to obtain the correlations data, the database to store those results and the dashboard app to show them.  We will now focus on these three blocks.


## Pipeline Description
The pipeline is configured to do the following process:
1. **Read the metadata file**: check the new GWAS entries, marked as status=P (pending), copy them to the gwas_list.txt from the GenoPred pipeline and marke them as status=C (completed). Command:
```sh
python pipeline/gwas_list_setup.py
```
2. **GenoPred dry run**: ensure that the pipeline detects the new GWAS, abort if not. Command:
```sh
output=$(snakemake -n --configfile=input_gcat/config.yaml --use-conda ${OUTPUT})
```
3. **Header check**: ensure that the headers of the GWAS will be correctly detected by GenoPred.  In case some mismatch is detected, the pipeline will show an error and stop execution. Command:
```sh
python pipeline/header_check.py
```
4. **GenoPred**: launch the GenoPred pipeline with the number of cores indicated as a parameter to the script. Command:
```sh
snakemake -j${NUM_CORES} --configfile=input_gcat/config.yaml --use-conda ${OUTPUT}
```
5. **Data processing**: process the output data to generate the information that will be displayed in the dashboard.  That information is stored in the   data directory of the app. Command:
```sh
python pipeline/load_data.py ../GenoPred/pipeline/results_gcat/output/GCATcore/pgs/EUR/megaprs
```
6. **Database update**: upload newly computed data to the database, which will make it visible on the dashboard app. Command:
```sh
python sqlitedb/db_upload_pipeline_data.py
```


### Configuration to Add a New GWAS
To compute the PRS for a new GWAS, some files and modifications are required.  They are listed here to simplify the process:
1. Add the GWAS .gz file to: `/imppc/labs/dnalab/studentdnabank/GenoPred/pipeline/results_gcat/gwas_data`
2. Add the GWAS metadata to the `gwas_metadata.csv` in: `/imppc/labs/dnalab/studentdnabank/Project/app/metadata/gwas_metadata.csv`. The details of the metadata file
    can be checked in the [Metadata Format](#metadata-format) chapter
3. Launch the pipeline with the `launch_pipeline.sh` script, indicating the number of cores you want to use and the output you want to get from the GenoPred pipeline.


### Metadata Format
You can check here the information required by each field:

| Field Name             | Description                                                                                                 |
|------------------------|-------------------------------------------------------------------------------------------------------------|
| **name**               | Name of the condition linked to the GWAS.                                                                   |
| **path**               | Path to the sumstats file.                                                                                  |
| **date**               | Date when the GWAS was downloaded.                                                                          |
| **label**              | A human readable name for the GWAS phenotype. Wrap in double quotes if multiple words.                      |
| **n_cases**            | Number of cases.                                                                                            |
| **n_controls**         | Number of controls.                                                                                         |
| **n**                  | The total sample size of the GWAS.                                                                          |
| **population**         | Reference population that the GWAS sample matches best.                                                     |
| **sampling**           | The proportion of the GWAS sample that were cases (if outcome is binary - otherwise specify NA)             |
| **prevalence**         | The prevalence of the phenotype in the general population (if outcome is binary - otherwise specify NA)     |
| **mean**               | The phenotype mean in the general population (if outcome is continuous, otherwise specify NA)               |
| **sd**                 | The phenotype sd in the general population (if outcome is continuous, otherwise specify NA)                 |
| **source**             | Link to the GWAS publication.                                                                               |
| **sumstat_source**     | Link to the GWAS summary statistics data.                                                                   |
| **prevalence_mean_source** | Link to the source of the prevalence and mean information. This is regarding the general population.    |


### Genotype Data Format
Genotype data is provided in plink2 format and divided by chromosomes.  This is the format set in the GenoPred pipeline configuration, but it can be changed to plink1 or other formats in case of need, modifying the [configuration file](../GenoPred/pipeline/input_gcat/target_list.txt).  The path to the target genotype data can also be modified in the same file.


### Execution Instructions
To correctly run the pipeline using the `launch_pipeline.sh`, two parameters need to be set: the number of cores you want to use for the tasks computed in parallel and 
the desired output you want to get from the GenoPred pipeline.  For the latter, it is recommended to use `target_pgs`, as it is the configuration that has been used
during the development of this software.  Other configurations have not beed tested and could lead to unexpected results.

Example command to run the pipeline:

```sh
./launch_pipeline.sh 8 target_pgs
```


## Database Description

### EER Model
#TODO


### Logical Model
The EER model detailed above translates to the following logical model:

```
GWAS(code, name, date, link_paper, link_sumstats, link_prevalence, n, population) 
        with    {code} as PK
                {name} as AK
        where   {population} references POPULATIONS
            
TARGETS(code, description, class, type) 
        with    {code} as PK

ICD_CODES(code) 
        with    {code} as PK
        where   {code} references TARGETS 

PHECODES(code) 
        with    {code} as PK
        where   {code} references TARGETS 

METABOLITES(code) 
        with    {code} as PK
        where   {code} references TARGETS 

INDIVIDUALS(iid, entity_id, gender, age, bmi, self_perceived_hs) 
        with    {iid} as PK
                {entity_id} as AK

DIVISIONS(type) 
        with    {type} as PK

REFERENCES(type) 
        with    {type} as PK

PERCENTILES(number) 
        with    {number} as PK

CORRELATIONS(gwas, target, reference, division, odds_ratio, CI_5, CI_95, P, R2, logpxdir)
        with    {gwas, target, reference, division} as PK
        where   {gwas} references GWAS 
                {target} references TARGETS 
                {reference} references REFERENCES 
                {division} references DIVISION 

PHENOTYPES(indiv_id, target_id) 
        with    {indiv_id, target_id} as PK
        where   {indiv_id} references INDIVIDUALS 
                {target_id} references TARGETS 

PRS(indiv_id, gwas_id, prs_score, prs_percentile_all, prs_percentile_female, prs_percentile_male) 
        with    {indiv_id, gwas_id} as PK
        where   {indiv_id} references INDIVIDUALS 
                {gwas_id} references GWAS 

PREVALENCES(gwas_id, target_id, percentile, prevalence_all, prevalence_female, prevalence_male) 
        with    {gwas_id, target_id, percentile} as PK
        where   {gwas_id} references GWAS 
                {target_id} references TARGETS 
                {percentile} references PERCENTILES 
```

### Physical Model - SQLite
The final schema of the DB can be found in the [schema.sql](/sqlitedb/schema.sql) script. The DB has been implemented in SQLite for simplicity, as it allows for easy backups and does not require server maintenance. This implementation is expected to work well while the system is only read-based, and the writes are only performed locally using the `db_manager.sh script`.
In case that the system evolves to be more write-heavy, especially if concurrent writes are expected, migrating to a more complex DBMS such as MySQL or PostreSQL is recommended. Given that all the logic for the DB is managed from the [DBHandler](/sqlitedb/db_handler.py), such migration should be relatively easy, as only the database schema and the DBHandler would need to be modified.


### How to Use the db_manager
#TODO


## Web Deployment
The version of the app on production runs on a Debian 12 based virtual machine (VM), and is configured with `Gunicorn` as the WSGI server and `Nginx` as the proxy and load balancer.  `Systemd` is used to manage the app process, to ensure it is restarted in case it crashes or the server stops for external reasons.

Each time new information is computed, the server needs to be updated to show it.  The steps are the following:
 
 1. Log into the VM and stop the running service
 ```sh
 ssh user@address
 sudo systemctl stop polygenie_app.service
 exit
 ```
 2. Update the data on the server running the deploy_app.sh script
 ```sh
 ./deploy_app.sh
 ```
 3. Log into the VM again and restart the service
 ```sh
 ssh user@address
 sudo systemctl start polygenie_app
 sudo systemctl enable polygenie_app
 exit
 ```

 ## Credits

- **Rafael de Cid Ibeas** - *Project Manager*
- **Xavier Farré Ramon** - *Project Supervisor*
- **Mireia Gasco Agorreta** - *Developer*
- **Gabriel Plana Gavaldà** - *DB Design Supervisor*
- **Marc Vidal Hernàndez** - *Graphic Design Supervisor*
