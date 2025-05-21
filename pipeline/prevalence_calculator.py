import itertools
from multiprocessing import Manager
from joblib import Parallel, delayed
import pandas as pd

from pipeline.singleton_metaclass import SingletonMetaclass

class PrevalenceCalculator(metaclass=SingletonMetaclass):
    """
    A class used to calculate prevalences of phenotypes based on polygenic risk scores (PRS).
    
    Methods
    -------
    calculate_prevalences(db_handler, gwas_name)
        Calculates the prevalences of phenotypes for a given GWAS name using data from the database handler.
    _process_batch(batch, prs, phenotypes)
        Processes a batch of tasks to compute prevalences for each GWAS and target combination.
    _compute_prevalences(gwas, target, ori_prs, phenotypes)
        Computes the prevalences of a target phenotype for each percentile of PRS, stratified by gender.
    """

    def calculate_prevalences(self, db_handler, gwas_name):
        """
        Calculate the prevalences of various phenotypes based on the provided GWAS name.
        This method retrieves target codes and phenotypes from the database handler, 
        generates tasks for each target, processes them in parallel, and returns the 
        calculated prevalences.
        Args:
            db_handler (DatabaseHandler): An instance of the database handler to interact with the database.
            gwas_name (str): The name of the GWAS (Genome-Wide Association Study) to be used for calculating prevalences.
        Returns:
            list: A list of calculated prevalences for the given GWAS name.
        """
        targets_list = db_handler.get_target_codes('Phecode') + \
                        db_handler.get_target_codes('ICD code') + \
                        db_handler.get_target_codes('binary')

        prs = db_handler.get_all_risk_scores()
        consulta = db_handler.get_all_indiv()
        gen2ent = dict(zip(consulta['iid'], consulta['entity_id']))
        prs['entity_id'] = prs['indiv_id'].map(gen2ent)
        phenotypes = pd.concat([
            db_handler.get_phenotypes('ICD code'),
            db_handler.get_phenotypes('Phecode')
        ], axis=0)

        # Generate all tasks
        batch_size = 100
        tasks = [(gwas_name, target) for target in targets_list]
        task_batches = [tasks[i:i + batch_size] for i in range(0, len(tasks), batch_size)]

        # Process batches in parallel
        results = Parallel(n_jobs=-1)(
            delayed(self._process_batch)(batch, prs, phenotypes) for batch in task_batches
        )

        # Flatten results
        prevalences = list(itertools.chain.from_iterable(results))
        return prevalences

    def _process_batch(self, batch, prs, phenotypes):
        """
        Processes a batch of GWAS and target pairs to compute prevalences.

        Args:
            batch (list of tuples): A list of tuples where each tuple contains a GWAS ID and a target ID.
            prs (pandas.DataFrame): A DataFrame containing PRS data with a 'gwas_id' column.
            phenotypes (pandas.DataFrame): A DataFrame containing phenotype data with a 'target_id' column.

        Returns:
            list: A list of computed prevalences for the given batch.
        """
        results = []
        for gwas, target in batch:
            results.extend(self._compute_prevalences(gwas, target, prs[prs['gwas_id'] == gwas], phenotypes[phenotypes['target_id'] == target]))
        return results
    
    def _compute_prevalences(self, gwas, target, ori_prs, phenotypes):
        """
        Compute the prevalence of a target phenotype across different percentiles of PRS (Polygenic Risk Scores).

        Args:
            gwas (str): The identifier for the GWAS (Genome-Wide Association Study) dataset.
            target (str): The target phenotype to compute prevalence for.
            ori_prs (pd.DataFrame): Original DataFrame containing PRS information for individuals.
            phenotypes (pd.DataFrame): DataFrame containing phenotype information for individuals.

        Returns:
            list: A list of tuples, each containing:
                - gwas (str): The identifier for the GWAS dataset.
                - target (str): The target phenotype.
                - percentile (int): The PRS percentile.
                - percentage_with_target['all'] (float): The percentage of individuals with the target phenotype in the given percentile for all genders.
                - percentage_with_target['female'] (float): The percentage of individuals with the target phenotype in the given percentile for females.
                - percentage_with_target['male'] (float): The percentage of individuals with the target phenotype in the given percentile for males.
        """
        prs = ori_prs.copy()
        prs.loc[:, 'has_target'] = prs['entity_id'].isin(phenotypes['indiv_id'])
        print(prs)

        results = []
        for percentile in range(100):
            percentage_with_target = {}
            for gender in ['all', 'male', 'female']:
                df_all = prs[(prs[f'prs_percentile_{gender}'] == percentile) & (prs[f'prs_percentile_{gender}'].notna())]
                total_all = len(df_all)
                total_with_target_all = int(df_all['has_target'].sum())
                percentage_with_target[gender] = (total_with_target_all / total_all) * 100 if total_all > 0 else 0

            results.append((
                gwas, 
                target, 
                percentile,
                percentage_with_target['all'], 
                percentage_with_target['female'], 
                percentage_with_target['male']
            ))
        return results
