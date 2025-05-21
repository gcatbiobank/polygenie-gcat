import os
import time
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
import re
import statsmodels.formula.api as smf
import lithops
from lithops.storage import Storage
import cloudpickle as pickle
from pipeline.correlation_calculator import CorrelationCalculator


class LinearCorrelationCalculatorLocal(CorrelationCalculator):
    """
    Linear regression calculator. This class will be a Singleton.
    """

    @staticmethod
    def clean_column_name(name):
        """Helper function to clean column names."""
        return re.sub(r'\W|^(?=\d)', '_', name)


    def compute_correlations(self, ehr, df, gwas_name, target_type, db_handler, n):
        """
        Computes linear regression for each metabolite in the dataset, using a parallel processing strategy.

        Parameters
        ----------
        ehr : DataFrame
            A DataFrame containing electronic health records data. Check the technical documentation for more information on the format.
        df : DataFrame
            A DataFrame containing GWAS and pc data. Check the technical documentation for more information on the format.
        gwas_name : str
            The name of the GWAS score.

        Returns
        -------
        DataFrame
            A DataFrame containing the linear regression results. Check the technical documentation for more information on the format.
        """
        res = []

        merged_df = pd.merge(ehr, df, on='entity_id')  
        exclude = ['IID', 'entity_id']

        # Ensure categorical columns are of type 'category'
        categorical_columns = ['GENDER', 'Chylomicrons']
        for column in categorical_columns:
            merged_df[column] = merged_df[column].astype('category')

        # Parallel tasks
        first_metabolite = next(iter(ehr.columns.difference(exclude)))
        first_division = next(iter(['quartile', 'decile']))
        first_reference = next(iter(['low', 'rest']))

        tasks = [(first_metabolite, first_division, first_reference)] * n
        
        results = Parallel(n_jobs=-1)(delayed(self._compute_regression)(metabolite, division, reference, merged_df, gwas_name, db_handler.get_target_scope(metabolite)) for metabolite, division, reference in tasks)
        flattened_results = [item for sublist in results for item in sublist]        
        res = pd.DataFrame(flattened_results)
        res['logpxdir'] = np.where(res['OR'] <1, np.log10(res['P']), -np.log10(res['P']))
        res.rename(columns={'outcome': 'target', 'score': 'gwas'}, inplace=True)
        
        return res

    def _compute_regression(self, metabolite, division, reference, merged_df, gwas_name, target_scope):
        """
        Performs linear regression for the given metabolite with specified division and reference.
        This function is only for internal use.
        """
        results = []
        local_df = merged_df.copy()

        metabolite_clean = self.clean_column_name(metabolite)

        if metabolite != metabolite_clean:
            local_df.rename(columns={metabolite: metabolite_clean}, inplace=True)

        # Ensure metabolite column is numeric
        local_df[metabolite_clean] = pd.to_numeric(local_df[metabolite_clean], errors='coerce')

        if target_scope == 'female':
            merged_df = merged_df[merged_df['GENDER'] == 'female']
        elif target_scope == 'male':
            merged_df = merged_df[merged_df['GENDER'] == 'male']
        else:
            pass

        # Prepare the formula
        formula = f"{metabolite_clean} ~ C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}')) + CURRENT_AGE + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Chylomicrons"

        # Fit the linear regression model
        try:
            start_time = time.time()
            model = smf.ols(formula=formula, data=local_df)
            fit = model.fit()
            end_time = time.time() - start_time
            # Extract the parameters, confidence intervals, and odds ratios
            params = fit.params
            conf = fit.conf_int()
            conf['Odds Ratio'] = params
            conf.columns = ['5%', '95%', 'Odds Ratio']
            conf = np.exp(conf)
            if reference == 'low':
                comps = ['high', 'intermediate']
            else:
                comps = ['high']
            for comp in comps:
                    results.append(self._process_results(fit, comp, gwas_name, division, reference, metabolite, conf))
            
            return results
        except Exception as e:
            print(f"Error fitting model for {metabolite}: {e}")        
            return []    
        
    
    def _process_results(self, fit, comp, gwas_name, division, reference, metabolite, conf):
        """
        Helper method to process the regression results.
        """
        return {
            'outcome': metabolite,
            'OR': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", 'Odds Ratio'],
            'CI_5': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", '5%'],
            'CI_95': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", '95%'],
            'Quartile': comp,
            'reference': reference,
            'division': division,
            'score': gwas_name,
            'P': fit.pvalues[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]"],
            'R2': fit.rsquared
        }

lithops_config = {
    'lithops': {
        'storage': 'minio',
        'backend': 'k8s'
    },
    'minio': {
        'endpoint': 'http://localhost:44444',
        'access_key_id': 'lab144',
        'secret_access_key': 'astl1a4b4'
    }
}

storage = lithops.Storage(config=lithops_config)
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

class LinearCorrelationCalculatorCloud(CorrelationCalculator):
    """
    Linear regression calculator. This class will be a Singleton.
    """

    @staticmethod
    def clean_column_name(name):
        """Helper function to clean column names."""
        return re.sub(r'\W|^(?=\d)', '_', name)


    def compute_correlations(self, ehr, df, gwas_name, target_type, db_handler, n, rep):
        """
        Computes linear regression for each metabolite in the dataset, using a parallel processing strategy.

        Parameters
        ----------
        ehr : DataFrame
            A DataFrame containing electronic health records data. Check the technical documentation for more information on the format.
        df : DataFrame
            A DataFrame containing GWAS and pc data. Check the technical documentation for more information on the format.
        gwas_name : str
            The name of the GWAS score.

        Returns
        -------
        DataFrame
            A DataFrame containing the linear regression results. Check the technical documentation for more information on the format.
        """

        res = []

        merged_df = pd.merge(ehr, df, on='entity_id')  

        if target_type == "Metabolite": variables = ehr.columns[1:-15]
        else: variables = ehr.columns[1:-14]

        # Ensure categorical columns are of type 'category'
        categorical_columns = ['GENDER', 'Chylomicrons']
        for column in categorical_columns:
            merged_df[column] = merged_df[column].astype('category')

        serialized_df = pickle.dumps(merged_df)
        merged_df_cobj = storage.put_cloudobject(serialized_df)

        # Parallel tasks
        if n != -1:
            first_metabolite = next(iter(variables))
            first_division = next(iter(['quartile', 'decile']))
            first_reference = next(iter(['low', 'rest']))
            args = [(first_metabolite, first_division, first_reference, merged_df_cobj, gwas_name, 'both')] * n
        else:
            args = [(metabolite, division, reference, merged_df_cobj, gwas_name, 'both') 
                for metabolite in variables
                for division in ['quartile', 'decile']
                for reference in ['low', 'rest']]

        fexec = lithops.FunctionExecutor()
        futures = fexec.map(LinearCorrelationCalculatorCloud._compute_regression, args)
        results = fexec.get_result(futures)

        if rep == 1:
            fexec.plot(dst='tests/test_results/plots')

        flattened_results = [item for sublist in results for item in sublist]        
        res = pd.DataFrame(flattened_results)
        res['logpxdir'] = np.where(res['OR'] <1, np.log10(res['P']), -np.log10(res['P']))
        res.rename(columns={'outcome': 'target', 'score': 'gwas'}, inplace=True)
        
        return res

    @staticmethod
    def _compute_regression(metabolite, division, reference, merged_df_cobj, gwas_name, target_scope):
        """
        Performs linear regression for the given metabolite with specified division and reference.
        This function is only for internal use.
        """
        results = []

        serialized_df = storage.get_cloudobject(merged_df_cobj)
        merged_df = pickle.loads(serialized_df)

        local_df = merged_df.copy()

        metabolite_clean = LinearCorrelationCalculatorCloud.clean_column_name(metabolite)

        if metabolite != metabolite_clean:
            local_df.rename(columns={metabolite: metabolite_clean}, inplace=True)

        # Ensure metabolite column is numeric
        local_df[metabolite_clean] = pd.to_numeric(local_df[metabolite_clean], errors='coerce')

        if target_scope == 'female':
            local_df = local_df[local_df['GENDER'] == 'female']
        elif target_scope == 'male':
            local_df = local_df[local_df['GENDER'] == 'male']
        else:
            pass

        # Prepare the formula
        formula = f"{metabolite_clean} ~ C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}')) + CURRENT_AGE + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Chylomicrons"

        # Fit the linear regression model
        try:
            model = smf.ols(formula=formula, data=local_df)
            fit = model.fit()
            # Extract the parameters, confidence intervals, and odds ratios
            params = fit.params
            conf = fit.conf_int()
            conf['Odds Ratio'] = params
            conf.columns = ['5%', '95%', 'Odds Ratio']
            conf = np.exp(conf)
            if reference == 'low':
                comps = ['high', 'intermediate']
            else:
                comps = ['high']
            for comp in comps:
                    results.append({
                        'outcome': metabolite,
                        'OR': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", 'Odds Ratio'],
                        'CI_5': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", '5%'],
                        'CI_95': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", '95%'],
                        'Quartile': comp,
                        'reference': reference,
                        'division': division,
                        'score': gwas_name,
                        'P': fit.pvalues[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]"],
                        'R2': fit.rsquared
                    })
            return results
        except Exception as e:
            print(f"Error fitting model for {metabolite}: {e}")        
            return []