"""
Module: correlation_calculator
------------------------------
This module defines an abstract base class and concrete implementations for regression calculations.
It provides a common interface for different types of regression analyses, such as logistic or linear regression.

Classes:
--------
- CorrelationCalculator: Abstract base class defining a common interface for correlation calculations.
- LinearCorrelationCalculator: Concrete class for performing linear regression calculations.
- LogisticCorrelationCalculator: Concrete class for performing logistic regression calculations.
- CorrelationCalculatorFactory: Factory class to instantiate the appropriate correlation calculator based on regression type.

Usage:
------
Subclasses of CorrelationCalculator must implement the `compute_correlations` and `_compute_regression` methods to perform specific regression analyses.

Example:
--------
# Instantiate a linear regression calculator
linear_calculator = CorrelationCalculatorFactory.get_calculator('linear')
results = linear_calculator.compute_correlations(ehr_data, gwas_data, 'gwas_name', 'Metabolite', db_handler)

# Instantiate a logistic regression calculator
logistic_calculator = CorrelationCalculatorFactory.get_calculator('logistic')
results = logistic_calculator.compute_correlations(discounts_data, features_data, 'gwas_name', 'Phecode', db_handler)
"""

from abc import ABC, abstractmethod
import re
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from pipeline.singleton_metaclass import SingletonMetaclass

class CorrelationCalculator(ABC, metaclass=SingletonMetaclass):
    """
    Abstract base class defining a common interface for correlation calculations.

    This class uses the Singleton design pattern to ensure that only one instance of the class exists.
    Subclasses must implement the `compute_correlations` and `_compute_regression` methods to perform
    specific regression analyses.

    Methods
    -------
    compute_correlations(*args, **kwargs):
        Public method to perform the correlations calculation. Subclasses must implement this method
        with the appropriate parameters for their specific regression type (e.g., logistic or linear regression).

    _compute_regression(*args, **kwargs):
        Protected method to perform the regression calculation. Subclasses must implement this method
        with the appropriate parameters for their specific regression type.
    """

    @abstractmethod
    def compute_correlations(self, *args, **kwargs):
        """
        Public method to perform the correlations calculation.

        Subclasses must implement this method with the appropriate parameters for their
        specific regression type (e.g., logistic or linear regression).

        Parameters
        ----------
        *args : tuple
            Positional arguments, vary by subclass.
        **kwargs : dict
            Keyword arguments, vary by subclass.

        Raises
        ------
        NotImplementedError
            If the method is not implemented by the subclass.
        """
        raise NotImplementedError("Subclasses must implement this method")
    
    @abstractmethod
    def _compute_regression(self, *args, **kwargs):
        """
        Protected method to perform the regression calculation.

        Subclasses must implement this method with the appropriate parameters for their
        specific regression type.

        Parameters
        ----------
        *args : tuple
            Positional arguments, vary by subclass.
        **kwargs : dict
            Keyword arguments, vary by subclass.

        Raises
        ------
        NotImplementedError
            If the method is not implemented by the subclass.
        """
        raise NotImplementedError("Subclasses must implement this method")
    
    @staticmethod
    def _process_results(fit, outcome, comps, gwas_name, division, reference, type):
        """
        Helper function to format the logistic regression results.
        """
        results = []
        params = fit.params
        conf = fit.conf_int()
        conf['Odds Ratio'] = params
        conf.columns = ['5%', '95%', 'Odds Ratio']
        conf = np.exp(conf)

        r = fit.rsquared if type == 'linear' else fit.prsquared

        for comp in comps:
            results.append({
                'outcome': outcome,
                'OR': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", 'Odds Ratio'],
                'CI_5': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", '5%'],
                'CI_95': conf.loc[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]", '95%'],
                'Quartile': comp,
                'reference': reference,
                'division': division,
                'score': gwas_name,
                'P': fit.pvalues[f"C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}'))[T.{comp}]"],
                'R2': r
            })
        return results

class LinearCorrelationCalculator(CorrelationCalculator):
    """
    Concrete class for performing linear regression calculations.

    This class implements the `compute_correlations` and `_compute_regression` methods
    to perform linear regression analyses.

    Methods
    -------
    compute_correlations(data, formula):
        Perform linear regression calculations on the provided data using the specified formula.

    _compute_regression(data, formula):
        Perform the actual linear regression calculation using statsmodels.
    """

    @staticmethod
    def clean_column_name(name):
        """Helper function to clean column names."""
        return re.sub(r'\W|^(?=\d)', '_', name)


    def compute_correlations(self, ehr, df, gwas_name, target_type, db_handler):
        """
        Perform linear regression calculations on the provided data using the specified formula.

        Parameters
        ----------
        data : pandas.DataFrame
            The data on which to perform the regression.
        formula : str
            The formula specifying the regression model.

        Returns
        -------
        results : statsmodels.regression.linear_model.RegressionResultsWrapper
            The results of the linear regression.
        """
        res = []

        merged_df = pd.merge(ehr, df, on='entity_id')  

        # Ensure categorical columns are of type 'category'
        categorical_columns = ['GENDER', 'Chylomicrons']
        for column in categorical_columns:
            merged_df[column] = merged_df[column].astype('category')
        
        if target_type == "Metabolite": variables = ehr.columns[1:-15]
        else: variables = ehr.columns[1:-14]

        # Check all the targets are valid
        if target_type == 'Questionaire_lin':
            db_targets = set(db_handler.get_target_codes('continuous'))
        else:
            db_targets = set(db_handler.get_target_codes(target_type))
        variables = [var for var in variables if var in db_targets]

        # Parallel tasks
        tasks = [(metabolite, division, reference) 
            for metabolite in variables
            for division in ['quartile', 'decile']
            for reference in ['low', 'rest']]
        
        results = Parallel(n_jobs=-1)(delayed(self._compute_regression)(metabolite, division, reference, merged_df, gwas_name, db_handler.get_target_scope(metabolite)) for metabolite, division, reference in tasks)
        flattened_results = [d for sublist1 in results for sublist2 in sublist1 for d in sublist2 if d is not None]     
        res = pd.DataFrame(flattened_results)
        res['logpxdir'] = np.where(res['OR'] <1, np.log10(res['P']), -np.log10(res['P']))
        res.rename(columns={'outcome': 'target', 'score': 'gwas'}, inplace=True)

        return res

    def _compute_regression(self, metabolite, division, reference, merged_df, gwas_name, target_scope):
        """
        Perform the actual linear regression calculation using statsmodels.

        Parameters
        ----------
        data : pandas.DataFrame
            The data on which to perform the regression.
        formula : str
            The formula specifying the regression model.

        Returns
        -------
        results : statsmodels.regression.linear_model.RegressionResultsWrapper
            The results of the linear regression.
        """

        results = []
        local_df = merged_df.copy()

        metabolite_clean = self.clean_column_name(metabolite)

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
            comps = ['high', 'intermediate'] if reference == 'low' else ['high']
            results.append(self._process_results(fit, metabolite, comps, gwas_name, division, reference, 'linear'))
            return results
        
        except Exception as e:
            print(f"Error fitting model for {metabolite}: {e}")        
            return []    

class LogisticCorrelationCalculator(CorrelationCalculator):
    """
    Concrete class for performing logistic regression calculations.

    This class implements the `compute_correlations` and `_compute_regression` methods
    to perform logistic regression analyses.

    Methods
    -------
    compute_correlations(data, formula):
        Perform logistic regression calculations on the provided data using the specified formula.

    _compute_regression(data, formula):
        Perform the actual logistic regression calculation using statsmodels.
    """

    def compute_correlations(self, df, discounts, gwas_name, target_type, db_handler):
        """
        Perform logistic regression calculations on the provided data using the specified formula.

        Parameters
        ----------
        data : pandas.DataFrame
            The data on which to perform the regression.
        formula : str
            The formula specifying the regression model.

        Returns
        -------
        results : statsmodels.discrete.discrete_model.BinaryResultsWrapper
            The results of the logistic regression.
        """
        res = []

        if target_type == 'Phecode':
            discounts['target_id'] = discounts['target_id'].apply(lambda x: 'A' + str(x).replace('.', 'B'))
            cols = df.columns.tolist()
            middle_cols = cols[1:-19]
            df.rename(columns={col: 'A' + str(col).replace('.', 'B') for col in middle_cols}, inplace=True)

        # Parallel tasks
        tasks = [(outcome, division, reference) 
            for outcome in discounts['target_id'].unique()
            for division in ['quartile', 'decile']
            for reference in ['low', 'rest']]
        
        results = Parallel(n_jobs=-1)(delayed(self._compute_regression)(outcome, division, reference, df, gwas_name, db_handler.get_target_scope(outcome)) for outcome, division, reference in tasks)
        flattened_results = [d for sublist1 in results for sublist2 in sublist1 for d in sublist2 if d is not None]

        # Finish formatting
        res = pd.DataFrame(flattened_results)
        res['logpxdir'] = np.where(res['OR'] <1, np.log10(res['P']), -np.log10(res['P']))
        res['gwas'] = gwas_name
        res.rename(columns={'outcome':'target'}, inplace=True)
        
        if (target_type == 'Phecode'):
            res['target'] = res['target'].apply(lambda x: x.lstrip('A').replace('B', '.')) # Undo the changes to the phecodes

        return res

    def _compute_regression(self, outcome, division, reference, df, gwas_name, target_scope):
        """
        Perform the actual logistic regression calculation using statsmodels.

        Parameters
        ----------
        data : pandas.DataFrame
            The data on which to perform the regression.
        formula : str
            The formula specifying the regression model.

        Returns
        -------
        results : statsmodels.discrete.discrete_model.BinaryResultsWrapper
            The results of the logistic regression.
        """

        results = []

        
        if target_scope == 'female':
            df = df[df['GENDER'] == 'female']
        elif target_scope == 'male':
            df = df[df['GENDER'] == 'male']
        else:
            pass
        
        try:
            fit = smf.logit(formula=f"{outcome} ~ C({gwas_name}_{division}_{reference}, Treatment(reference='{reference}')) + CURRENT_AGE + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", data=df).fit(method='bfgs', maxiter=300, disp=False)
            comps = ['high', 'intermediate'] if reference == 'low' else ['high']
            results.append(self._process_results(fit, outcome, comps, gwas_name, division, reference, 'logistic'))
            return results
        
        except Exception as e:
            print(f"Error in compute_logit: {e}")
            return []
    
class CorrelationCalculatorFactory:
    """
    Factory class to instantiate the appropriate correlation calculator based on regression type.

    Methods
    -------
    get_calculator(regression_type):
        Return an instance of the appropriate correlation calculator based on the specified regression type.
    """

    @staticmethod
    def get_calculator(regression_type: str) -> CorrelationCalculator:
        """
        Return an instance of the appropriate correlation calculator based on the specified regression type.

        Parameters
        ----------
        regression_type : str
            The type of regression ('linear' or 'logistic').

        Returns
        -------
        CorrelationCalculator
            An instance of the appropriate correlation calculator.

        Raises
        ------
        ValueError
            If the specified regression type is not supported.
        """
        if regression_type == 'logistic':
            return LogisticCorrelationCalculator()
        elif regression_type == 'linear':
            return LinearCorrelationCalculator()
        else:
            raise ValueError(f"Unsupported regression type: {regression_type}")
