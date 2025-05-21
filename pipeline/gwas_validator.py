import os
import sys
import gzip
import pandas as pd
from colorama import Fore, Style

class GWASValidator:
    """
    GWASValidator is a class designed to validate Genome-Wide Association Study (GWAS) summary statistics files. 
    It ensures that the required columns are present in the data files and logs the validation results.

    Attributes:
        ss_head_dict (dict): A dictionary mapping standard column names to their possible aliases.
        required_columns (dict): A dictionary specifying the required columns for validation.
        log_file (str): Path to the log file where validation results are recorded.
        reverse_dict (dict): A dictionary mapping aliases to their standard column names.

    Methods:
        __init__(log_file="pipeline/logs/header_check.log"):
            Initializes the GWASValidator with a log file and creates a reverse dictionary for column name mapping.
        
        create_reverse_dict(ss_head_dict):
            Creates a reverse dictionary mapping aliases to standard column names.
        
        map_columns(column_names):
            Maps a list of column names to their standard names using the reverse dictionary.
        
        check_required_columns(data, file_name):
            Checks if the provided data contains the required columns and logs the results.
        
        read_file(file_path):
            Reads a file into a DataFrame, handling different file formats and delimiters.
        
        validate_files(files_df):
            Validates multiple GWAS files and logs the results, returning a list of failed files.
    """
    ss_head_dict = {
        "SNP": [
            "SNP", "PREDICTOR", "SNPID", "MARKERNAME", "MARKER_NAME", "SNPTESTID",
            "ID_DBSNP49", "RSID", "ID", "RS_NUMBER", "MARKER", "RS", "RSNUMBER",
            "RS_NUMBERS", "SNP.NAME", "SNP ID", "SNP_ID", "LOCATIONALID", "ASSAY_NAME",
            "VARIANT_ID", "HM_RSID"
        ],
        "A1": [
            "A1", "ALLELE1", "ALLELE_1", "INC_ALLELE", "EA", "A1_EFFECT", "REF",
            "EFFECT_ALLELE", "RISK_ALLELE", "EFFECTALLELE", "EFFECT_ALL", "REFERENCE_ALLELE",
            "REF_ALLELE", "REFERENCEALLELE", "EA", "ALLELE_1", "INC_ALLELE", "ALLELE1",
            "A", "A_1", "CODED_ALLELE", "TESTED_ALLELE", "HM_EFFECT_ALLELE"
        ],
        "A2": [
            "A2", "ALLELE2", "ALLELE_2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "DEC_ALLELE",
            "OA", "NEA", "ALT", "A2_OTHER", "NONREF_ALLELE", "NEFFECT_ALLELE", "NEFFECTALLELE",
            "NONEFFECT_ALLELE", "OTHER_ALL", "OTHERALLELE", "NONEFFECTALLELE", "ALLELE0",
            "ALLELE_0", "ALT_ALLELE", "A_0", "NONCODED_ALLELE", "HM_OTHER_ALLELE"
        ],
        "BETA": [
            "BETA", "B", "EFFECT_BETA", "EFFECT", "EFFECTS", "SIGNED_SUMSTAT", "EST",
            "GWAS_BETA", "EFFECT_A1", "EFFECTA1", "EFFECT_NW", "HM_BETA", "LOG_ODDS"
        ],
        "OR": [
            "OR", "ODDS-RATIO", "ODDS_RATIO", "ODDSRATIO", "OR(MINALLELE)", "OR.LOGISTIC",
            "OR_RAN", "OR(A1)", "HM_ODDS_RATIO"
        ],
        "SE": [
            "SE", "STDER", "STDERR", "STD", "STANDARD_ERROR", "OR_SE", "STANDARDERROR",
            "STDERR_NW", "META.SE", "SE_DGC", "SE.2GC"
        ],
        "Z": [
            "Z", "ZSCORE", "Z-SCORE", "ZSTAT", "ZSTATISTIC", "GC_ZSCORE", "BETAZSCALE"
        ],
        "INFO": [
            "INFO", "IMPINFO", "IMPQUALITY", "INFO.PLINK", "INFO_UKBB", "INFO_UKB"
        ],
        "P": [
            "P", "PVALUE", "PVAL", "P_VALUE", "GC_PVALUE", "WALD_P", "P.VAL", "GWAS_P",
            "P-VALUE", "P-VAL", "FREQUENTIST_ADD_PVALUE", "P.VALUE", "P_VAL", "SCAN-P",
            "P.LMM", "META.PVAL", "P_RAN", "P.ADD", "P_BOLT_LMM"
        ],
        "N": [
            "N", "WEIGHT", "NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N",
            "N_COMPLETE_SAMPLES", "N_TOTAL", "N_SAMPLES", "N_ANALYZED", "NSAMPLES",
            "SAMPLESIZE", "SAMPLE_SIZE", "TOTAL_SAMPLE_SIZE", "TOTALSAMPLESIZE", "N_TOTAL_SUM"
        ],
        "N_CAS": [
            "N_CAS", "NCASE", "N_CASE", "N_CASES", "NCAS", "NCA", "NCASES", "CASES", "CASES_N"
        ],
        "N_CON": [
            "N_CON", "NCONTROL", "N_CONTROL", "N_CONTROLS", "NCON", "NCO", "N_CON",
            "NCONTROLS", "CONTROLS", "CONTROLS_N"
        ],
        "FREQ": [
            "FRQ", "FREQ", "MAF", "AF", "CEUAF", "FREQ1", "EAF", "FREQ1.HAPMAP",
            "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU", "EFFECT_ALLELE_FREQ", "FREQ.A1",
            "MA_FREQ", "MAF_NW", "FREQ_A1", "A1FREQ", "CODED_ALLELE_FREQUENCY",
            "FREQ_TESTED_ALLELE_IN_HRS", "EAF_HRC", "EAF_UKB", "EFFECT_ALLELE_FREQUENCY",
            "HM_EFFECT_ALLELE_FREQUENCY"
        ],
        "FRQ_A": ["FRQ_A", "FREQ_A", "F_A"],
        "FRQ_U": ["FRQ_U", "FREQ_U", "F_U"],
        "CHR": [
            "CHR", "CH", "CHROMOSOME", "CHROM", "CHR_BUILD38", "CHR_BUILD37", "CHR_BUILD36",
            "CHR_B38", "CHR_B37", "CHR_B36", "CHR_ID", "SCAFFOLD", "HG19CHR", "CHR.HG19",
            "CHR_HG19", "HG18CHR", "CHR.HG18", "CHR_HG18", "CHR_BP_HG19B37", "HG19CHRC",
            "HM_CHROM"
        ],
        "BP": [
            "ORIGBP", "BP", "POS", "POSITION", "LOCATION", "PHYSPOS", "GENPOS",
            "CHR_POSITION", "POS_B38", "POS_BUILD38", "POS_B37", "POS_BUILD37", "BP_HG19B37",
            "POS_B36", "POS_BUILD36", "POS.HG19", "POS.HG18", "POS_HG19", "POS_HG18",
            "BP_HG19", "BP_HG18", "BP.GRCH38", "BP.GRCH37", "POSITION(HG19)", "POSITION(HG18)",
            "POS(B38)", "POS(B37)", "BASE_PAIR_LOCATION", "HM_POS"
        ]
    }
    required_columns = {
        "identifier": ["SNP", "RSID", "CHR", "BP"],
        "effect_size": ["BETA", "OR", "LOG_ODDS", "SIGNED_SUMSTAT", "Z"],
        "p_or_se": ["P", "PVALUE", "PVAL", "P_VALUE", "SE", "STDER", "STDERR", "STD"],
        "effect_allele": ["A1"],
        "noneffect_allele": ["A2"],
        "allele_frequency": ["FREQ"]
    }

    def __init__(self, log_file="pipeline/logs/header_check.log"):
        self.reverse_dict = self.create_reverse_dict(self.ss_head_dict)
        self.log_file = log_file

    @staticmethod
    def create_reverse_dict(ss_head_dict):
        reverse_dict = {}
        for standard_name, aliases in ss_head_dict.items():
            for alias in aliases:
                reverse_dict[alias] = standard_name
        return reverse_dict

    def map_columns(self, column_names):
        """Maps column names to standard names."""
        return [self.reverse_dict.get(col.upper(), col).upper() for col in column_names]

    def check_required_columns(self, data, file_name):
        """
        Check if the data has the required columns.
        :param data: DataFrame containing the data
        :param file_name: name of the file being checked
        :return: boolean indicating success or failure
        """
        df = data
        df.columns = df.columns.str.upper()
        mapped_columns = self.map_columns(df.columns.tolist())
        column_set = set(mapped_columns)

        # Convert required column names to upper case for case-insensitive check
        required_identifier_columns = [col.upper() for col in self.required_columns["identifier"]]
        required_effect_size_columns = [col.upper() for col in self.required_columns["effect_size"]]
        required_p_or_se_columns = [col.upper() for col in self.required_columns["p_or_se"]]
        required_a1_columns = [col.upper() for col in self.required_columns["effect_allele"]]
        required_a2_columns = [col.upper() for col in self.required_columns["noneffect_allele"]]
        required_freq_columns = [col.upper() for col in self.required_columns["allele_frequency"]]

        # Check results
        has_identifier = any(col in column_set for col in required_identifier_columns)
        has_effect_size = any(col in column_set for col in required_effect_size_columns)
        has_p_or_se = any(col in column_set for col in required_p_or_se_columns)
        has_a1 = any(col in column_set for col in required_a1_columns)
        has_a2 = any(col in column_set for col in required_a2_columns)
        has_freq = any(col in column_set for col in required_freq_columns)

        # Log results
        with open(self.log_file, 'a') as log:
            log.write(f"File: {file_name} - ")
            if has_identifier and has_effect_size and has_p_or_se and has_a1 and has_a2:
                log.write("Status: PASS - All required columns are present.\n\n")
                return True
            else:
                log.write("Status: FAIL - Missing required columns.\n")
                if not has_identifier:
                    log.write("Missing: Identifier column (e.g., SNP, RSID, CHR, BP)\n")
                if not has_effect_size:
                    log.write("Missing: Effect size column (e.g., BETA, OR, LOG_ODDS, SIGNED_SUMSTAT, Z)\n")
                if not has_p_or_se:
                    log.write("Missing: P-value or SE column (e.g., P, SE)\n")
                if not has_a1:
                    log.write("Missing: Effect allele\n")
                if not has_a2:
                    log.write("Missing: Non effect allele\n")
                if not has_freq:
                    log.write("Missing: Allele frequency\n")
                log.write("\n")
                return False

    def read_file(self, file_path):
        """
        Read the file into a DataFrame.
        :param file_path: path to the file
        :return: DataFrame containing the data
        """
        if file_path.endswith(".gz"):
            with gzip.open(file_path, 'rt') as f:
                delimiter = '\t' if '.tsv' in file_path else ',' if '.csv' in file_path else ' '
                if delimiter == ' ':
                    return pd.read_csv(f, sep=r'\s+', nrows=1)  # Read only the header
                else:
                    return pd.read_csv(f, delimiter=delimiter, nrows=1)  # Read only the header
        else:
            delimiter = '\t' if '.tsv' in file_path else ',' if '.csv' in file_path else ' '
            if delimiter == ' ':
                return pd.read_csv(file_path, sep=r'\s+', nrows=1)  # Read only the header
            else:
                return pd.read_csv(file_path, delimiter=delimiter, nrows=1)  # Read only the header

    def validate_files(self, files_df):
        """
        Validate multiple GWAS files.
        :param files_df: DataFrame containing file paths
        :return: List of failed files
        """
        if os.path.exists(self.log_file):
            os.remove(self.log_file)  # Clear previous log file

        failed_files = []

        for file in files_df['path']:
            print(file)
            df_file = self.read_file(file)
            if not self.check_required_columns(df_file, file):
                failed_files.append(file)

        if failed_files:
            error_message = f"The following files did not pass validation: {', '.join(failed_files)}"
            print(Fore.RED + error_message + Style.RESET_ALL)
            sys.exit(1)

if __name__ == "__main__":
    validator = GWASValidator()
    data = "input_data/gwas_metadata.csv"
    df = pd.read_csv(data)
    validator.validate_files(df)
