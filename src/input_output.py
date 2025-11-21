import gzip
import pandas as pd
from pathlib import Path
import re
import pickle

from config import DATA_DIR, FINDOR_RESULTS_DIR, LEAD_SNP_COL, SNP_COL, SUMSTATS_DIR, TRACKLISTS_DIR, PHEN_FILENAMES_PATH, NLPVAL_COL, ALT_COL, BP_COL, CHR_COL, REF_COL, SUBSAMPLE_DIR, FULL_SUMSTATS_DIR
from logger import setup_logger
from utils import create_filename_tracklist

import numpy as np

logger = setup_logger(__name__)

# TODO rewrite this for efficient lookups (e.g with a shelve db)

def retrieve_sumstats_filename(phen_id):
    with open(PHEN_FILENAMES_PATH, 'rb') as f:
        my_dict = pickle.load(f)
    return str(my_dict.get(phen_id))

def retrieve_n_cases(phen_id):
    with open(DATA_DIR / 'phen_cases.pkl', 'rb') as f:
        my_dict = pickle.load(f)
    return my_dict.get(phen_id)

def retrieve_n_controls(phen_id):
    with open(DATA_DIR / 'phen_controls.pkl', 'rb') as f:
        my_dict = pickle.load(f)
    return my_dict.get(phen_id)

def read_in_sumstats(filename: str, column_names: dict, drop_na_pvals:bool=True, read_in_beta:bool=False):
    """
    Reads summary statistics from a file and returns a DataFrame with selected columns.
    Expects to find the columns 'chr', 'pos', 'ref', 'alt', 'neglog10_pval_EUR'.

    Parameters:
        path (str): Path to the summary statistics file. Supports both plain text and .tsv.bgz files.

    Returns:
        pd.DataFrame: DataFrame containing the selected columns from the summary statistics file.
    """
    
    dir_sumstats = Path(SUMSTATS_DIR)
    sumstats_path = dir_sumstats / filename

    default_columns = [CHR_COL, BP_COL, REF_COL, ALT_COL, NLPVAL_COL]

    # Rename according to column_names:
    reversed_column_names = {val:key for key, val in column_names.items()}

    selected_columns = []

    for default_name in default_columns:
        custom_name = reversed_column_names.get(default_name)
        if custom_name:
            selected_columns.append(custom_name)
        else:
            selected_columns.append(default_name)

    logger.info('Reading in sumstats')
    logger.debug(f'sumstats path {sumstats_path}')

    if(read_in_beta):
        selected_columns.append('beta_EUR')

    file = gzip.open(sumstats_path, "rt") if str(sumstats_path).endswith(".tsv.bgz") else open(sumstats_path, "r")
    try:
        sumstats = pd.concat(pd.read_csv(file, sep="\t", usecols=selected_columns, chunksize=100000), ignore_index=True)
    except ValueError as e:
        raise RuntimeError("Reading in the Summary Statistics failed. "
                        "By default, the columns 'chr', 'pos', 'ref', 'alt' and 'neglog10_pval_EUR'"
                        "are expected. If the input file columns deviate, specify them in the"
                        "flags --chr --bp --ref --alt --neglogp.") from e

    # logger.debug(f'Read in sumstats with column names: \n {sumstats.head()}')

    if column_names:
        # rename back to always have the same column names:
        sumstats.rename(columns=column_names, inplace=True)
        logger.debug(f'Changed sumstats column names to default: \n {sumstats.head()}')

    # TODO: Can I drop this? We'll see after running this a couple of times? Or should this always be here anyway?
    if sumstats[NLPVAL_COL].isna().any():
        logger.debug(f'Found {sumstats[NLPVAL_COL].isna().sum()} of {len(sumstats[NLPVAL_COL])} N/A rows in column {NLPVAL_COL}.')
        if drop_na_pvals:
            logger.info('Dropping N/A rows and resetting the index')
            sumstats.dropna(subset=[NLPVAL_COL], inplace=True)
            sumstats.reset_index(drop=True, inplace=True)
    else:
        logger.debug(f'Found no N/A rows in column {NLPVAL_COL}.')

    return sumstats

def read_in_findor_results(filename):
    
    filepath = FINDOR_RESULTS_DIR / filename

    columns = ['SNP', 'N', 'Z', 'A1', 'A2', 'P', 'P_weighted']

    logger.info(f'Reading in reweighted FINDOR sumstats for {filename}')
    logger.debug(f'sumstats path {filepath}')
    
    file = pd.concat(pd.read_csv(filepath, usecols=columns, sep='\s+', chunksize=100000), ignore_index=True) # type: ignore

    logger.debug(f'Read in sumstats with column names: \n {file.head()}')

    logger.debug(f'Found {file["P"].isna().sum()} of {len(file["P"])} N/A rows in column P.')
    logger.debug(f'Found {file["P_weighted"].isna().sum()} of {len(file["P_weighted"])} N/A rows in column P_weighted.')
    file.dropna(subset=["P", "P_weighted"], inplace=True)
    file.reset_index(drop=True, inplace=True)

    return file

def read_in_fastgwa_sumstats(filename: str, drop_na_pvals:bool=True):
    filepath = SUBSAMPLE_DIR / filename
    column_names = {'chr': 'CHR', 'pos': 'POS', 'ref': 'A2', 'alt': 'A1', 'neglog10_pval_EUR': 'P'}
    columns = list(column_names.values())

    logger.info(f'Reading in fastGWA sumstats for {filename}')
    logger.debug(f'sumstats path {filepath}')

    file = pd.read_csv(filepath, usecols=columns, sep='\t')

    sumstats = file.rename(columns={v: k for k, v in column_names.items()})

    # Actually converting the p-values to negative log10 p-values:
    sumstats[NLPVAL_COL] = -np.log10(sumstats[NLPVAL_COL].astype(float))

    if sumstats[NLPVAL_COL].isna().any():
        logger.debug(f'Found {sumstats[NLPVAL_COL].isna().sum()} of {len(sumstats[NLPVAL_COL])} N/A rows in column {NLPVAL_COL}.')
        if drop_na_pvals:
            logger.info('Dropping N/A rows and resetting the index')
            sumstats.dropna(subset=[NLPVAL_COL], inplace=True)
            sumstats.reset_index(drop=True, inplace=True)
    else:
        logger.debug(f'Found no N/A rows in column {NLPVAL_COL}.')

    return sumstats

def read_in_plink_sumstats(filename: str, drop_na_pvals:bool=True):

    filepath = SUBSAMPLE_DIR / filename
    column_names = {'chr': '#CHROM', 'pos': 'POS', 'ref': 'REF', 'alt': 'ALT', 'neglog10_pval_EUR': 'P'}
    columns = list(column_names.values())

    logger.info(f'Reading in plink sumstats for {filename}')
    logger.debug(f'sumstats path {filepath}')

    file = pd.read_csv(filepath, usecols=columns, sep='\t')
    
    sumstats = file.rename(columns={v: k for k, v in column_names.items()})

    sumstats[NLPVAL_COL] = -np.log10(sumstats[NLPVAL_COL].astype(float))

    if sumstats[NLPVAL_COL].isna().any():
        logger.debug(f'Found {sumstats[NLPVAL_COL].isna().sum()} of {len(sumstats[NLPVAL_COL])} N/A rows in column {NLPVAL_COL}.')
        if drop_na_pvals:
            logger.info('Dropping N/A rows and resetting the index')
            sumstats.dropna(subset=[NLPVAL_COL], inplace=True)
            sumstats.reset_index(drop=True, inplace=True)
    else:
        logger.debug(f'Found no N/A rows in column {NLPVAL_COL}.')

    return sumstats

def read_in_full_sumstats(filename: str, drop_na_pvals:bool=True):
    filepath = FULL_SUMSTATS_DIR / filename
    column_names = {'chr': 'CHR', 'pos': 'POS', 'ref': 'A2', 'alt': 'A1', 'neglog10_pval_EUR': 'P'}
    columns = list(column_names.values())

    logger.info(f'Reading in fastGWA sumstats for {filename}')
    logger.debug(f'sumstats path {filepath}')

    file = pd.read_csv(filepath, usecols=columns, sep='\t')

    sumstats = file.rename(columns={v: k for k, v in column_names.items()})

    # Actually converting the p-values to negative log10 p-values:
    sumstats[NLPVAL_COL] = -np.log10(sumstats[NLPVAL_COL].astype(float))

    if sumstats[NLPVAL_COL].isna().any():
        logger.debug(f'Found {sumstats[NLPVAL_COL].isna().sum()} of {len(sumstats[NLPVAL_COL])} N/A rows in column {NLPVAL_COL}.')
        if drop_na_pvals:
            logger.info('Dropping N/A rows and resetting the index')
            sumstats.dropna(subset=[NLPVAL_COL], inplace=True)
            sumstats.reset_index(drop=True, inplace=True)
    else:
        logger.debug(f'Found no N/A rows in column {NLPVAL_COL}.')

    return sumstats

def read_in_tracklist(phen_id: int):
    logger.info('Reading in tracklist.')
    dir_tracklists = Path(TRACKLISTS_DIR)
    filename_tracklist = create_filename_tracklist(phen_id)
    tracklist_path = dir_tracklists / filename_tracklist
    logger.debug(f'Tracklist path {tracklist_path}.')

    tracklist = pd.read_csv(tracklist_path, header=None)
    return tracklist.iloc[:, 0].tolist()

def read_in_tracklist_from_file(filename: str):
    logger.info('Reading in tracklist.')
    dir_tracklists = Path(TRACKLISTS_DIR)
    filename_tracklist = filename
    tracklist_path = dir_tracklists / filename_tracklist
    logger.debug(f'Tracklist path {tracklist_path}.')

    tracklist = pd.read_csv(tracklist_path, header=None)
    return tracklist.iloc[:, 0].tolist()

def read_missing_snps_from_plink_log(log_path):
    log_lines = Path(log_path).read_text().splitlines()
    missing_snps = []
    pattern = re.compile(r"Warning: '(\S+)' is missing from the main dataset, and is a top variant\.")
    
    for line in log_lines:
        match = pattern.search(line)
        if match:
            missing_snps.append(match.group(1))
    return missing_snps

def process_clumped(clumped_path):
    """
    Parses a PLINK .clumped file and maps each SNP (lead or clumped/tagging) 
    to its corresponding lead SNP.

    Parameters:
        clumped_path : str or Path
            Path to the PLINK .clumped file generated with `--clump`.
        snp_col : str, optional
            Name of the output column representing individual SNPs (default: 'snp').
        lead_snp_col : str, optional
            Name of the output column representing the lead SNP each SNP is assigned to (default: 'lead_snp').

    Returns:
        pandas.DataFrame
            A DataFrame with two columns:
            - SNP_COL: all SNPs that were retained or assigned to a clump
            - LEAD_SNP_COL: the lead SNP that each SNP belongs to

    Notes:
        - SNPs listed in the 'SP2' column are considered clumped/tagging SNPs and are assigned
        to the corresponding lead SNP from the same row.
        - SNPs marked as 'NONE' in 'SP2' are treated as lead SNPs with no associated tags.
    """
        
    clumped_df = pd.read_csv(clumped_path, sep='\s+') # type: ignore
    snp_rows = []
    lead_snp_rows = []

    for _, row in clumped_df.iterrows():
        snp = row['SNP']
        snp_rows.append(snp)
        lead_snp_rows.append(snp)
        if row['SP2'] != 'NONE':
            for tag in row['SP2'].split(','):
                tagging_snp = tag.split('(')[0]
                snp_rows.append(tagging_snp)
                lead_snp_rows.append(snp)
    
    return pd.DataFrame({SNP_COL: snp_rows, LEAD_SNP_COL: lead_snp_rows})