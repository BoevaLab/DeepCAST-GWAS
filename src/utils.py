import numpy as np
import pandas as pd

from config import LEAD_SNP_COL, NLPVAL_COL, SNP_COL
from logger import setup_logger

logger = setup_logger(__name__)

def sad_columns_from_tracks(tracklist):
    return [f"SAD{track}" for track in tracklist]

# ISSUE with this function: Some track lists have a underscore, most don't
def create_filename_tracklist(phen_id: int):
    return f'tracks_phen{phen_id}.csv' # f'tracks_phen{phen_id}.csv'

def compute_adjusted_pval_threshold(len_initial, len_updated, pval_threshold):
    """
    Computes an adjusted p-value threshold to account for changes in the number of statistical tests.

    This is useful when multiple testing correction (i.e. Bonferroni) needs to be updated
    after filtering or reducing the number of SNPs or hypotheses being tested.

    The adjusted threshold is calculated as:
        adjusted_threshold = (len_initial / len_updated) * pval_threshold
        since:
            pval_threshold = pval_threshold / len_initial
            adjusted_threshold = pval_threshold / len_updated
            adjusted_threshold = (pval_threshold / len_updated) * (len_initial / len_initial)
            adjusted_threshold = (pval_threshold * len_initial) / (len_updated * len_initial)
            adjusted_threshold = (pval_threshold / len_initial) * (len_initial / len_updated)

    Parameters:
        len_initial (int): The original number of statistical tests (e.g., all SNPs).
        len_updated (int): The reduced number of tests after filtering.
        pval_threshold (float): The original p-value threshold (e.g., 5e-8).

    Returns:
        float: The adjusted p-value threshold, scaled to the new number of tests.
               Returns pval_threshold if len_updated is 0 to avoid division by zero.

    Note:
        This assumes a linear adjustment (e.g., Bonferroni-style correction).
    """
    if len_updated:
        return (len_initial / len_updated) * pval_threshold
    else:
        return pval_threshold
    
def pval_from_neglog10(sumstats):
    return 10 ** (-sumstats[NLPVAL_COL])

def neglog10_from_pval(pval):
    return -np.log10(pval)

def assign_to_baseline_loci(deepcast_sumstats: pd.DataFrame, baseline_sumstats: pd.DataFrame):
    """
    Assigns SNPs from the deepcast summary statistics to baseline lead SNP loci if applicable.

    This function performs a left join of the `deepcast_sumstats` DataFrame with the
    `baseline_sumstats` DataFrame on the 'snp' column, appending the corresponding 'lead_snp'
    from the baseline set where available.

    Parameters
    ----------
    deepcast_sumstats : pd.DataFrame
        Summary statistics DataFrame containing at least a 'snp' column.
    baseline_sumstats : pd.DataFrame
        Baseline clumped summary statistics containing 'snp' and 'lead_snp' columns.

    Returns
    -------
    pd.DataFrame
        A copy of `deepcast_sumstats` with an additional column 'lead_snp_baseline',
        indicating the lead SNP assignment from the baseline loci (or NaN if not matched).
    """
    return deepcast_sumstats.merge(baseline_sumstats[[SNP_COL, LEAD_SNP_COL]], how='left', on=SNP_COL, suffixes=[None, '_baseline'])