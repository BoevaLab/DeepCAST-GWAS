# from tqdm import tqdm
import concurrent.futures
from functools import partial
import json
import pandas as pd
import statsmodels.stats.multitest as mt

from config import ADJPVAL_COL, ALT_COL, BIN_COL, BP_COL, CHR_COL, FDR_THRESHOLD, PVAL_COL, REF_COL, REF_DICT_PATH, CHROMOSOMES
from logger import setup_logger
from merge_chr_sad_scores import merge_chr_sad_scores
from utils import pval_from_neglog10, sad_columns_from_tracks

logger = setup_logger(__name__)

def run_sfdr(df):
    df[PVAL_COL] = pval_from_neglog10(df)
    # TODO why is it sorted?
    _, adjusted_pvals, _, _ = mt.multipletests(pvals=df[PVAL_COL], method='fdr_bh', alpha=FDR_THRESHOLD)
    df[ADJPVAL_COL] = adjusted_pvals
    return df

def perform_sfdr(sumstats: pd.DataFrame, bin_col=BIN_COL):
    # TODO: do this everywhere or nowhere?
    # Check if Dataframe has correct columns:
    if not all(col in sumstats.columns for col in [bin_col]):
        # logger.critical(f"Dataframe columns misspecified, should contain 'chr', 'pos', 'ref', 'alt', got {sumstats.columns} instead.")
        raise ValueError(f"Dataframe columns misspecified, should contain {bin_col}, got {sumstats.columns} instead.")

    logger.info('Performing sFDR')

    # Pre-group summary stats by bin for parallelization
    strata_dict = {
        bin: sub_df
        for bin, sub_df in sumstats.groupby(bin_col)
    }

    strata=strata_dict.values()

    # Process files in parallel
    dataframes = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        try:
            logger.info(f'Perform FDR by stratum:')
            results = list(executor.map(run_sfdr, strata))
        except Exception as e:
            raise Exception(f'Exception merging sumstats for chromosome with reference file: {e}')
    
    # Combine results
    logger.info(f'Collect non-empty chromosome dataframes.')
    for res in results:
        if res is not None and not res.empty:
            dataframes.append(res)
            
    if not dataframes:
        logger.debug("Merge of SAD data with sumstats came back empty!")
        return pd.DataFrame()
        
    logger.info('Merge chromosome dataframes.')
    df_summary_stats_result = pd.concat(dataframes, ignore_index=True)
    logger.info(f"Merged data contains {len(df_summary_stats_result)} SNPs")
    return df_summary_stats_result