from scipy.stats import norm
import numpy as np
import concurrent.futures
from functools import partial
import json
import pandas as pd
from pathlib import Path

from config import ALT_COL, BP_COL, CADD_COL, CHR_COL, DATA_DIR, NLPVAL_COL, REF_COL, SNP_COL, REF_DIR_1KG, REF_DICT_PATH, CHROMOSOMES
from logger import setup_logger

logger = setup_logger(__name__)

def merge_chr_reference(reference_filename: str, sumstats_chr: pd.DataFrame, chr: int, cadd_scores=False):
    """
    Process a single CSV file:
      - Reads the CSV and extracts the chromosome using the precompiled regex pattern.
      - Retrieves the corresponding summary stats from summary_stats_by_chr.
      - Sets the index on the merge keys and performs a join.
      - Returns the merged DataFrame for that CSV.
    """

    dir_ref = DATA_DIR / 'cadd_reference' if cadd_scores else DATA_DIR / '1kg_reference_genome'
    reference_file_path = dir_ref / f'cadd_reference_chr{chr}.csv' if cadd_scores else dir_ref / reference_filename

    cols = [CHR_COL, BP_COL, REF_COL, ALT_COL, SNP_COL, CADD_COL] if cadd_scores else [CHR_COL, BP_COL, REF_COL, ALT_COL, SNP_COL]

    # logger.info(f'Reading in {reference_file_path}:')
    try:
        chunk = pd.read_csv(reference_file_path, usecols=cols)
    except Exception as e:
        logger.error(f"Error reading {reference_file_path}: {e}")
        return None
    
    # Set index on chunk to match summary stats index
    chunk = chunk.set_index([CHR_COL, BP_COL, REF_COL, ALT_COL])
    # logger.info(f'Successfully read in reference. Merging sumstats with reference for Chromosome {chr}.')
    merged_result = sumstats_chr.join(chunk, how='left', rsuffix='_ref').reset_index()
    
    # Keep only rows with valid SAD and p_value values
    merged_result = merged_result[merged_result[NLPVAL_COL].notna()]
    # logger.info(f'Successfully merged files. Merged sumstats for Chromosome {chr} contains {len(merged_result)} columns after merge.')
    logger.info(f"Processed file {reference_filename}")
    return merged_result

def call_chr_reference_merge(chr: int, reference_filename_by_chr, summary_stats_by_chr, cadd_scores=False):
    return merge_chr_reference(
        reference_filename=reference_filename_by_chr[chr],
        sumstats_chr=summary_stats_by_chr[chr],
        chr=chr,
        cadd_scores=cadd_scores
        )

def merge_sumstats_reference(sumstats: pd.DataFrame, cadd_scores=False):
    if not all(col in sumstats.columns for col in [CHR_COL, BP_COL, REF_COL, ALT_COL]):
    # if not all(col in sumstats.columns for col in [CHR_COL, BP_COL, REF_COL, ALT_COL, 'beta_EUR']):
        logger.critical(f"Dataframe columns misspecified, should contain 'chr', 'pos', 'ref', 'alt', got {sumstats.columns} instead.")
        # raise ValueError(f"Dataframe columns misspecified, should contain {CHR_COL, BP_COL, REF_COL, ALT_COL, 'beta_EUR'}, got {sumstats.columns} instead.")

    # Pre-group summary stats by chromosome for faster joins
    summary_stats_by_chr = {
        chr_val: sub_df.copy().set_index([CHR_COL, BP_COL, REF_COL, ALT_COL])
        for chr_val, sub_df in sumstats.groupby(CHR_COL)
    }

    # open dictionary of merged reference file by chromosome:
    with open(REF_DICT_PATH, 'r') as f:
        reference_filename_by_chr = json.load(f)
        reference_filename_by_chr = {int(k): v for k, v in reference_filename_by_chr.items()}

    # Process files in parallel
    dataframes = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        try:
            logger.info(f'Merge summary stats with reference genome{" with CADD scores" if cadd_scores else ""} by chromosome:')
            func = partial(call_chr_reference_merge,
                      reference_filename_by_chr=reference_filename_by_chr,
                      summary_stats_by_chr=summary_stats_by_chr,
                      cadd_scores=True)
            # results = list(tqdm(executor.map(func, CHROMOSOMES))) tqdm progress bar does not seem to work (at least for console output)
            results = list(executor.map(func, CHROMOSOMES))
        except Exception as e:
            raise Exception(f'Exception merging sumstats for chromosome with reference file: {e}')
    
    # Combine results
    logger.info(f'Collect non-empty chromosome dataframes.')
    for res in results:
        if res is not None and not res.empty:
            dataframes.append(res)
            
    if not dataframes:
        logger.debug("Merge of reference with sumstats came back empty!")
        return pd.DataFrame()
        
    logger.info('Merge chromosome dataframes.')
    df_summary_stats_result = pd.concat(dataframes, ignore_index=True)
    logger.info(f"Merged data contains {len(df_summary_stats_result)} SNPs")
    return df_summary_stats_result




