# from tqdm import tqdm
import concurrent.futures
from functools import partial
import json
import pandas as pd


import pandas as pd
from pathlib import Path

from config import ALT_COL, BP_COL, CHR_COL, NLPVAL_COL, REF_COL, SNP_COL, REF_DIR_1KG

from config import ALT_COL, BP_COL, CHR_COL, REF_COL, REF_DICT_PATH, CHROMOSOMES
from logger import setup_logger
from merge_chr_sad_scores import merge_chr_sad_scores
from utils import sad_columns_from_tracks

logger = setup_logger(__name__)

def call_chr_merge(chr: int, reference_filename_by_chr, sad_columns, summary_stats_by_chr):
    return merge_chr_sad_scores(
        reference_filename=reference_filename_by_chr[chr],
        sad_columns=sad_columns,
        sumstats_chr=summary_stats_by_chr[chr]
        )

def generate_reference_genome():
    # # TODO: do this everywhere or nowhere?
    # # Check if Dataframe has correct columns:
    # if not all(col in sumstats.columns for col in [CHR_COL, BP_COL, REF_COL, ALT_COL]):
    #     # logger.critical(f"Dataframe columns misspecified, should contain 'chr', 'pos', 'ref', 'alt', got {sumstats.columns} instead.")
    #     raise ValueError(f"Dataframe columns misspecified, should contain {CHR_COL, BP_COL, REF_COL, ALT_COL}, got {sumstats.columns} instead.")

    # logger.info('Merging sumstats with enformer SAD tracks')
    # # Initialize SAD columns with track names
    # sad_columns = sad_columns_from_tracks(tracklist)
    # logger.debug(f'SAD columns: {sad_columns}')

    # Pre-group summary stats by chromosome for faster joins
    # summary_stats_by_chr = {
    #     chr_val: sub_df.copy().set_index([CHR_COL, BP_COL, REF_COL, ALT_COL])
    #     for chr_val, sub_df in sumstats.groupby(CHR_COL)
    # }

    # open dictionary of merged reference file by chromosome:
    with open(REF_DICT_PATH, 'r') as f:
        reference_filename_by_chr = json.load(f)
        reference_filename_by_chr = {int(k): v for k, v in reference_filename_by_chr.items()}

    # Process files in parallel
    dataframes = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        try:
            logger.info(f'Merge summary stats with enformer tracks by chromosome:')
            func = partial(call_chr_merge,
                      reference_filename_by_chr=reference_filename_by_chr,
                      sad_columns=sad_columns,
                      summary_stats_by_chr=summary_stats_by_chr)
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
        logger.debug("Merge of SAD data with sumstats came back empty!")
        return pd.DataFrame()
        
    logger.info('Merge chromosome dataframes.')
    df_summary_stats_result = pd.concat(dataframes, ignore_index=True)
    logger.info(f"Merged data contains {len(df_summary_stats_result)} SNPs")
    return df_summary_stats_result


def merge_chr_sad_scores(reference_filename: str, sad_columns: list[str], sumstats_chr: pd.DataFrame):
    """
    Process a single CSV file:
      - Reads the CSV and extracts the chromosome using the precompiled regex pattern.
      - Retrieves the corresponding summary stats from summary_stats_by_chr.
      - Sets the index on the merge keys and performs a join.
      - Returns the merged DataFrame for that CSV.
    """
    dir_1kg = Path(REF_DIR_1KG)
    reference_file_path = dir_1kg / reference_filename

    try:
        chunk = pd.read_csv(reference_file_path, usecols=[CHR_COL, BP_COL, REF_COL, ALT_COL, SNP_COL] + sad_columns)
    except Exception as e:
        logger.error(f"Error reading {reference_file_path}: {e}")
        return None
    
    # Set index on chunk to match summary stats index
    chunk = chunk.set_index([CHR_COL, BP_COL, REF_COL, ALT_COL])
    merged_result = sumstats_chr.join(chunk, how='left').reset_index()
    
    # Keep only rows with valid SAD and p_value values
    merged_result = merged_result[merged_result[sad_columns[0]].notna() & merged_result[NLPVAL_COL].notna()]
    logger.info(f"Processed file {reference_filename}")
    return merged_result