# from tqdm import tqdm
import concurrent.futures
from functools import partial
import json
import pandas as pd

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

def merge_sumstats_enformer_tracks(sumstats: pd.DataFrame, tracklist: list[int]):
    # TODO: do this everywhere or nowhere?
    # Check if Dataframe has correct columns:
    if not all(col in sumstats.columns for col in [CHR_COL, BP_COL, REF_COL, ALT_COL]):
        # logger.critical(f"Dataframe columns misspecified, should contain 'chr', 'pos', 'ref', 'alt', got {sumstats.columns} instead.")
        raise ValueError(f"Dataframe columns misspecified, should contain {CHR_COL, BP_COL, REF_COL, ALT_COL}, got {sumstats.columns} instead.")

    logger.info('Merging sumstats with enformer SAD tracks')
    # Initialize SAD columns with track names
    sad_columns = sad_columns_from_tracks(tracklist)
    logger.debug(f'SAD columns: {sad_columns}')

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