import pandas as pd
import numpy as np

from config import REF_DIR_1KG
from logger import setup_logger

logger = setup_logger(__name__)

def coding_snp_mask_chr(chr:int, reference_file_path, coding_regions_chr:pd.DataFrame):
    # sumstats_chr.sort_values('pos')
    # coding_regions_chr.sort_values
    # idx_right = 0
    # for idx_left, row in sumst
    # 1) build an IntervalIndex from the range-table

    try:
        chunk = pd.read_csv(reference_file_path, usecols=['chr', 'pos', 'snp'])
    except Exception as e:
        logger.error(f"Error reading {reference_file_path}: {e}")
        return None

    # logger.debug(f'For chromosome {chr_val}: Creating index.')
    iv = pd.IntervalIndex.from_arrays(coding_regions_chr['start'],
                                    coding_regions_chr['end'],
                                    closed='both')   # [start, end] inclusive
    # logger.debug(f'For chromosome {chr_val}: Created interval index.')

    # 2) look up every number at once
    _, missing = iv.get_indexer_non_unique(chunk['pos']) # type: ignore
    # -1 → not contained
    # logger.debug(f"Found {len(missing)} SNPs for chromosome {chr_val} not contained in a coding region from {len(sumstats_chr)} SNPs")

    mask = np.ones(len(chunk), dtype=bool)

    mask[missing] = False

    logger.debug(f"Mask for chromosome {chr} contains {mask.sum()} coding SNPs out of {len(chunk)} SNPs")

    return chunk[mask]['snp']