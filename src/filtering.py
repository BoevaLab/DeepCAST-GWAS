import numpy as np
import pandas as pd

from config import EXON_REGIONS_PATH, NLPVAL_COL, PVAL_THRESHOLD, SNP_COL
from config import SAD_DEVIATION_FACTOR as SD_RANGE
from utils import neglog10_from_pval, sad_columns_from_tracks
from logger import setup_logger

logger = setup_logger(__name__)

def identify_coding_snps(sumstats: pd.DataFrame, snp_col_reference = 'snp'):
    exon_regions: pd.DataFrame = pd.read_csv(EXON_REGIONS_PATH)
    coding_region_set = exon_regions[snp_col_reference].to_list()
    return sumstats[SNP_COL].isin(coding_region_set)

def identify_sad_relevant_snps(sumstats: pd.DataFrame, tracklist: list[int]):
    sad_columns = sad_columns_from_tracks(tracklist)
    track_data = sumstats[sad_columns].to_numpy()
    
    # normalize by column instead:
    track_data_norm = (track_data - track_data.mean(axis=0)) / track_data.std(axis=0)
    in_sad_range_mask = ((track_data_norm < - SD_RANGE) + (track_data > SD_RANGE)).any(axis=1)
    
    # Before: overall mean and std:
    # mean = np.mean(track_data)
    # sd = np.std(track_data)
    
    # in_sad_range_mask = ((track_data < mean - SD_RANGE * sd) + (track_data > mean + SD_RANGE * sd)).any(axis=1)
    return in_sad_range_mask

# def get_mean_sad_score(sumstats: pd.DataFrame, tracklist: list[int]):
#     sad_columns = sad_columns_from_tracks(tracklist)
#     track_data = sumstats[sad_columns].to_numpy()
#     logger.debug(f'track_data dimension: {track_data.shape}')
#     row_mean = np.mean(track_data, axis=1)
#     logger.debug(f'row_mean dimension: {row_mean.shape}')
#     logger.debug(f'row_mean length: {len(row_mean):.2e}')
    
#     return row_mean

# with column-wise normalization first:
def get_mean_sad_z_score(sumstats: pd.DataFrame, tracklist: list[int]):
    sad_columns = sad_columns_from_tracks(tracklist)
    track_data = sumstats[sad_columns].to_numpy()

    track_data_norm = (track_data - track_data.mean(axis=0)) / track_data.std(axis=0)
    row_mean_z = np.mean(track_data_norm, axis=1)
    
    return row_mean_z

def get_mean_absolute_sad_z_score(sumstats: pd.DataFrame, tracklist: list[int]):
    sad_columns = sad_columns_from_tracks(tracklist)
    track_data = sumstats[sad_columns].to_numpy()
    track_data_norm = (track_data - track_data.mean(axis=0)) / track_data.std(axis=0)

    # Average of absolute z-scores per row
    row_mean_abs_z = np.mean(np.abs(track_data_norm), axis=1)

    return row_mean_abs_z


def filter_significant(sumstats):
    """
    Identifies SNPs that exceed a significance threshold based on -log10(p-value).

    Parameters:
        sumstats (pd.DataFrame): A DataFrame containing summary statistics, including a column with -log10(p-values).

    Returns:
        pd.DataFrame: A filtered DataFrame containing only SNPs with -log10(p-value) >= -log10(pval_threshold).
    
    Note:
        The function assumes that the specified column contains -log10 transformed p-values.
    """
    neglog10_threshold = neglog10_from_pval(PVAL_THRESHOLD)
    mask = sumstats[NLPVAL_COL] >= neglog10_threshold
    return mask
    # changed code, see if it still works
    # return sumstats[sumstats[NLPVAL_COL] >= neglog10_threshold]