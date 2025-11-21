import pandas as pd
from pathlib import Path

from config import ALT_COL, BP_COL, CHR_COL, NLPVAL_COL, REF_COL, SNP_COL, REF_DIR_1KG
from logger import setup_logger

logger = setup_logger(__name__)

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