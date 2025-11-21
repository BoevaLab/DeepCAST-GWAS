import json
from pathlib import Path

import pandas as pd
from config import CHR_COL, BP_COL, REF_COL, ALT_COL, DATA_DIR, SNP_COL
from config import REF_DIR_1KG, CHROMOSOMES
from config import ALT_COL, BP_COL, CHR_COL, DATA_DIR, REF_COL, REF_DICT_PATH, REF_DIR_1KG
from logger import setup_logger
import concurrent.futures
from functools import partial

logger = setup_logger(__name__)


def create_cadd_reference_chr(chr):
    out_folder=DATA_DIR / 'cadd_reference'
    cadd_dir = DATA_DIR / "cadd_scores"
    dir_1kg = Path(REF_DIR_1KG) / '..' / '1kg_reference_genome'

    reference_file_path = dir_1kg / f'1000G.MAF_threshold=0.005.{chr}_combined.csv'
    
    logger.info(f'Read in cadd scores for chromosome {chr} from {reference_file_path}.')
    col_names = [CHR_COL, BP_COL, REF_COL, ALT_COL, 'raw', 'phred']
    cadd_scores_chr = pd.read_csv(cadd_dir / f"chr{chr}_cadd.tsv", sep="\t", header=None, names=col_names)
    
    logger.info(f'Read in reference for chromosome {chr} from {reference_file_path}.')
    reference_chr = pd.read_csv(reference_file_path, usecols=[CHR_COL, BP_COL, REF_COL, ALT_COL, SNP_COL])

    logger.info(f'Merging reference and caddd scores for chromosome {chr}.')
    merged_reference = cadd_scores_chr.merge(reference_chr, how='left', on=[CHR_COL, BP_COL, REF_COL, ALT_COL])
    logger.info(f'{len(merged_reference)} SNPs after merge.')
    merged_reference = merged_reference.dropna().reset_index()
    logger.info(f'{len(merged_reference)} SNPs after dropping NA rsid rows.')

    logger.info(f'Saving reference.')
    out_name = f'cadd_reference_chr{chr}.csv'
    merged_reference.to_csv(out_folder / out_name)
    logger.info(f'Saved chromosome {chr} cadd score reference to {(out_name / out_folder).resolve()}')

def create_cadd_reference_files():

    with open(REF_DICT_PATH, 'r') as f:
        reference_filename_by_chr = json.load(f)
        reference_filename_by_chr = {int(k): v for k, v in reference_filename_by_chr.items()}

    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
        try:
            logger.info(f'Merge summary stats with enformer tracks by chromosome:')
            list(executor.map(create_cadd_reference_chr, CHROMOSOMES))
        except Exception as e:
            raise Exception(f'Exception merging sumstats for chromosome with reference file: {e}')

    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     try:
    #         logger.info(f'Merge summary stats with enformer tracks by chromosome:')
    #         results = list(executor.map(create_cadd_reference_chr, CHROMOSOMES))
    #     except Exception as e:
    #         raise Exception(f'Exception merging sumstats for chromosome with reference file: {e}')
    #     logger.info(f'Saving reference.')
    #     out_folder=DATA_DIR / 'cadd_reference'
    #     for (chr, df) in zip(CHROMOSOMES, results): 
    #         out_name = f'cadd_reference_chr{chr}.csv'
    #         df.to_csv(out_folder / out_name)
    #         logger.info(f'Saved chromosome {chr} cadd score reference to {(out_name / out_folder).resolve()}')


if __name__ == "__main__":
    create_cadd_reference_files()

    