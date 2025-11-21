from functools import partial
import json
from pathlib import Path
import pandas as pd
import numpy as np
import concurrent.futures

from config import BP_COL, CHR_COL, CHROMOSOMES, REF_DICT_PATH, REF_DIR_1KG, SNP_COL
from data_preparation.coding_snps.in_exon_region_chr import coding_snp_mask_chr
from logger import setup_logger

logger = setup_logger(__name__)

def call_coding_snp_mask_chr(chr, ref_dir_path: Path, reference_filenames_by_chr, coding_regions_by_chr):
    return coding_snp_mask_chr(chr, ref_dir_path / reference_filenames_by_chr[chr], coding_regions_by_chr[chr])

def find_coding_snps(coding_regions):
    # summary_stats_by_chr = {
    #     chr_val: sub_df.copy()
    #     for chr_val, sub_df in sumstats.groupby(CHR_COL)
    # }

    # open dictionary of merged reference file by chromosome:
    with open(REF_DICT_PATH, 'r') as f:
        reference_filename_by_chr = json.load(f)
        reference_filename_by_chr = {int(k): v for k, v in reference_filename_by_chr.items()}
        
    coding_regions_by_chr = {
        chr_val: sub_df.copy()
        for chr_val, sub_df in coding_regions.groupby('chr')
    }
    
    coding_snps: list[str] = []
    # Without limiting the workers, the kernel crashes
    with concurrent.futures.ProcessPoolExecutor(max_workers=7) as executor:
        try:
            func = partial(call_coding_snp_mask_chr,
                           ref_dir_path=REF_DIR_1KG,
                           reference_filenames_by_chr=reference_filename_by_chr,
                           coding_regions_by_chr=coding_regions_by_chr)
            results = list(executor.map(func, CHROMOSOMES))
        except Exception as e:
            raise Exception(f'Exception filtering coding snps with reference file: {e}')
    
    # Combine results
    for res in results:
        if res is not None:
            coding_snps.extend(res)
    
    logger.debug(f'Concatenated lists. Found {len(coding_snps)} coding snps')

    return pd.DataFrame({'snp':coding_snps})