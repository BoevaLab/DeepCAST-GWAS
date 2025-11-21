from functools import partial
import pandas as pd
# from tqdm import tqdm
import concurrent.futures

from config import ALT_COL, BP_COL, CHR_COL, CHROMOSOMES, LEAD_SNP_COL, NLPVAL_COL, PVAL_COL, REF_COL, SNP_COL
from logger import setup_logger
from clump_snps_chr import clump_snps_chr

logger = setup_logger(__name__)

def call_chr_clumping(chr: int, pval_threshold_lead, pval_threshold_clump, summary_stats_by_chr, pval_col, log_transformed):
    return clump_snps_chr(
        pval_threshold_lead=pval_threshold_lead,
        pval_threshold_clump=pval_threshold_clump,
        sumstats=summary_stats_by_chr[chr],
        chr=chr,
        pval_col=pval_col,
        log_transformed=log_transformed
        )

def get_significant_loci(sumstats: pd.DataFrame, pval_threshold_lead, pval_threshold_clump, pval_col=NLPVAL_COL, log_transformed=True):

    # Pre-group summary stats by chromosome for faster joins
    summary_stats_by_chr = {
        chr_val: sub_df.copy()
        for chr_val, sub_df in sumstats.groupby(CHR_COL)
    }

    chromosomes = list(summary_stats_by_chr.keys())

    # Process files in parallel
    dataframes = []
    lost_snp_dataframes = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        try:
            logger.info(f'Identify lead SNPs and clump into loci by chromosome:')
            # TODO looks like this is not parallelized correctly?
            func = partial(call_chr_clumping,
                           summary_stats_by_chr=summary_stats_by_chr,
                           pval_threshold_lead=pval_threshold_lead,
                           pval_threshold_clump=pval_threshold_clump,
                           pval_col=pval_col,
                           log_transformed=log_transformed
                           )
            # results = list(tqdm(executor.map(func, CHROMOSOMES))) tqdm progress bar does not seem to work (at least for console output)
            # results = list(executor.map(func, CHROMOSOMES)) TODO be sure this works
            results = list(executor.map(func, chromosomes))
        except Exception as e:
            raise Exception(f'Exception performing ld clumping: {e}')
    
    # Combine results
    logger.info(f'Collect non-empty chromosome dataframes.')
    for res in results:
        if res[0] is not None and not res[0].empty:
            dataframes.append(res[0])
            lost_snp_dataframes.append(res[1])
            
    if not dataframes:
        logger.debug("LD clumping came back empty!")
        clumped_result = pd.DataFrame(columns=[CHR_COL,BP_COL,REF_COL,ALT_COL,NLPVAL_COL,SNP_COL,'in_coding_region','sad_relevant',PVAL_COL,LEAD_SNP_COL])
    else:
        logger.info('Concatenate chromosome dataframes.')
        clumped_result = pd.concat(dataframes, ignore_index=True)

    if not lost_snp_dataframes:
        lost_snps = pd.DataFrame(columns=[SNP_COL,CHR_COL,BP_COL,REF_COL,ALT_COL,NLPVAL_COL,'in_coding_region','sad_relevant',PVAL_COL])
    else:
        lost_snps = pd.concat(lost_snp_dataframes, ignore_index=True)
    logger.info(f"Merged data contains {len(clumped_result)} SNPs.")
    return clumped_result, lost_snps
    
