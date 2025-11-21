from pathlib import Path
import subprocess
import tempfile
import pandas as pd

from config import BP_COL, CHR_COL, PVAL_COL, SNP_COL, PLINK_PATH, LD_REFERENCE_PATH, KB_RADIUS, R2_THRESHOLD
from input_output import process_clumped, read_missing_snps_from_plink_log
from utils import pval_from_neglog10
from logger import setup_logger

logger = setup_logger(__name__)

def clump_snps_chr(sumstats: pd.DataFrame, pval_threshold_lead, pval_threshold_clump, chr, pval_col, log_transformed):

    # TODO drop na values reset index? Try it out, see if I get any na values
    # logger.debug(f"Sumstats dataframe has columns {sumstats.columns} for chromosome {chr}.")
    # Added if condition to check if PVAL column is already there (need it for findor weighted p-values)
    # if not sumstats[PVAL_COL].isna().all():
    # sumstats[PVAL_COL] = pval_from_neglog10(sumstats)

    if(log_transformed):
        sumstats[PVAL_COL] = pval_from_neglog10(sumstats)
    else:
        sumstats[PVAL_COL] = sumstats[pval_col]

    if sumstats[PVAL_COL].isna().any():
        logger.debug(f"Found {sumstats[PVAL_COL].isna().sum()} N/A values in new pval column for chromosome {chr}.")

    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp()

    # Prepare input file for PLINK
    # temp_dir = Path(PLINK_OUTPUT_DIR) / f"chr{chr}"
    # temp_dir.mkdir(exist_ok=True)
    assoc_file = Path(temp_dir) / "plink_input.assoc"
    clumped_output = Path(temp_dir) / "plink_output"
    
    # Format the association file for PLINK
    plink_input = sumstats[[SNP_COL, CHR_COL, BP_COL, PVAL_COL]].copy()
    plink_input.columns = ['SNP', 'CHR', 'BP', 'P']
    plink_input.to_csv(assoc_file, sep='\t', index=False)

    # Run PLINK clumping
    # logger.info(f"Running PLINK clumping for chromosome {chr}")
    cmd_str = f"{PLINK_PATH} --bfile {LD_REFERENCE_PATH} --clump {assoc_file} " \
                f"--clump-p1 {pval_threshold_lead} --clump-p2 {pval_threshold_clump} " \
                f"--clump-r2 {R2_THRESHOLD} --clump-kb {KB_RADIUS} " \
                f"--out {clumped_output}"
    
    try:
        subprocess.run(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as e:
        logger.error(f"Error executing PLINK: {e}")
        return pd.DataFrame(), []
    
    # Process PLINK results
    # logger.info(f"Processing PLINK results for chromosome {chr}")
    clumped_path = Path(temp_dir) / f"{clumped_output}.clumped"

    if not clumped_path.exists() or clumped_path.stat().st_size == 0:
        logger.info(f"PLINK clumping found no significant SNPs for chromosome {chr}.")
        return pd.DataFrame(), []
    
    # Get lead SNPs from clumped file
    lead_snps: pd.DataFrame = process_clumped(clumped_path)
    # if lead_snps.empty:
    #     logger.info(f"Warning: No lead SNPs identified for chromosome {chr}.")
        
    # logger.info(f"Found {len(lead_snps)} lead SNPs after clumping for chromosome {chr}")
    
    # Add the clumping information to the original data:
    sumstats_merged = sumstats.merge(lead_snps, how='inner', on=SNP_COL)

    log_path = Path(f"{clumped_output}.log")

    lost_snp_ids = pd.DataFrame(read_missing_snps_from_plink_log(log_path), columns=[SNP_COL])

    # logger.info(f"{len(lost_snp_ids)} significant SNPs were not found in reference data for chromosome {chr}.")
    logger.info(f"Found {len(lead_snps)} lead SNPs after clumping for chromosome {chr}, {len(lost_snp_ids)} significant SNPs were not found in reference data.")

    lost_snps = lost_snp_ids.merge(sumstats, how='left', on=SNP_COL)

    return sumstats_merged, lost_snps