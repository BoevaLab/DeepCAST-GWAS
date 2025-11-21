import pandas as pd

from config import KB_RADIUS, LD_BASED, PVAL_THRESHOLD, R2_THRESHOLD, LEAD_SNP_COL, LEAD_SNP_BL_COL
from logger import setup_logger

logger = setup_logger(__name__)

def generate_metadata_dict(sumstats_baseline, sumstats_deepcast, n_cases, adjusted_p_value_threshold=PVAL_THRESHOLD):
    # num loci baseline
    num_loci_baseline = sumstats_baseline[LEAD_SNP_COL].nunique()
    # num loci deepcast
    num_loci_deepcast = sumstats_deepcast[LEAD_SNP_COL].nunique()
    # overlapping loci:
    num_recovered_loci = sumstats_deepcast[LEAD_SNP_BL_COL].nunique()
    num_overlapping_loci = 0
    for (_, df) in sumstats_deepcast.groupby(LEAD_SNP_COL):
        if not df[LEAD_SNP_BL_COL].isna().all():
            num_overlapping_loci += 1
    num_newly_discovered_loci =  num_loci_deepcast - num_overlapping_loci
    # changed this because it was wrong:
    # num_overlapping_loci = len(set(sumstats_deepcast[LEAD_SNP_BL_COL].dropna()))
    # num_newly_discovered_loci = len(set(sumstats_deepcast[sumstats_deepcast[LEAD_SNP_BL_COL].isna()][LEAD_SNP_COL]))
    # lead_snps in coding regions
    num_coding_lead_snps_deepcast = sumstats_deepcast['in_coding_region'].sum() # len(sumstats_deepcast[sumstats_deepcast['in_coding_region']])
    num_coding_lead_snps_baseline = sumstats_baseline['in_coding_region'].sum() # len(sumstats_baseline[sumstats_baseline['in_coding_region']])
    # loci recovered
    percentage_loci_recovered = num_recovered_loci/num_loci_baseline if num_loci_baseline != 0 else 1


    return pd.DataFrame({
        "num_loci_baseline": [num_loci_baseline],
        "num_loci_deepcast": [num_loci_deepcast],
        "num_overlapping_loci": [num_overlapping_loci],
        "num_recovered_loci": [num_recovered_loci],
        "percentage_loci_recovered": [percentage_loci_recovered],
        "num_newly_discovered_loci": [num_newly_discovered_loci],
        "num_coding_lead_snps_deepcast": [num_coding_lead_snps_deepcast],
        "num_coding_lead_snps_baseline": [num_coding_lead_snps_baseline],
        "adjusted_p_value_threshold": [adjusted_p_value_threshold],
        "n_cases": [n_cases]
    })