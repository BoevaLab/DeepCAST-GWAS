import pandas as pd

from data_preparation.coding_snps.exon_regions import generate_coding_regions
from data_preparation.coding_snps.in_exon_region import find_coding_snps
from config import EXON_REGIONS_PATH, SNP_COL
from logger import setup_logger

logger = setup_logger(__name__)

def generate_coding_snps_list():
    coding_regions = generate_coding_regions().drop_duplicates(subset=['chr', 'start', 'end'])
    # Process files in parallel
    # merged_sumstats = merge_sumstats_with_reference('phecode-772.1-both_sexes.tsv.bgz')
    coding_snps: pd.DataFrame = find_coding_snps(coding_regions).drop_duplicates()
    # coding_snps = merged_sumstats[coding_snp_mask][['snp']].dropna().drop_duplicates(subset=['snp'])
    return coding_snps

if __name__ == "__main__":
    coding_snps = generate_coding_snps_list()
    pd.DataFrame({SNP_COL:coding_snps}).to_csv(EXON_REGIONS_PATH)