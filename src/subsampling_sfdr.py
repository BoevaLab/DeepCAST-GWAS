from datetime import date
import pandas as pd
from pathlib import Path
import argparse
import re
import json

from extract_metadata import generate_metadata_dict
from filtering import filter_significant, get_mean_sad_z_score, identify_coding_snps, identify_sad_relevant_snps
from input_output import read_in_fastgwa_sumstats, read_in_plink_sumstats, read_in_tracklist_from_file
from ld_based_clumping import get_significant_loci
from logger import setup_logger
from merge_sumstats_enformer_tracks import merge_sumstats_enformer_tracks
from config import ADJPVAL_COL, BIN_COL, FDR_THRESHOLD, NLPVAL_COL, ALT_COL, BP_COL, CHR_COL, REF_COL, RESULTS_DIR, SAD_DEVIATION_FACTOR as SD_RANGE, SNP_COL, STD_BIN_LABELS, STD_BINS, SUBSAMPLE_DIR
from config import PVAL_THRESHOLD
from fdr_adjust_pvals import perform_sfdr
from utils import compute_adjusted_pval_threshold, assign_to_baseline_loci
from extract_metadata import generate_metadata_dict


tracklist_MAPPING_FILE = SUBSAMPLE_DIR / "_tracklist_mapping.json"

PAT = re.compile(
    r'^(?P<pheno>.+?)_'                    # phenotype
    r'(?:fastgwa(?:_full)?|plink)_'        # method
    r'(?P<n>\d+)_'                         # cohort size
    r'(?P<replicate>\d+)\.'                # replicate ID before the dot
)

def parse_filename(filename: str) -> tuple[str, int, int]:
    m = PAT.search(filename)
    if not m:
        raise ValueError(f"Unrecognized filename format: {filename}")
    return m.group('pheno'), int(m.group('n')), int(m.group('replicate'))

logger = setup_logger(__name__)

def sfdr_main(sumstats_filename, tracklist_mapping_file = tracklist_MAPPING_FILE):
    logger.info(f'Starting Deepcast with sFDR on {date.today()}')

    logger.info(f'Sumstats file: {sumstats_filename}.')
    phenotype, n, replicate_id = parse_filename(sumstats_filename)

    if sumstats_filename.endswith(".fastGWA"):
        sumstats = read_in_fastgwa_sumstats(sumstats_filename)
    elif "_plink_" in sumstats_filename:
        sumstats = read_in_plink_sumstats(sumstats_filename)
    else:
        raise ValueError(f'Unsupported sumstats file format: {sumstats_filename}')

    with open(tracklist_mapping_file, 'r') as f:
        tracklist_mapping = json.load(f)
    if phenotype not in tracklist_mapping:
        raise ValueError(f'Phenotype {phenotype} not found in tracklist mapping file {tracklist_mapping_file}')

    tracklist_file = tracklist_mapping[phenotype]

    tracklist: list[int] = read_in_tracklist_from_file(tracklist_file)
    logger.info(f'Tracklist: {tracklist}')

    logger.info('Merging Sumstats with tracklists.')
    merged_sumstats = merge_sumstats_enformer_tracks(sumstats, tracklist)

    logger.info('Identifying deepcast-relevant SNPs:')
    merged_sumstats['in_coding_region'] = identify_coding_snps(merged_sumstats)
    coding_snps_deepcast: int = merged_sumstats['in_coding_region'].sum()
    logger.info(f'Found {coding_snps_deepcast} coding SNPs.')

    logger.info('Making SNP buckets based on SAD scores:')
    merged_sumstats['mean_sad_z_score'] = get_mean_sad_z_score(merged_sumstats, tracklist)
    merged_sumstats[BIN_COL] = pd.cut(merged_sumstats[~merged_sumstats['in_coding_region']]['mean_sad_z_score'], bins=STD_BINS, labels=STD_BIN_LABELS)

    logger.info(f"{len(STD_BINS)}  bins with counts:")
    merged_sumstats[BIN_COL] = merged_sumstats[BIN_COL].cat.add_categories(['coding'])
    merged_sumstats[BIN_COL] = merged_sumstats[BIN_COL].fillna('coding')

    logger.info(f"{len(STD_BINS)} bins with counts:") # since it's interval borders, it counts one too muchj, but we also need to add 1 for 'coding'
    logger.info(f"{merged_sumstats[BIN_COL].value_counts().sort_index()}")
    logger.info(f"Perform sFDR per stratum at threshold {FDR_THRESHOLD}:")
    fdr_adjusted_sumstats = perform_sfdr(merged_sumstats)
    # logger.info(f"Result has columns {fdr_adjusted_sumstats.columns}:")
    logger.info(f"Number of SNPs passing at threshold {FDR_THRESHOLD} {len(fdr_adjusted_sumstats[fdr_adjusted_sumstats[ADJPVAL_COL]<FDR_THRESHOLD])}:")

    logger.info('Get significant loci for baseline snps (retain fdr-significant as tagging snps):')
    baseline_clumping_snps = fdr_adjusted_sumstats[(fdr_adjusted_sumstats[ADJPVAL_COL]<FDR_THRESHOLD)|filter_significant(fdr_adjusted_sumstats)].copy()
    logger.debug(baseline_clumping_snps['chr'].value_counts(dropna=False))
    # TODO deal with missing chromosomes?
    baseline_snps_clumped, lost_snps_baseline = get_significant_loci(baseline_clumping_snps, pval_threshold_lead=PVAL_THRESHOLD, pval_threshold_clump=1)
    logger.info('Get significant loci under adjusted p-values:')
    fdr_snps_clumped, lost_fdr_snps = get_significant_loci(fdr_adjusted_sumstats, pval_threshold_lead=FDR_THRESHOLD, pval_threshold_clump=FDR_THRESHOLD, pval_col=ADJPVAL_COL, log_transformed=False)

    logger.info('Identify deepcast lead snps contained in baseline loci:')
    fdr_snps_with_baseline_loci = assign_to_baseline_loci(fdr_snps_clumped, baseline_snps_clumped)

    logger.info('Discard baseline tagging snps at adjusted threshold:')
    # TODO changed the code to use mask instead of returning filtered dataframe, check if it still works
    # baseline_snps_clumped = filter_significant(baseline_snps_clumped)
    baseline_snps_clumped = baseline_snps_clumped[filter_significant(baseline_snps_clumped)]

    lost_snps = lost_snps_baseline.merge(lost_fdr_snps, how='outer', on=None) if not lost_snps_baseline.empty else lost_snps_baseline

    metadata: pd.DataFrame = generate_metadata_dict(sumstats_baseline=baseline_snps_clumped, sumstats_deepcast=fdr_snps_with_baseline_loci, n_cases = -1)

    return baseline_snps_clumped, fdr_snps_with_baseline_loci, lost_snps, metadata, phenotype, n, replicate_id

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--run_id', type=str, default='unspecified', help='An identifier for this run')
    parser.add_argument('--filename', type=str, default=None, help='Sumstats filename (will try to retrieve from phen_id if not specified)')

    args = parser.parse_args()

    run_id: str = args.run_id
    filename = args.filename

    significant_snps_baseline, fdr_significant_snps, snps_missing_from_reference, metadata , phenotype, n, replicate_id = sfdr_main(filename)

    # Save results in deepcast_phenotypes folder under a subfolder for the run and the phenotype.
    phenotype_name = f'phen_{phenotype}'
    run_name = f'run_{run_id}'
    output_folder = Path(RESULTS_DIR) / run_name/ 'subsampling_sfdr' / phenotype_name / f'n_{n}' / f'replicate_{replicate_id}'
    output_folder.mkdir(parents=True, exist_ok=True)
    
    significant_snps_baseline.to_csv(output_folder / "significant_snps_baseline.csv", index=False)
    fdr_significant_snps.to_csv(output_folder / "fdr_significant_snps.csv", index=False)
    snps_missing_from_reference.to_csv(output_folder / "snps_missing_from_reference.csv", index=False)
    metadata.to_csv(output_folder / "metadata.csv", index=False)
    
    logger.info(f"Results saved in folder: {output_folder}.")
    