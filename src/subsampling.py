from datetime import date
import pandas as pd
from pathlib import Path
import argparse
import re
import json

from extract_metadata import generate_metadata_dict
from filtering import filter_significant, identify_coding_snps, identify_sad_relevant_snps
from input_output import read_in_fastgwa_sumstats, read_in_plink_sumstats, read_in_tracklist_from_file
from ld_based_clumping import get_significant_loci
from logger import setup_logger
from merge_sumstats_enformer_tracks import merge_sumstats_enformer_tracks
from config import RESULTS_DIR, SAD_DEVIATION_FACTOR as SD_RANGE, SNP_COL, SUBSAMPLE_DIR
from config import PVAL_THRESHOLD
from utils import compute_adjusted_pval_threshold, assign_to_baseline_loci

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

def main(sumstats_filename, tracklist_mapping_file = tracklist_MAPPING_FILE):
    logger.info(f'Starting Method on {date.today()}')

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
    merged_sumstats = merge_sumstats_enformer_tracks(sumstats, tracklist) #  pd.read_csv(Path(DEEPCAST_DIR) / "merged_sumstats_5517.csv") # Remove later ########################

    logger.info('Identifying deepcast-relevant SNPs:')
    merged_sumstats['in_coding_region'] = identify_coding_snps(merged_sumstats)
    coding_snps_deepcast: int = merged_sumstats['in_coding_region'].sum()
    logger.info(f'Found {coding_snps_deepcast} coding SNPs.')

    merged_sumstats['sad_relevant'] = identify_sad_relevant_snps(merged_sumstats, tracklist)
    sad_relevant_snps_deepcast: int = merged_sumstats['sad_relevant'].sum()
    logger.info(f'Found {sad_relevant_snps_deepcast} SNPs exceeding {SD_RANGE} standard deviations around the mean in the distribution of SAD values.')

    deepcast_snps = merged_sumstats[merged_sumstats['sad_relevant'] | merged_sumstats['in_coding_region']]
    logger.info(f"{len(deepcast_snps)} of {len(merged_sumstats)} snps are deepcast relevant")
    adjusted_pval_threshold = compute_adjusted_pval_threshold(len(merged_sumstats), len(deepcast_snps), PVAL_THRESHOLD)
    logger.info(f"Adjusted p-value threshold: {adjusted_pval_threshold}")

    logger.info('Get significant loci for baseline snps (retain tagging snps at adjusted threshold):')
    baseline_snps_clumped, lost_snps_baseline = get_significant_loci(merged_sumstats, pval_threshold_lead=PVAL_THRESHOLD, pval_threshold_clump=adjusted_pval_threshold)
    logger.info('Get significant loci for deepcast snps:')
    deepcast_snps_clumped, lost_snps_deepcast = get_significant_loci(deepcast_snps, pval_threshold_lead=adjusted_pval_threshold, pval_threshold_clump=adjusted_pval_threshold)

    logger.info('Identify deepcast lead snps contained in baseline loci:')
    deepcast_snps_clumped_with_overlap = assign_to_baseline_loci(deepcast_snps_clumped, baseline_snps_clumped)

    logger.info('Discard baseline tagging snps at adjusted threshold:')
    # TODO changed the code to use mask instead of returning filtered dataframe, check if it still works
    # baseline_snps_clumped = filter_significant(baseline_snps_clumped)
    baseline_snps_clumped = baseline_snps_clumped[filter_significant(baseline_snps_clumped)]

    lost_snps = lost_snps_baseline.merge(lost_snps_deepcast, how='outer', on=None) if not lost_snps_baseline.empty else lost_snps_baseline


    metadata: pd.DataFrame = generate_metadata_dict(sumstats_baseline=baseline_snps_clumped, sumstats_deepcast=deepcast_snps_clumped_with_overlap, n_cases = -1, adjusted_p_value_threshold=adjusted_pval_threshold)

    return baseline_snps_clumped, deepcast_snps_clumped_with_overlap, lost_snps, metadata, phenotype, n, replicate_id

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--run_id', type=str, default='unspecified', help='An identifier for this run')
    parser.add_argument('--filename', type=str, default=None, help='Sumstats filename (will try to retrieve from phen_id if not specified)')

    args = parser.parse_args()

    run_id: str = args.run_id
    filename = args.filename

    significant_snps_baseline, significant_snps_deepcast, snps_missing_from_reference, metadata, phenotype, n, replicate_id = main(filename)

    # Save results in deepcast_phenotypes folder under a subfolder for the run and the phenotype.
    phenotype_name = f'phen_{phenotype}'
    run_name = f'run_{run_id}'
    output_folder = Path(RESULTS_DIR) / run_name / phenotype_name / f'n_{n}' / f'replicate_{replicate_id}'
    output_folder.mkdir(parents=True, exist_ok=True)
    
    significant_snps_baseline.to_csv(output_folder / "significant_snps_baseline.csv", index=False)
    significant_snps_deepcast.to_csv(output_folder / "significant_snps_deepcast.csv", index=False)
    snps_missing_from_reference.to_csv(output_folder / "snps_missing_from_reference.csv", index=False)
    metadata.to_csv(output_folder / "metadata.csv", index=False)
    
    logger.info(f"Results saved in folder: {output_folder}.")
    