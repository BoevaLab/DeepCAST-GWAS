from datetime import date
import pandas as pd
from pathlib import Path
import argparse

from extract_metadata import generate_metadata_dict
from filtering import filter_significant, identify_coding_snps, identify_sad_relevant_snps
from input_output import read_in_sumstats, read_in_tracklist, retrieve_n_cases, retrieve_sumstats_filename
from ld_based_clumping import get_significant_loci
from logger import setup_logger
from merge_sumstats_enformer_tracks import merge_sumstats_enformer_tracks
from config import NLPVAL_COL, ALT_COL, BP_COL, CHR_COL, REF_COL, RESULTS_DIR, SAD_DEVIATION_FACTOR as SD_RANGE, SNP_COL
from config import PVAL_THRESHOLD
from utils import compute_adjusted_pval_threshold, assign_to_baseline_loci

logger = setup_logger(__name__)

def main(phen_id, filename:str|None=None, column_names={}):
    logger.info(f'Starting Method on {date.today()}')

    logger.info(f'Retrieving sumstats file name for Phenotype {phen_id}')
    sumstats_file: str = filename if filename is not None else retrieve_sumstats_filename(phen_id)
    n_cases = retrieve_n_cases(phen_id)

    logger.info(f'Sumstats file: {sumstats_file}.')
    sumstats: pd.DataFrame = read_in_sumstats(sumstats_file, column_names)

    tracklist: list[int] = read_in_tracklist(phen_id)
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


    metadata: pd.DataFrame = generate_metadata_dict(sumstats_baseline=baseline_snps_clumped, sumstats_deepcast=deepcast_snps_clumped_with_overlap, adjusted_p_value_threshold=adjusted_pval_threshold, n_cases = n_cases)

    return baseline_snps_clumped, deepcast_snps_clumped_with_overlap, lost_snps, metadata

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--phen', type=int, required=True, help='Phenotype id (index in phenotype manuscript)')
    parser.add_argument('--run_id', type=str, default='unspecified', help='An identifier for this run')
    parser.add_argument('--filename', type=str, default=None, help='Sumstats filename (will try to retrieve from phen_id if not specified)')
    parser.add_argument('--chr', type=str, help='Name of the chromosome column in the sumstats files')
    parser.add_argument('--bp', type=str, help='Name of the base pair position column in the sumstats files')
    parser.add_argument('--ref', type=str, help='Name of the reference allele column in the sumstats files')
    parser.add_argument('--alt', type=str, help='Name of the alternative allele column in the sumstats files')
    parser.add_argument('--neglogp', type=str, help='Name of the negative logarithmic p-value (log10) column in the sumstats files')

    args = parser.parse_args()

    phen_id: int = args.phen
    run_id: str = args.run_id
    filename = args.filename

    column_names = {}
    if args.chr: column_names[args.chr] = CHR_COL
    if args.bp: column_names[args.bp] = BP_COL
    if args.ref: column_names[args.ref] = REF_COL
    if args.alt: column_names[args.alt] = ALT_COL
    if args.neglogp: column_names[args.neglogp] = NLPVAL_COL

    significant_snps_baseline, significant_snps_deepcast, snps_missing_from_reference, metadata = main(phen_id, filename, column_names)

    # Save results in deepcast_phenotypes folder under a subfolder for the run and the phenotype.
    phenotype_name = f'phen{phen_id}'
    run_name = f'run_{run_id}'
    output_folder = Path(RESULTS_DIR) / run_name / phenotype_name
    output_folder.mkdir(parents=True, exist_ok=True)
    
    significant_snps_baseline.to_csv(output_folder / "significant_snps_baseline.csv", index=False)
    significant_snps_deepcast.to_csv(output_folder / "significant_snps_deepcast.csv", index=False)
    snps_missing_from_reference.to_csv(output_folder / "snps_missing_from_reference.csv", index=False)
    metadata.to_csv(output_folder / "metadata.csv", index=False)
    
    logger.info(f"Results saved in folder: {output_folder}.")
    