from datetime import date
import pandas as pd
from pathlib import Path
import argparse

from extract_metadata import generate_metadata_dict
from filtering import filter_significant, get_mean_sad_z_score, identify_coding_snps, identify_sad_relevant_snps
from input_output import read_in_sumstats, read_in_tracklist, retrieve_n_cases, retrieve_sumstats_filename
from ld_based_clumping import get_significant_loci
from logger import setup_logger
from merge_sumstats_enformer_tracks import merge_sumstats_enformer_tracks
from config import ADJPVAL_COL, BIN_COL, FDR_THRESHOLD, NLPVAL_COL, ALT_COL, BP_COL, CHR_COL, REF_COL, RESULTS_DIR, SAD_DEVIATION_FACTOR as SD_RANGE, SNP_COL, STD_BIN_LABELS, STD_BINS
from config import PVAL_THRESHOLD
from fdr_adjust_pvals import perform_sfdr
from utils import compute_adjusted_pval_threshold, assign_to_baseline_loci
from extract_metadata import generate_metadata_dict

logger = setup_logger(__name__)

def sfdr_main(phen_id, filename:str|None=None, column_names={}):
    logger.info(f'Starting Deepcast with sFDR on {date.today()}')

    logger.info(f'Retrieving sumstats file name for Phenotype {phen_id}')
    sumstats_file: str = filename if filename is not None else retrieve_sumstats_filename(phen_id)
    n_cases = retrieve_n_cases(phen_id)

    logger.info(f'Sumstats file: {sumstats_file}.')
    sumstats: pd.DataFrame = read_in_sumstats(sumstats_file, column_names)

    tracklist: list[int] = read_in_tracklist(phen_id)
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

    metadata: pd.DataFrame = generate_metadata_dict(sumstats_baseline=baseline_snps_clumped, sumstats_deepcast=fdr_snps_with_baseline_loci, n_cases = n_cases)

    return baseline_snps_clumped, fdr_snps_with_baseline_loci, lost_snps, metadata

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

    significant_snps_baseline, fdr_significant_snps, snps_missing_from_reference, metadata = sfdr_main(phen_id, filename, column_names)

    # Save results in deepcast_phenotypes folder under a subfolder for the run and the phenotype.
    phenotype_name = f'phen{phen_id}'
    run_name = f'run_{run_id}'
    output_folder = Path(RESULTS_DIR) / 'deepcast_sfdr' / run_name / phenotype_name
    output_folder.mkdir(parents=True, exist_ok=True)
    
    significant_snps_baseline.to_csv(output_folder / "significant_snps_baseline.csv", index=False)
    fdr_significant_snps.to_csv(output_folder / "fdr_significant_snps.csv", index=False)
    snps_missing_from_reference.to_csv(output_folder / "snps_missing_from_reference.csv", index=False)
    metadata.to_csv(output_folder / "metadata.csv", index=False)
    
    logger.info(f"Results saved in folder: {output_folder}.")
    