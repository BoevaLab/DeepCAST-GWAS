import numpy as np
from scipy.stats import norm
from scipy.stats import chi2
import argparse

from utils import pval_from_neglog10
from data_preparation.ldsc_files.merge_sumstats import merge_sumstats_reference
from input_output import read_in_sumstats, retrieve_n_cases, retrieve_sumstats_filename
from logger import setup_logger
from config import NLPVAL_COL, ALT_COL, BP_COL, CHR_COL, REF_COL, REFORMATTED_SUMSTATS_DIR

logger = setup_logger(__name__)

def reformat_sumstats(phen_id, filename:str|None=None, column_names={}):
    logger.info(f'Read in sumstats for phenotype {phen_id}.')
    sumstats_filename=retrieve_sumstats_filename(phen_id) if not filename else filename
    sumstats = read_in_sumstats(filename=sumstats_filename, column_names=column_names, read_in_beta=True)

    logger.info(f'Merge sumstats with reference genome to get SNP id and filter irrelevant.')
    merged_sumstats = merge_sumstats_reference(sumstats)

    logger.info(f'Retrieve number of cases for GWAS.')
    merged_sumstats['N'] = retrieve_n_cases(phen_id)

    snps = merged_sumstats.dropna(subset=['snp']).reset_index(drop=True)

    logger.info(f'{len(snps)} SNPs remaining after merge.')

    # hopefully fine but I'm not sure:
    logger.info(f'Compute Z score from p-value.')
    snps['P'] = pval_from_neglog10(snps)
    # snps["Z"] = np.sign(snps['beta_EUR']) * norm.isf(snps["P"] / 2)
    snps["my_Z"] = np.sign(snps['beta_EUR']) * norm.isf(snps["P"] / 2)
    snps['Z'] = np.sqrt(chi2.isf(snps["P"], 1))

    logger.info('Rename columns and drop irrelevant to match ldsc format.')
    sumstats_formatted = snps.rename(columns={'snp':'SNP', 'ref':'A1', 'alt':'A2'})[['Z', 'SNP', 'N', 'P', 'A1', 'A2']]

    return sumstats_formatted, sumstats_filename

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

    ldsc_sumstats, filename = reformat_sumstats(phen_id, filename, column_names)

    # Save results in deepcast_phenotypes folder under a subfolder for the run and the phenotype.
    phenotype_name = f'phen{phen_id}'
    run_name = f'run_{run_id}'
    output_folder = REFORMATTED_SUMSTATS_DIR / run_name
    output_folder.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Saving sumstats.")
    ldsc_sumstats.to_csv(output_folder/filename, sep="\t", index=False, compression="gzip")
    
    logger.info(f"Results saved in folder: {output_folder} as {filename}.")
    