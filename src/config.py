from pathlib import Path
import numpy as np

## Method Parameters:
SAD_DEVIATION_FACTOR = 1
PVAL_THRESHOLD = 5e-8
FDR_THRESHOLD = 0.0025

R2_THRESHOLD = 0.2
KB_RADIUS = 500
WINDOW_SIZE = 500000

CHROMOSOMES = list(range(1, 23))

STD_BINS = [-np.inf, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, np.inf]
STD_BIN_LABELS = ["<-2σ", "[-2σ, -1σ]", "[-1σ, -0.5σ]", "[-0.5σ, -0.25σ]", "[-0.25σ, 0σ]", "[0σ, 0.25σ]", "[0.25σ, 0.5σ]", "[0.5σ, 1σ]", "[1σ, 2σ]", ">2σ"]


# Folder paths: (Adjust them to match your folder structure!)
# TODO delete this for submission?
DEEPCAST_DIR = Path('/cluster/work/boeva/lrabuzin/deepcast_gwas')

LOGS_DIR = DEEPCAST_DIR / 'logs'
RESULTS_DIR = DEEPCAST_DIR / 'results'
DATA_DIR = DEEPCAST_DIR / 'data'

SUMSTATS_DIR = DATA_DIR / 'ukbb_phens/icd10_sumstats_'
TRACKLISTS_DIR = DATA_DIR / 'tracklists'
REF_DIR_1KG = DATA_DIR / '1000genomes_as_csv'

# File paths:
LD_REFERENCE_PATH = DATA_DIR / '1kg_ld_reference/EUR'
EXON_REGIONS_PATH = DATA_DIR / 'genome_assembly/coding_snps.csv'
REF_DICT_PATH = DATA_DIR / 'reference_files_by_chr.json'
# TODO delete this
PLINK_PATH = '/cluster/apps/biomed/boeva/lrabuzin/conda/envs/deepcast_gwas/bin/plink'

# TODO do I need to change those?
ICD_PHENOTYPES_PATH = DATA_DIR / 'icd10_indices_n100.json'
PHEN_FILENAMES_PATH = DATA_DIR / 'phen_filenames.pkl'
PHEN_MANIFEST_PATH = DATA_DIR / 'phenotype_manifest.csv'

# Column names: (These can generally stay like this)
CHR_COL = 'chr'
BP_COL ='pos'
REF_COL ='ref'
ALT_COL ='alt'
NLPVAL_COL ='neglog10_pval_EUR'
PVAL_COL = 'pval'
ADJPVAL_COL = 'adjusted_pvals'
SNP_COL = 'snp'
LEAD_SNP_COL = 'lead_snp'
LEAD_SNP_BL_COL = LEAD_SNP_COL + '_baseline'
BIN_COL = 'sad_z_score_bin'