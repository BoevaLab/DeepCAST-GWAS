"""
DeepCAST-GWAS configuration.

The pipeline is configured via this module. For portability, defaults are relative
to the repository root, and you can override any path via environment variables.

Common overrides:
  - DEEPCAST_GWAS_DIR: repository root (default: inferred from this file location)
  - DEEPCAST_GWAS_DATA_DIR: data directory (default: ${DEEPCAST_GWAS_DIR}/data)
  - DEEPCAST_GWAS_RESULTS_DIR: results directory (default: ${DEEPCAST_GWAS_DIR}/results)
  - DEEPCAST_GWAS_LOGS_DIR: logs directory (default: ${DEEPCAST_GWAS_DIR}/logs)
  - DEEPCAST_GWAS_SUMSTATS_DIR: GWAS sumstats directory (default: ${DATA_DIR}/sumstats)
  - DEEPCAST_GWAS_TRACKLISTS_DIR: tracklists directory (default: ${DATA_DIR}/tracklists)
  - DEEPCAST_GWAS_REF_1KG_SAD_DIR: 1KG+SAD reference directory (default: ${DATA_DIR}/1000genomes_as_csv)
  - DEEPCAST_GWAS_REF_DICT_PATH: chromosome->reference filename mapping (default: ${DATA_DIR}/reference_files_by_chr.json)
  - DEEPCAST_GWAS_LD_REFERENCE_PREFIX: PLINK --bfile prefix for LD reference (default: ${DATA_DIR}/1kg_ld_reference/EUR)
  - DEEPCAST_GWAS_PLINK: path to plink binary (default: "plink" on PATH)
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

## -----------------------------
## Method parameters
## -----------------------------
SAD_DEVIATION_FACTOR = 1
PVAL_THRESHOLD = 5e-8
FDR_THRESHOLD = 0.0025

R2_THRESHOLD = 0.2
KB_RADIUS = 500
WINDOW_SIZE = 500000

CHROMOSOMES = list(range(1, 23))

STD_BINS = [-np.inf, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, np.inf]
STD_BIN_LABELS = [
    "<-2sd",
    "[-2sd, -1sd]",
    "[-1sd, -0.5sd]",
    "[-0.5sd, -0.25sd]",
    "[-0.25sd, 0sd]",
    "[0sd, 0.25sd]",
    "[0.25sd, 0.5sd]",
    "[0.5sd, 1sd]",
    "[1sd, 2sd]",
    ">2sd",
]

## -----------------------------
## Repository / filesystem layout
## -----------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
DEEPCAST_DIR = Path(os.environ.get("DEEPCAST_GWAS_DIR", str(REPO_ROOT))).expanduser().resolve()

LOGS_DIR = Path(os.environ.get("DEEPCAST_GWAS_LOGS_DIR", str(DEEPCAST_DIR / "logs"))).expanduser().resolve()
RESULTS_DIR = Path(os.environ.get("DEEPCAST_GWAS_RESULTS_DIR", str(DEEPCAST_DIR / "results"))).expanduser().resolve()
DATA_DIR = Path(os.environ.get("DEEPCAST_GWAS_DATA_DIR", str(DEEPCAST_DIR / "data"))).expanduser().resolve()

# Primary inputs
SUMSTATS_DIR = Path(os.environ.get("DEEPCAST_GWAS_SUMSTATS_DIR", str(DATA_DIR / "sumstats"))).expanduser().resolve()
TRACKLISTS_DIR = Path(os.environ.get("DEEPCAST_GWAS_TRACKLISTS_DIR", str(DATA_DIR / "tracklists"))).expanduser().resolve()
REF_DIR_1KG = Path(os.environ.get("DEEPCAST_GWAS_REF_1KG_SAD_DIR", str(DATA_DIR / "1000genomes_as_csv"))).expanduser().resolve()
REF_DICT_PATH = Path(os.environ.get("DEEPCAST_GWAS_REF_DICT_PATH", str(DATA_DIR / "reference_files_by_chr.json"))).expanduser().resolve()

# Variant annotation / preprocessing inputs
EXON_REGIONS_PATH = Path(
    os.environ.get("DEEPCAST_GWAS_CODING_SNPS_PATH", str(DATA_DIR / "genome_assembly" / "coding_snps.csv"))
).expanduser().resolve()
GFF_FILE_PATH = Path(
    os.environ.get("DEEPCAST_GWAS_GFF_PATH", str(DATA_DIR / "genome_assembly" / "genes.gff3"))
).expanduser().resolve()

# LD clumping inputs (PLINK)
LD_REFERENCE_PATH = Path(
    os.environ.get("DEEPCAST_GWAS_LD_REFERENCE_PREFIX", str(DATA_DIR / "1kg_ld_reference" / "EUR"))
).expanduser().resolve()
PLINK_PATH = os.environ.get("DEEPCAST_GWAS_PLINK", "plink")

# Pan-UKBB helpers (optional; used when resolving filenames / N cases from a phenotype id)
PHEN_FILENAMES_PATH = Path(os.environ.get("DEEPCAST_GWAS_PHEN_FILENAMES_PATH", str(DATA_DIR / "phen_filenames.pkl"))).expanduser().resolve()
PHEN_MANIFEST_PATH = Path(os.environ.get("DEEPCAST_GWAS_PHEN_MANIFEST_PATH", str(DATA_DIR / "phenotype_manifest.csv"))).expanduser().resolve()
ICD_PHENOTYPES_PATH = Path(os.environ.get("DEEPCAST_GWAS_ICD_PHENOTYPES_PATH", str(DATA_DIR / "icd10_indices_n100.json"))).expanduser().resolve()

# Optional directories used by auxiliary scripts (not required for basic / sFDR pipeline runs)
SUBSAMPLE_DIR = Path(os.environ.get("DEEPCAST_GWAS_SUBSAMPLE_DIR", str(DATA_DIR / "subsamples"))).expanduser().resolve()
FULL_SUMSTATS_DIR = Path(os.environ.get("DEEPCAST_GWAS_FULL_SUMSTATS_DIR", str(DATA_DIR / "full_sumstats"))).expanduser().resolve()
FINDOR_RESULTS_DIR = Path(os.environ.get("DEEPCAST_GWAS_FINDOR_RESULTS_DIR", str(DATA_DIR / "findor_results"))).expanduser().resolve()
REFORMATTED_SUMSTATS_DIR = Path(
    os.environ.get("DEEPCAST_GWAS_REFORMATTED_SUMSTATS_DIR", str(RESULTS_DIR / "ldsc_sumstats"))
).expanduser().resolve()

## -----------------------------
## Column names
## -----------------------------
CHR_COL = "chr"
BP_COL = "pos"
REF_COL = "ref"
ALT_COL = "alt"
NLPVAL_COL = "neglog10_pval_EUR"

PVAL_COL = "pval"
ADJPVAL_COL = "adjusted_pvals"

SNP_COL = "snp"
LEAD_SNP_COL = "lead_snp"
LEAD_SNP_BL_COL = LEAD_SNP_COL + "_baseline"
BIN_COL = "sad_z_score_bin"

# Optional annotation columns used by data prep utilities
CADD_COL = "phred"

# Misc flags referenced in some helper modules / notebooks
LD_BASED = True