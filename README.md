# DeepCAST-GWAS (camera-ready)

DeepCAST-GWAS is a GWAS post-processing pipeline that integrates **Enformer-derived SAD scores** and coding-variant annotations to prioritize variants and loci.

This repository provides two “production” entrypoints:
- **Baseline + DeepCAST (FWER-style thresholding)**: `src/main.py`
- **DeepCAST with stratified FDR (sFDR)**: `src/deepcast_sfdr.py`

The pipeline covers the complete workflow:
- **Data preprocessing** (reference formatting; coding SNP list generation; optional LDSC formatting utilities)
- **Integration of Enformer-derived SAD scores** into GWAS summary statistics
- **Variant filtering / prioritization** (coding SNPs and SAD outliers)
- **LD-based locus definition** via PLINK clumping
- **Stratified FDR** (binning by SAD z-score + coding stratum)

Example SLURM scripts for HPC are provided below (copy/paste templates).

## Quickstart

1) **Configure paths** in `src/config.py` (portable defaults are relative to the repo). Most users only need to set:
- `SUMSTATS_DIR` (GWAS summary stats)
- `TRACKLISTS_DIR` (phenotype-specific Enformer tracklists)
- `REF_DIR_1KG` + `REF_DICT_PATH` (1KG reference augmented with Enformer SAD tracks)
- `LD_REFERENCE_PATH` + `PLINK_PATH` (for clumping)

You can override most paths via environment variables (see docstring at the top of `src/config.py`).

2) **Run baseline + DeepCAST (FWER-style)**

```bash
python3 -u src/main.py --phen 123 --run_id demo_run
```

3) **Run DeepCAST with stratified FDR (sFDR)**

```bash
python3 -u src/deepcast_sfdr.py --phen 123 --run_id demo_run
```

Outputs are written under `results/` (see “Outputs” below).

## Installation

DeepCAST-GWAS is pure Python, but requires a working PLINK installation for LD clumping.

### Python dependencies

At minimum you need:
- `python>=3.10`
- `pandas`, `numpy`, `scipy`, `statsmodels`

Optional (only for some data-prep utilities):
- `openai`, `pydantic`, `python-dotenv` (used by `src/data_preparation/curate_sad_tracks.py`)

### PLINK

Ensure `plink` is available on your `$PATH`, or set `DEEPCAST_GWAS_PLINK` to an explicit binary path.

## Data layout (inputs)

DeepCAST-GWAS expects a small set of well-defined inputs. The default layout is:

```text
deepcast_gwas/
  data/
    sumstats/                       # GWAS summary stats files
    tracklists/                     # phenotype-specific Enformer track indices
    1000genomes_as_csv/             # 1KG reference files augmented with SAD columns
    reference_files_by_chr.json     # chr -> reference filename mapping
    genome_assembly/
      coding_snps.csv               # list of coding SNP rsIDs (generated)
  results/
  logs/
  src/
```

### GWAS summary statistics

The pipeline reads a single GWAS summary-statistics file for a phenotype. The reader supports plain text and `.tsv.bgz` (bgzip) files.

Required columns (default names):
- `chr`: chromosome (int-like)
- `pos`: base-pair position (int-like)
- `ref`: reference allele
- `alt`: alternate allele
- `neglog10_pval_EUR`: \(-\log_{10}(p)\)

If your columns differ, you can map them via CLI flags:
- `--chr`, `--bp`, `--ref`, `--alt`, `--neglogp`

For non-indexed phenotypes, pass the file directly with `--filename` (relative to `SUMSTATS_DIR`).

### Tracklists (per phenotype)

Each phenotype needs a tracklist: a one-column CSV of Enformer track indices, one per line.

By default, the filename is:
- `data/tracklists/tracks_phen<PHEN_ID>.csv`

(This naming is controlled by `src/utils.py:create_filename_tracklist`.)

### 1KG reference augmented with Enformer SAD scores

To integrate Enformer-derived SAD scores, DeepCAST-GWAS joins your GWAS variants onto per-chromosome reference CSVs.

Inputs:
- `data/reference_files_by_chr.json`: mapping `"1" -> "<filename_for_chr1>.csv"`, …, `"22" -> ...`
- `data/1000genomes_as_csv/<that filename>`: reference CSVs

Each reference CSV must contain the merge keys and annotation columns:
- **Merge keys**: `chr`, `pos`, `ref`, `alt`
- **Required**: `snp` (rsID)
- **SAD columns**: `SAD<track_index>` for every track listed in the phenotype’s tracklist

You can download the reference genome merged with the Enformer tracks [here](www.provide_download_link.com).

### Coding SNP list (exon-derived)

DeepCAST-GWAS marks variants as “coding” by membership in `data/genome_assembly/coding_snps.csv` (rsIDs).

You can download exon regions as pdf [here](www.provide_download_link.com).

To (re)generate the coding SNP list:
- Set `GFF_FILE_PATH` in `src/config.py` to a gene annotation GFF/GFF3
- Run:

```bash
python3 -u src/data_preparation/coding_snps/generate_coding_snp_list.py
```

## Running the pipeline

### Baseline + DeepCAST (FWER-style thresholding)

Run:

```bash
python3 -u src/main.py --phen 123 --run_id my_run
```

High-level steps (as implemented in `src/main.py`):
- Read GWAS sumstats
- Load phenotype tracklist
- Merge sumstats with SAD tracks (per chromosome) using `reference_files_by_chr.json`
- Identify DeepCAST-relevant SNPs:
  - coding SNPs (from `coding_snps.csv`)
  - SAD outliers (based on `SAD_DEVIATION_FACTOR`)
- Compute an adjusted p-value threshold for the DeepCAST subset
- Define loci via LD clumping (PLINK) for:
  - baseline (genome-wide significant lead SNPs; keeps DeepCAST-adjusted tagging variants)
  - DeepCAST-relevant variants (lead SNPs at the adjusted threshold)

### DeepCAST with stratified FDR (sFDR)

Run:

```bash
python3 -u src/deepcast_sfdr.py --phen 123 --run_id my_run
```

High-level steps (as implemented in `src/deepcast_sfdr.py`):
- Read GWAS sumstats
- Merge SAD tracks
- Compute per-variant mean SAD z-score (across the phenotype’s tracks)
- Assign each SNP to a **stratum**:
  - non-coding SNPs: binned by SAD z-score (`STD_BINS` / `STD_BIN_LABELS`)
  - coding SNPs: forced into a dedicated `"coding"` stratum
- Perform BH-FDR within each stratum, then clump to loci with PLINK

## Outputs

### Baseline + DeepCAST (`src/main.py`)

Output directory:
- `results/run_<RUN_ID>/phen<PHEN_ID>/`

Files:
- `significant_snps_baseline.csv`
- `significant_snps_deepcast.csv`
- `snps_missing_from_reference.csv`
- `metadata.csv`

### sFDR (`src/deepcast_sfdr.py`)

Output directory:
- `results/deepcast_sfdr/run_<RUN_ID>/phen<PHEN_ID>/`

Files:
- `significant_snps_baseline.csv`
- `fdr_significant_snps.csv`
- `snps_missing_from_reference.csv`
- `metadata.csv`

## Running on HPC (SLURM)

The scripts below are written as **templates** (placeholders are intentional). You typically only need to edit:
- conda/module activation
- resource requests (`--cpus-per-task`, `--mem`, `--time`)
- the phenotype list file used for arrays
- environment overrides for `src/config.py` paths (optional but recommended on shared filesystems)

### Baseline + DeepCAST (single phenotype)

```bash
#!/bin/bash
#SBATCH --job-name=deepcast_baseline
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/dev/null

set -euo pipefail

DEEPCAST_GWAS_DIR="${DEEPCAST_GWAS_DIR:-/path/to/deepcast_gwas}"
PHEN_ID="${PHEN_ID:-123}"
RUN_ID="${RUN_ID:-${SLURM_JOB_ID}}"

source ~/.bashrc
# conda activate deepcast_gwas

cd "${DEEPCAST_GWAS_DIR}"
log_dir="logs/baseline/run_${RUN_ID}"
mkdir -p "${log_dir}"
exec > "${log_dir}/job_${SLURM_JOB_ID}.out" 2>&1

python3 -u src/main.py \
  --phen "${PHEN_ID}" \
  --run_id "${RUN_ID}"
```

### Baseline + DeepCAST (phenotype array)

```bash
#!/bin/bash
#SBATCH --job-name=deepcast_baseline_arr
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=1-100
#SBATCH --output=/dev/null

set -euo pipefail

DEEPCAST_GWAS_DIR="${DEEPCAST_GWAS_DIR:-/path/to/deepcast_gwas}"
PHEN_LIST="${PHEN_LIST:-/path/to/phen_ids.txt}"  # one integer phenotype id per line
RUN_ID="${RUN_ID:-${SLURM_ARRAY_JOB_ID}}"

source ~/.bashrc
# conda activate deepcast_gwas

cd "${DEEPCAST_GWAS_DIR}"
log_dir="logs/baseline/run_${RUN_ID}"
mkdir -p "${log_dir}"
exec > "${log_dir}/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" 2>&1

PHEN_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PHEN_LIST}")"
echo "PHEN_ID=${PHEN_ID} RUN_ID=${RUN_ID}"

python3 -u src/main.py \
  --phen "${PHEN_ID}" \
  --run_id "${RUN_ID}"
```

### sFDR (single phenotype)

```bash
#!/bin/bash
#SBATCH --job-name=deepcast_sfdr
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/dev/null

set -euo pipefail

DEEPCAST_GWAS_DIR="${DEEPCAST_GWAS_DIR:-/path/to/deepcast_gwas}"
PHEN_ID="${PHEN_ID:-123}"
RUN_ID="${RUN_ID:-${SLURM_JOB_ID}}"

source ~/.bashrc
# conda activate deepcast_gwas

cd "${DEEPCAST_GWAS_DIR}"
log_dir="logs/sfdr/run_${RUN_ID}"
mkdir -p "${log_dir}"
exec > "${log_dir}/job_${SLURM_JOB_ID}.out" 2>&1

python3 -u src/deepcast_sfdr.py \
  --phen "${PHEN_ID}" \
  --run_id "${RUN_ID}"
```

### sFDR (phenotype array)

```bash
#!/bin/bash
#SBATCH --job-name=deepcast_sfdr_arr
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=1-100
#SBATCH --output=/dev/null

set -euo pipefail

DEEPCAST_GWAS_DIR="${DEEPCAST_GWAS_DIR:-/path/to/deepcast_gwas}"
PHEN_LIST="${PHEN_LIST:-/path/to/phen_ids.txt}"  # one integer phenotype id per line
RUN_ID="${RUN_ID:-${SLURM_ARRAY_JOB_ID}}"

source ~/.bashrc
# conda activate deepcast_gwas

cd "${DEEPCAST_GWAS_DIR}"
log_dir="logs/sfdr/run_${RUN_ID}"
mkdir -p "${log_dir}"
exec > "${log_dir}/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" 2>&1

PHEN_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PHEN_LIST}")"
echo "PHEN_ID=${PHEN_ID} RUN_ID=${RUN_ID}"

python3 -u src/deepcast_sfdr.py \
  --phen "${PHEN_ID}" \
  --run_id "${RUN_ID}"
```

## Configuration reference

Method thresholds:
- **`SAD_DEVIATION_FACTOR`**: SAD outlier cutoff (in SD units) for DeepCAST relevance
- **`PVAL_THRESHOLD`**: baseline genome-wide significance threshold
- **`FDR_THRESHOLD`**: sFDR threshold for the adjusted p-values
- **`STD_BINS` / `STD_BIN_LABELS`**: bins used for SAD-based strata in sFDR

LD clumping:
- **`R2_THRESHOLD`**, **`KB_RADIUS`**: PLINK clumping parameters
- **`LD_REFERENCE_PATH`**: `plink --bfile` prefix for LD reference (e.g., 1KG EUR)
- **`PLINK_PATH`**: plink binary (or `"plink"`)

## Troubleshooting

- **Merge returns empty**: check that your reference CSVs contain the correct merge keys (`chr`, `pos`, `ref`, `alt`) and the SAD columns (`SAD<track>`), and that `reference_files_by_chr.json` points to the correct filenames.
- **PLINK errors / missing variants**: ensure `LD_REFERENCE_PATH` points to a valid PLINK binary dataset prefix and that alleles/rsIDs are compatible.
- **Missing phenotype mappings**: if you don’t use phenotype-id-based lookup, prefer passing `--filename` and keep `--phen` only for tracklist naming.