# DeepCAST-GWAS
## Running DeepCAST-FWER
### Preparing the data
You will need to load your data and then adjust the configuration file `src/config.py`.

You will also need to provide a tracklist for the phenotype you are looking at. The tracklist name has to match the filename and should be located in the `TRACKLISTS_DIR` as specified in the configuration file.
Support for phenotypes that are not in the PanUKBB sumstats is currently limited.

#### config parameters explained:
> SAD_DEVIATION_FACTOR  =  1

Cutoff for filtering SNPs based on their SAD scores for FWER control method.

>PVAL_THRESHOLD  =  5e-8

Conventional p-value threshold.

>FDR_THRESHOLD  =  0.0025

  FDR threshold for FDR control method.

>R2_THRESHOLD  =  0.2
KB_RADIUS  =  500
WINDOW_SIZE  =  500000

Parameters for PLINK clumping.
  
> CHROMOSOMES  =  list(range(1,  23))

  Chromosomes to consider.

>STD_BINS  =  [-np.inf,  -2,  -1,  -0.5,  -0.25,  0,  0.25,  0.5,  1,  2,  np.inf]
STD_BIN_LABELS  =  ["<-2σ",  "[-2σ, -1σ]",  "[-1σ, -0.5σ]",  "[-0.5σ, -0.25σ]",  "[-0.25σ, 0σ]",  "[0σ, 0.25σ]",  "[0.25σ, 0.5σ]",  "[0.5σ, 1σ]",  "[1σ, 2σ]",  ">2σ"]

Bins for stratified FDR control method.

You can download the reference genome merged with the enformer tracks [here](www.provide_download_link.com).

You can download the exon regions as pdf [here](www.provide_download_link.com).

### Calling the method
DeepCAST can be run by calling `src/main.py` (FWER control) or `src/deepcast_sfdr.py` (FDR control).

It has the following parameters:
>--phen

Phenotype id for recomputing a PanUKBB phenotype with the provided phenotype indexing (`data/phenotype_identifier.csv`)
>--run_id

provide a run_id to to create a result subfolder.
>--filename

When running phenotypes that are not indexed, the method will look for this filename in the `SUMSTATS_DIR` specified in the configuration file.
>--chr
--bp
--ref
--alt
--neglogp

Column names for the sumstats file can be specified in case they deviate from the deafault names specified in the configuration file.



> Written with [StackEdit](https://stackedit.io/).