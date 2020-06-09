# GSMR CpG-mRNA pipeline

This pipeline researches the causal relation between
two molecular phenotypes. 
For example, the relation between methylation and gene expression
from CpG beta-values and RNAseq counts.

## Pipeline steps

  - Calculation of eQTLs using a patched QTLtools v1.2
  - Calculation of mQTLs using a patched QTLtools v1.2
  - Conversion of eQTLs to SMR-COJO files & threshold MAF
  - Conversion of mQTLs to SMR-COJO files & threshold MAF
  - Optionally: Finding the optimal configuration of GSMR calls
  - Creating the GSMR exposure & outcomes files
  - Running GSMR in bi-GSMR mode
  - Summarize GSMR output & Plotting

## Usage:

### Edit `pipeline.config`.

It is recommended to keep your config under a separate profile, instead of editing the `DEFAULT` settings:

```
[peer]
    covariance_exposure = /hpc/dhl_ec/llandsmeer/_atheroexpress/peer/full40/COVARIATES40MFULL
    covariance_outcome = /hpc/dhl_ec/llandsmeer/_atheroexpress/peer_eqtl/full70/PEER70RNASEQ
    job_directory = 2020-06-09-peer-opt
    qtl_mode = nominal
    exclude_covariates_exposure =
    exclude_covariates_outcome =
```

### Execute `python3 main.py <your profile>`

This will generate the jobfiles, and print bash commands to execute submit those jobs
on the command line. Leave `<your profile>` empty to create the `DEFAULT` job.

### Submit jobs

Execute the output of `main.py`. For a test run, you might want to run only a single
line instead of all output.
