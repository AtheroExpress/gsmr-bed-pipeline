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

 - Edit `pipeline.config`.
 - `python3 main.py`
 - Execute printed bash statements to submit jobs
