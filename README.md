# GSMR CpG-mRNA pipeline

This pipeline researches the relation between methylation
and mRNA expression using GSMR.

## Pipeline steps

  - Calculation of eQTLs using a patched QTLtools v1.2
  - Calculation of mQTLs using a patched QTLtools v1.2
  - Conversion of eQTLs to COJO files & threshold MAF
  - Conversion of mQTLs to COJO files
  - Optionally: Finding the optimal configuration of GSMR calls
  - Creating the GSMR exposure & outcomes files
  - Running GSMR in bi- mode
  - Tabulating GSMR output + SNPs & Plotting
