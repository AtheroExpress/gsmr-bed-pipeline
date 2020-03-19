# GSMR CpG-mRNA pipeline

This pipeline researches the relation between methylation
and mRNA expression using GSMR.

The B allele will be

## Pipeline steps

  - Calculation of eQTLs using a patched QTLtools v1.2
  - Calculation of mQTLs using a patched QTLtools v1.2
  - Conversion of eQTLs to COJO files & threshold MAF
  - Conversion of mQTLs to COJO files
  - Optionally: Finding the optimal configuration of GSMR calls
  - Creating the GSMR exposure & outcomes files
  - Running GSMR in bi- mode
  - Tabulating GSMR output + SNPs & Plotting

## Todo:

  - Implement MAF filter
  - Some freedom in sumstats specification
  - Find correct p-value threshold
  - Implement user-specific configs
  - Implement covariates exclusion

## eQTL discovery

```
#!/usr/bin/env bash

#SBATCH --time 10:00:00
#SBATCH --mem 40G

CHR=21
MODE=eqtl

BED_RNASEQ="/hpc/dhl_ec/llandsmeer/_atheroexpress/_integrative_v2/_data/2020-03-18-rnaseq-15710x624.qc.bed.gz"
BED_AEMS450K1="/hpc/dhl_ec/llandsmeer/_atheroexpress/_integrative_v2/_data/2020-03-18-aems450k1.bvalues.plaque.qtl.bed.gz"
COV="/hpc/dhl_ec/llandsmeer/_atheroexpress/_integrative_v2/_data/phenocov.aedb"
VCF="/hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11/aegs.qc.1kgp3hrcr11.chr${CHR}.vcf.gz"
# VCF=./region.vcf.gz

if [ "$MODE" = eqtl ]; then
    BED="$BED_RNASEQ"
    FILEOUT="$(date -Id).eqtl.chr${CHR}.se.nom"
elif [ "$MODE" = mqtl ]; then
    BED="$BED_AEMS450K1"
    FILEOUT="$(date -Id).mqtl.chr${CHR}.se.nom"
else
    echo "unknown mode"
    exit 2
fi

FILEOUT="${FILEOUT}.region"

qtltools_v1.2-stderr cis \
        --nominal   0.0001 \
        --vcf       "${VCF}" \
        --bed       "${BED}" \
        --cov       "${COV}" \
        --out       "${FILEOUT}" \
        --window    1000000 \
        --seed      12421 \
        --std-err \
        --region    21:9500000-40500000
```

## Split to COJO

```python
#!/usr/bin/env python3

import sys

COJO_HEADER = 'SNP A1 A2 freq b se p n'.split()

sumstats_chr21 = '/hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11/aegs.qc.1kgp3hrcr11.chr21.sumstats.txt.gz'

def load_snp_info(filename):
    import gzip
    '''load SNP summary statistics info

    Reads the summary statistics file <filename> into
    the dictionary { RSID => (Effect allele, Other allele, EAF, Nsamples) }

    B will be the effect allele.

    Assumes that the relevant headers in the summary
    statistics file are called:
        - 'rsid', 'alleleA', 'alleleB' and 'all_maf'
        - variant call check: 'all_AA', 'all_BB', 'all_AB'

    '''
    snp_info = {}
    with gzip.open(sumstats_chr21, 'rt') as f:
        header = None
        for line in f:
            if line.startswith('#'):
                continue
            elif header is None:
                header = line.split()
                continue
            row = dict(zip(header, line.split()))
            AA_calls = float(row['all_AA'])
            AB_calls = float(row['all_AB'])
            BB_calls = float(row['all_BB'])
            B_freq = (2*BB_calls + AB_calls) / (2*AA_calls + 2*AB_calls + 2*BB_calls)
            rsid = row['rsid']
            effect_allele = row['alleleB']
            other_allele = row['alleleA']
            eaf = B_freq
            n = int(AA_calls + AB_calls + BB_calls + 0.5)
            # see SMR documenation: https://cnsgenomics.com/software/smr/#SMR&HEIDIanalysis
            snp_info[rsid] = (effect_allele, other_allele, eaf, n)
    return snp_info

def main():
    qtltools_nom_stderr_output = sys.argv[1]
    outdir = sys.argv[2]

    snp_info = load_snp_info(sumstats_chr21)

    current = None
    lines = []

    # 0          1  2        3        4 5     6      7            8  9       10      11          12       13           14          15
    # phenotype  region               S nsnp  dist   rsid         region             p_nom       r2       beta         se          best
    # cg17035109 21 10882030 10882030 - 15372 999554 21:9882476,. 21 9882476 9882476 1.31758e-19 0.170322 -2.02227e-15 2.12781e-16 0
    # cg17035109 21 10882030 10882030 - 15372 999499 21:9882531,. 21 9882531 9882531 1.31758e-19 0.170322 -1.90539e-15 2.00483e-16 0
    # cg17035109 21 10882030 10882030 - 15372 999474 21:9882556,. 21 9882556 9882556 1.31758e-19 0.170322 -1.8952e-15  1.99411e-16 0

    nfiles = 0

    def submit():
        nonlocal nfiles
        #if nfiles > 50: exit(0)
        if lines:
            nfiles += 1
            with open(outdir + '/' + current, 'w') as f:
                print(*COJO_HEADER, file=f)
                for line in lines:
                    parts = line.split()
                    assert len(parts) == 16
                    snp_id = parts[7]
                    if snp_id.endswith(',.'):
                        snp_id = snp_id[:-2]
                    if snp_id not in snp_info:
                        #print('ERR NOT FOUND', snp_id)
                        continue
                    snp_chr = parts[8]
                    snp_pos = parts[9]
                    snp_p = parts[11]
                    snp_b = parts[13]
                    snp_se = parts[14]
                    A1, A2, A1freq, n = snp_info[snp_id]
                    print(snp_id+',.', A1, A2, A1freq, snp_b, snp_se, snp_p, n, file=f)
        del lines[:]

    with open(qtltools_nom_stderr_output) as f:
        for line in f:
            probe, rest = line.split(None, 1)
            if probe != current:
                submit()
                current = probe
            lines.append(line.rstrip())
    submit()

if __name__ == '__main__':
    main()

```

## Create BED files

```
plink --vcf /hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11/aegs.qc.1kgp3hrcr11.chr21.vcf.gz --out chr21 --make-bed
```

## GSMR Exposure & Outcomes files

```
find _mqtl/ -type f | awk -F / '{print $2 " " $0}' > _data/gsmr_exposure.txt
find _eqtl/ -type f | awk -F / '{print $2 " " $0}' > _data/gsmr_outcome.txt
```

## Run GSMR

```
gcta_1.92.1b6 \
    --bfile _data/chr21 \
    --gsmr-file _data/gsmr_exposure.txt _data/gsmr_outcome.txt \
    --gsmr-direction 2 \
    --out gsmr_result \
    --gwas-thresh 0.01 \
    --effect-plot \
    --clump-r2 0.1
```

## Summarize results

```
$ cat gsmr_result.gsmr  | grep -v 'nan.*nan.*nan.*nan' | column -t
Exposure    Outcome          bxy       se       p            nsnp
cg01850767  ENSG00000159082  -1769.85  143.517  6.09103e-35  8
cg01850767  ENSG00000142168  -42610.7  2767.6   1.73445e-53  9
```

## Plot

```
source("gsmr_plot.r")
gsmr_data = read_gsmr_data("./gsmr_result.eff_plot.gz")
gsmr_summary(gsmr_data)
plot_gsmr_effect(gsmr_data, "cg01850767", "ENSG00000159082")
plot_gsmr_effect(gsmr_data, "cg01850767", "ENSG00000142168")
```
