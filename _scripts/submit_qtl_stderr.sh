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
