#!/usr/bin/env bash

#python3 to_cojo.py 2020-03-18.mqtl.se.nom _mqtl/
#python3 to_cojo.py 2020-03-18.eqtl.se.nom _eqtl/

# python3 to_cojo.py 2020-03-18.mqtl.chr21.se.nom _mqtl/
#python3 to_cojo.py 2020-03-18.eqtl.chr21.se.nom _eqtl/

# plink --vcf /hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11/aegs.qc.1kgp3hrcr11.chr21.vcf.gz --out chr21 --make-bed

find _mqtl/ -type f | awk -F / '{print $2 " " $0}' > _data/gsmr_exposure.txt
find _eqtl/ -type f | awk -F / '{print $2 " " $0}' > _data/gsmr_outcome.txt

gcta_1.92.1b6 \
    --bfile _data/chr21 \
    --gsmr-file _data/gsmr_exposure.txt _data/gsmr_outcome.txt \
    --gsmr-direction 0 \
    --out gsmr_result
