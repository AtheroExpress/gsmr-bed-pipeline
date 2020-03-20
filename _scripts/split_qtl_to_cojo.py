#!/usr/bin/env python3

import sys

COJO_HEADER = 'SNP A1 A2 freq b se p n'.split()

sumstats_chr21 = '/hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_EAGLE2_1000Gp3v5HRCr11/aegs.qc.1kgp3hrcr11.chr21.sumstats.txt.gz'

def load_snp_info(filename, MAF_threshold=0.05):
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
            if 0 < MAF_threshold < 0.5 and (B_freq < MAF_threshold or B_freq > 1 - MAF_threshold):
                continue
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
    MAF_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05

    snp_info = load_snp_info(sumstats_chr21, MAF_threshold)

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
