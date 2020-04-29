#!/usr/bin/env python3

import sys
import glob

COJO_HEADER = 'SNP A1 A2 freq b se p n'.split()

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
    print('MAF threshold:', MAF_threshold)
    print('reading snp info from', filename)
    snp_info = {}
    n_skip = 0
    with gzip.open(filename, 'rt') as f:
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
            if 0 < MAF_threshold < 0.5 and (B_freq < MAF_threshold or B_freq > 1 - MAF_threshold):
                n_skip += 1
                continue
            if rsid.endswith(',.'):
                rsid = rsid[:-2]
            effect_allele = row['alleleB']
            other_allele = row['alleleA']
            eaf = B_freq
            n = int(AA_calls + AB_calls + BB_calls + 0.5)
            # see SMR documenation: https://cnsgenomics.com/software/smr/#SMR&HEIDIanalysis
            snp_info[rsid] = (effect_allele, other_allele, eaf, n)
    print(' -> discarded', n_skip, 'snps due to MAF threshold')
    return snp_info

def main():
    qtltools_nom_stderr_output_pattern = sys.argv[1]
    sumstats_file = sys.argv[2]
    outdir = sys.argv[3]
    MAF_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05

    snp_info = load_snp_info(sumstats_file, MAF_threshold)

    print('read', len(snp_info), 'snps')

    current = None
    lines = []

    # nominal:
    # 0          1  2        3        4 5     6      7            8  9       10      11          12       13           14          15
    # phenotype  region               S nsnp  dist   rsid         region             p_nom       r2       beta         se          best
    # cg17035109 21 10882030 10882030 - 15372 999554 21:9882476,. 21 9882476 9882476 1.31758e-19 0.170322 -2.02227e-15 2.12781e-16 0
    # cg17035109 21 10882030 10882030 - 15372 999499 21:9882531,. 21 9882531 9882531 1.31758e-19 0.170322 -1.90539e-15 2.00483e-16 0
    # cg17035109 21 10882030 10882030 - 15372 999474 21:9882556,. 21 9882556 9882556 1.31758e-19 0.170322 -1.8952e-15  1.99411e-16 0

    # permutatino:
    # 0                1   2         3         4  5      6        7            8   9         10        11   12       13        14           15           16         17       18       19         20
    # pid              pc  pb        pe        ps nsnps  dist     vid          vc  vb        ve        df   dummy    par1      par2(nindep) pnom         r2         slope    se       pval_emp   pval_bml     
    # ENSG00000099991  22  24407643  24574596  +  50113  -457933  22:23949710  22  23949710  23949710  590  106.909  1.2459    64582.1      1.97667e-23  0.155288   34481.6  3310.91  0.766234   0.683884
    # ENSG00000099998  22  24615623  24641110  -  47107  435173   22:24180450  22  24180450  24180450  590  45.4023  3.1737    57319.5      9.19196e-49  0.306051   224015   13887.3  0.681319   0.517512
    # ENSG00000128262  22  24647797  24661493  +  46806  -467347  22:24180450  22  24180450  24180450  590  29.5432  0.998721  516.086      9.01812e-45  0.284178   34774.3  2272.16  0.687313   0.611355

    nfiles = 0
    nlines = 0
    nmissing = 0

    def submit():
        nonlocal nfiles, nlines, nmissing
        if lines:
            nfiles += 1
            with open(outdir + '/' + current, 'w') as f:
                print(*COJO_HEADER, file=f)
                for line in lines:
                    parts = line.split()
                    snp_id = parts[7]
                    if snp_id.endswith(',.'):
                        snp_id_test = snp_id[:-2]
                    else:
                        snp_id_test = snp_id
                    if snp_id_test not in snp_info:
                        nmissing += 1
                        continue
                    if len(parts) == 16: # nominal
                        snp_chr = parts[8]
                        snp_pos = parts[9]
                        snp_p = parts[11]
                        snp_b = parts[13]
                        snp_se = parts[14]
                    elif len(parts) == 21: # permutation
                        snp_chr = parts[8]
                        snp_pos = parts[9]
                        snp_p = parts[20]
                        snp_b = parts[17]
                        snp_se = parts[18]
                    else:
                        print('unknown qtltools output format')
                        exit(1)
                    A1, A2, A1freq, n = snp_info[snp_id_test]
                    print(snp_id, A1, A2, A1freq, snp_b, snp_se, snp_p, n, file=f)
                    nlines += 1
        del lines[:]

    for filename in glob.glob(qtltools_nom_stderr_output_pattern):
        with open(filename) as f:
            for line in f:
                probe, rest = line.split(None, 1)
                if probe != current:
                    submit()
                    current = probe
                lines.append(line.rstrip())
        submit()
        print(filename, 'wrote', nlines, 'qtls in', nlines, 'cojo files', nmissing, 'discarded(maf)')

if __name__ == '__main__':
    main()
