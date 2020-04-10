#!/usr/bin/env python3

import sys

COJO_HEADER = 'SNP A1 A2 freq b se p n'.split()

exposure_file = sys.argv[1] # '2020-03-25-first-gw-run/9/gsmr/gsmr_exposure.txt'
outcome_file = sys.argv[2] # '2020-03-25-first-gw-run/9/gsmr/gsmr_outcome.txt'

def load(filename):
    pheno_to_rsids = {}
    with open(filename) as f:
        for sz, line in enumerate(f):
            pass
    lastprog = -10
    try:
        with open(filename) as f:
            for idx, line in enumerate(f):
                phenotype, cojo_file = line.split()
                rsids = []
                with open(cojo_file) as g:
                    assert(next(g).split() == COJO_HEADER)
                    for line in g:
                        rsid = line.split(None, maxsplit=1)[0]
                        rsids.append(rsid)
                pheno_to_rsids[phenotype] = (cojo_file, set(rsids))
                prog = int(round(idx*100/sz))
                if prog >= lastprog + 5:
                    print(' - read', idx, 'lines', prog, '%', file=sys.stderr)
                    lastprog = prog
    except KeyboardInterrupt:
        print(' - interrupted', file=sys.stderr)
    return pheno_to_rsids

print('loading exposure', file=sys.stderr)
E = load(exposure_file)

print('loading outcome', file=sys.stderr)
O = load(outcome_file)

print('done loading', file=sys.stderr)

for exposure_name, (exposure_cojo, exposure_rsids) in E.items():
    for outcome_name, (outcome_cojo, outcome_rsids) in O.items():
        common_rsids = exposure_rsids & outcome_rsids
        if len(common_rsids) > 7:
            print(exposure_name, outcome_name, exposure_cojo, outcome_cojo, sep='\t')
