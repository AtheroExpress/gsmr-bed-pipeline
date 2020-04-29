import os
import sys
import textwrap

import job_definitions
import config

def main():

    config.load_section('DEFAULT')
    for section in sys.argv[1:]:
        config.load_section(section)

    os.makedirs(os.path.join(config.job_directory, 'excl_cov'), exist_ok=True)

    config.excl_cov_exposure_file = os.path.join(config.job_directory, 'excl_cov', 'exposure.excl')
    config.excl_cov_outcome_file = os.path.join(config.job_directory, 'excl_cov', 'outcome.excl')

    with open(config.excl_cov_exposure_file, 'w') as f:
        for cov in config.exclude_covariates_exposure.split():
            print(cov, file=f)
    with open(config.excl_cov_outcome_file, 'w') as f:
        for cov in config.exclude_covariates_outcome.split():
            print(cov, file=f)

    for region in config.regions.split():
        gsmr_filt_file = os.path.join(config.job_directory, region + '.gsmr')
        if os.path.exists(gsmr_filt_file):
            print('# already exists', region)
        region_jobs = job_definitions.jobs_for_region(region)
        regiondir = os.path.join(config.job_directory, region)
        os.makedirs(regiondir, exist_ok=True)
        submit_file = os.path.join(config.job_directory, region, 'submit.sh')

        with open(submit_file, 'w') as f:
            print(textwrap.dedent(rf'''
            #!/usr/bin/env bash
            set -ex
            ''').strip(), file=f)
            for job in region_jobs:
                jobname = job['jobname']
                if job['depends']:
                    afterok = ':'.join('${'+dep+'}' for dep in job['depends'])
                    dep_arg = '-d afterok:' + afterok
                else:
                    dep_arg = ''
                jobfile = submit_file + '.' + jobname
                print(f'output=$(sbatch {dep_arg} -J {jobname}{region} -e {regiondir}/{jobname}.submit.out -o {regiondir}/{jobname}.submit.out {jobfile})', file=f)
                with open(jobfile, 'w') as g:
                    print(job['job'], file=g)
                print('echo $output', file=f)
                print(f'{jobname}=${{output##* }}', file=f)
                print(file=f)

        print('bash', submit_file)

    # awk 'NR==1{print "FILENAME", $0} /^cg/{if($5<1e-5){print FILENAME, $0}} ' 2020-04-06-max-power/*.gsmr | cut -f 2 -d / | sort -n

if __name__ == '__main__':
    main()
