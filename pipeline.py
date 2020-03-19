import os

import job_definitions
import config

region_jobs = job_definitions.jobs_for_region(config.qtl_region)

os.makedirs(config.job_directory, exist_ok=True)
submit_file = os.path.join(config.job_directory, 'submit.sh')

with open(submit_file, 'w') as f:
    print(r'''#!/usr/bin/env bash
    #SBATCH --time 15:00:00
    #SBATCH --mem 40G

    set -ex
    ''', file=f)

    for job in region_jobs:
        print(job, file=f)

print('sbatch', submit_file)
