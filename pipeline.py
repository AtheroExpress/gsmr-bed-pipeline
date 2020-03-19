import job_definitions
import config

region_jobs = job_definitions.jobs_for_region(config.qtl_region)

print(r'''#!/usr/bin/env bash
#SBATCH --time 15:00:00
#SBATCH --mem 40G

set -ex
''')

for job in region_jobs:
    print(job)
