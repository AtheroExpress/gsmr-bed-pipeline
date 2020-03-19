import job_definitions
import config

region_jobs = job_definitions.jobs_for_region(config.qtl_region)

print('set -ex')

for job in region_jobs:
    print(job)
