import os
import textwrap

import job_definitions
import config

def main():
    for region in config.regions.split():
        print('# Generating jobs for region', region)
        region_jobs = job_definitions.jobs_for_region(region)
        os.makedirs(os.path.join(config.job_directory, region), exist_ok=True)
        submit_file = os.path.join(config.job_directory, region, 'submit.sh')

        with open(submit_file, 'w') as f:
            print(textwrap.dedent(r'''
            #!/usr/bin/env bash
            #SBATCH --time 15:00:00
            #SBATCH --mem 40G

            set -ex
            ''').strip(), file=f)
            for job in region_jobs:
                print(job, file=f)

        print('sbatch', submit_file)

if __name__ == '__main__':
    main()
