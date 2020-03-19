import os
import textwrap
import config

def inject(paths, fmt):
    chrom = config.qtl_region.split(':')[0]
    def register_path(fmt):
        path = fmt.format(config.job_directory, chrom)
        if path not in paths:
            paths.append(path)
        return path
    fmt = textwrap.dedent(fmt).lstrip('\n')
    return fmt.format(
        exposure_bed = config.exposure_bed,
        outcome_bed = config.outcome_bed,
        exposure_qtl = register_path('{0}/{1}/exposure/exposure.nom.txt'),
        outcome_qtl = register_path('{0}/{1}/outcome/outcome.nom.txt'),
        exposure_cojo_dir = register_path('{0}/{1}/exposure/exposure.cojo/'),
        outcome_cojo_dir = register_path('{0}/{1}/outcome/outcome.cojo/'),
        gsmr_exposure = register_path('{0}/{1}/gsmr/gsmr_exposure.txt'),
        gsmr_outcome = register_path('{0}/{1}/gsmr/gsmr_outcome.txt'),
        gsmr_out = register_path('{0}/{1}/gsmr/gsmr.txt'),
        gsmr_out_filtered = register_path('{0}/{1}.gsmr'),
        gen_bed = register_path('{0}/{1}/gsmr/bed'),
        covariance = config.covariance,
        vcf = config.vcf_per_chr.format(chr=chrom),
        job_directory = config.job_directory,
        qtl_nom_pvalue = config.qtl_nom_pvalue,
        qtl_window = config.qtl_window,
        qtl_seed = config.qtl_seed,
        qtl_region = config.qtl_region)


def jobs_for_region(region):
    paths = []

    qtltools_exposure = inject(paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 10:00:00
    #SBATCH --mem 40G

    vendor/qtltools_v1.2-stderr cis \
            --nominal   "{qtl_nom_pvalue}" \
            --vcf       "{vcf}" \
            --bed       "{exposure_bed}" \
            --cov       "{covariance}" \
            --out       "{exposure_qtl}" \
            --window    "{qtl_window}" \
            --seed      "{qtl_seed}" \
            --std-err \
            --region    "{qtl_region}"

    python3 _scripts/split_qtl_to_cojo.py \
            "{exposure_qtl}" \
            "{exposure_cojo_dir}"

    find "{exposure_cojo_dir}" -type f \
            | awk -F / '{{print $NF " " $0}}' \
            > "{gsmr_exposure}"
    ''')


    qtltools_outcome = inject(paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 10:00:00
    #SBATCH --mem 40G

    vendor/qtltools_v1.2-stderr cis \
            --nominal   "{qtl_nom_pvalue}" \
            --vcf       "{vcf}" \
            --bed       "{outcome_bed}" \
            --cov       "{covariance}" \
            --out       "{outcome_qtl}" \
            --window    "{qtl_window}" \
            --seed      "{qtl_seed}" \
            --std-err \
            --region    "{qtl_region}"

    python3 _scripts/split_qtl_to_cojo.py \
            "{outcome_qtl}" \
            "{outcome_cojo_dir}"

    find "{outcome_cojo_dir}" -type f \
            | awk -F / '{{print $NF " " $0}}' \
            > "{gsmr_outcome}"
    ''')

    create_bed = inject(paths, r'''
    #SBATCH --time 1:00:00
    #SBATCH --mem 40G

    #!/usr/bin/env bash
    plink2 --make-bed \
            --vcf "{vcf}" \
            --out "{gen_bed}"
    ''')

    run_gsmr = inject(paths, r'''
    #SBATCH --time 6:00:00
    #SBATCH --mem 40G

    #!/usr/bin/env bash
    gcta_1.92.1b6 \
        --bfile "{gen_bed}" \
        --gsmr-file "{gsmr_exposure}" "{gsmr_outcome}" \
        --gsmr-direction 2 \
        --out "{gsmr_out}" \
        --gwas-thresh 0.01 \
        --effect-plot \
        --clump-r2 0.1

    cat "{gsmr_out"} \
            | grep -v 'nan.*nan.*nan.*nan' \
            | column -t \
            > "{gsmr_out_filtered}"

    ''')

    basedirs = sorted(list(set(map(os.path.dirname, paths))))

    create_dirs = '\n'.join('mkdir -p {0}'.format(basedir) for basedir in basedirs) + '\n'
    create_dirs = textwrap.dedent(r'''
    #!/usr/bin/env bash
    {0}
    ''').format(create_dirs).lstrip('\n')

    return [
        create_dirs,
        qtltools_exposure,
        qtltools_outcome,
        create_bed,
        run_gsmr
    ]
