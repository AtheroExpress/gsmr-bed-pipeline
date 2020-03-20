import os
import textwrap
import config

def inject(region, paths, fmt):
    chrom = region.split(':')[0]
    def register_path(fmt):
        path = fmt.format(config.job_directory, region)
        if path not in paths:
            paths.append(path)
        return path
    fmt = textwrap.dedent(fmt).lstrip('\n')
    return fmt.format(
        exposure_bed = config.exposure_bed.format(chr=chrom),
        outcome_bed = config.outcome_bed.format(chr=chrom),
        exposure_qtl = register_path('{0}/{1}/exposure/exposure.nom.txt'),
        outcome_qtl = register_path('{0}/{1}/outcome/outcome.nom.txt'),
        exposure_cojo_dir = register_path('{0}/{1}/exposure/exposure.cojo/'),
        outcome_cojo_dir = register_path('{0}/{1}/outcome/outcome.cojo/'),
        gsmr_exposure = register_path('{0}/{1}/gsmr/gsmr_exposure.txt'),
        gsmr_outcome = register_path('{0}/{1}/gsmr/gsmr_outcome.txt'),
        gsmr_out = register_path('{0}/{1}/gsmr/gsmr.txt'),
        gsmr_out_filtered = register_path('{0}/{1}.gsmr'),
        gsmr_plot_dir = register_path('{0}/{1}/plot/'),
        gen_bed = register_path('{0}/{1}/gsmr/bed'),
        covariance = config.covariance,
        vcf = config.vcf_per_chr.format(chr=chrom),
        job_directory = config.job_directory,
        qtl_nom_pvalue = config.qtl_nom_pvalue,
        qtl_window = config.qtl_window,
        qtl_seed = config.qtl_seed,
        region = region if ':' in region else False,
        software_rscript = config.software_rscript,
        arg_QTLtools_region = '--region "{0}"'.format(region) if ':' in region else '',
        MAF_threshold = float(config.maf_threshold)
        )

def jobs_for_region(region):
    paths = []

    qtltools_exposure = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 10:00:00
    #SBATCH --mem 40G

    if [ ! -f "{exposure_qtl}" ]; then
        vendor/qtltools_v1.2-stderr cis \
                --nominal   "{qtl_nom_pvalue}" \
                --vcf       "{vcf}" \
                --bed       "{exposure_bed}" \
                --cov       "{covariance}" \
                --out       "{exposure_qtl}" \
                --window    "{qtl_window}" \
                --seed      "{qtl_seed}" \
                --std-err {arg_QTLtools_region}
    fi

    python3 _scripts/split_qtl_to_cojo.py \
            "{exposure_qtl}" \
            "{exposure_cojo_dir}" \
            "{MAF_threshold}"

    find "{exposure_cojo_dir}" -type f \
            | awk -F / '{{print $NF " " $0}}' \
            > "{gsmr_exposure}"
    ''')


    qtltools_outcome = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 10:00:00
    #SBATCH --mem 40G

    if [ ! -f "{outcome_qtl}" ]; then
        vendor/qtltools_v1.2-stderr cis \
                --nominal   "{qtl_nom_pvalue}" \
                --vcf       "{vcf}" \
                --bed       "{outcome_bed}" \
                --cov       "{covariance}" \
                --out       "{outcome_qtl}" \
                --window    "{qtl_window}" \
                --seed      "{qtl_seed}" \
                --std-err {arg_QTLtools_region}
    fi

    python3 _scripts/split_qtl_to_cojo.py \
            "{outcome_qtl}" \
            "{outcome_cojo_dir}" \
            "{MAF_threshold}"

    find "{outcome_cojo_dir}" -type f \
            | awk -F / '{{print $NF " " $0}}' \
            > "{gsmr_outcome}"
    ''')

    create_bed = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 1:00:00
    #SBATCH --mem 40G

    #!/usr/bin/env bash
    plink2 --make-bed \
            --vcf "{vcf}" \
            --out "{gen_bed}"
    ''')

    run_gsmr = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 6:00:00
    #SBATCH --mem 40G

    if [ ! -f "{gsmr_out}.gsmr" ]; then
        gcta_1.92.1b6 \
            --bfile "{gen_bed}" \
            --gsmr-file "{gsmr_exposure}" "{gsmr_outcome}" \
            --gsmr-direction 2 \
            --out "{gsmr_out}" \
            --gwas-thresh 0.01 \
            --effect-plot \
            --clump-r2 0.1
    fi

    cat "{gsmr_out}.gsmr" \
            | grep -v 'nan.*nan.*nan.*nan' \
            | column -t \
            > "{gsmr_out_filtered}"

    python3 _scripts/summarize_gsmr.py \
            "{gsmr_out_filtered}" \
            "{exposure_qtl}" \
            "{outcome_qtl}"\
            | column -t \
            > "{gsmr_out_filtered}.summary"
    ''')

    plot_gsmr = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 8:00:00
    #SBATCH --mem 40G

    "{software_rscript}" _scripts/run_gsmr_plot.r \
        "{gsmr_out}.eff_plot.gz" \
        "{gsmr_out_filtered}" \
        "{gsmr_plot_dir}"
    ''')

    basedirs = sorted(list(set(map(os.path.dirname, paths))))

    create_dirs = '\n'.join('mkdir -p {0}'.format(basedir) for basedir in basedirs) + '\n'
    create_dirs = textwrap.dedent(r'''
    #!/usr/bin/env bash
    {0}
    '''.lstrip('\n')).format(create_dirs)

    return [
        create_dirs,
        qtltools_exposure,
        qtltools_outcome,
        create_bed,
        run_gsmr,
        plot_gsmr
    ]
