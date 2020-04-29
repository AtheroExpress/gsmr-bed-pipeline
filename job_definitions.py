import os
import textwrap
import config

class DoNotUseThisSetting():
    def __repr__(self):
        raise Exception('DoNotUseThisSetting')

def inject(region, paths, fmt):
    assert config.qtl_mode in ['nominal', 'permutation']
    chrom = region.split(':')[0]
    def register_path(fmt):
        path = fmt.format(config.job_directory, region)
        if path not in paths:
            paths.append(path)
        return path
    def intge(var, n):
        assert int(var) >= n
        return int(var)
    fmt = textwrap.dedent(fmt).lstrip('\n')
    return fmt.format(
        exposure_bed = config.exposure_bed.format(chr=chrom),
        outcome_bed = config.outcome_bed.format(chr=chrom),
        exposure_qtl = register_path('{0}/{1}/exposure/exposure.nom.txt'),
        outcome_qtl = register_path('{0}/{1}/outcome/outcome.nom.txt'),
        exposure_cojo_dir = register_path('{0}/{1}/exposure/exposure.cojo/'),
        outcome_cojo_dir = register_path('{0}/{1}/outcome/outcome.cojo/'),
        gsmr_combinations = register_path('{0}/{1}/gsmr/gsmr_combinations.txt'),
        gsmr_exposure = register_path('{0}/{1}/gsmr/gsmr_exposure.txt'),
        gsmr_outcome = register_path('{0}/{1}/gsmr/gsmr_outcome.txt'),
        gsmr_out_dir = register_path('{0}/{1}/gsmr/combinations/'),
        gsmr_out_filtered = register_path('{0}/{1}.gsmr'),
        gsmr_out=DoNotUseThisSetting(),
        gsmr_plot_dir = register_path('{0}/{1}/plot/'),
        gen_bed = register_path('{0}/{1}/gsmr/bed'),
        covariance = config.covariance,
        vcf = config.vcf_per_chr.format(chr=chrom),
        sumstats = config.sumstats_per_chr.format(chr=chrom),
        job_directory = config.job_directory,
        qtl_nom_pvalue = config.qtl_nom_pvalue if config.qtl_mode == 'nominal' else DoNotUseThisSetting(),
        qtl_permutations = config.qtl_permutations if config.qtl_mode == 'permutation' else DoNotUseThisSetting(),
        qtl_mode = config.qtl_mode,
        qtl_window = config.qtl_window,
        qtl_seed = config.qtl_seed,
        #region = region if ':' in region else '{0:02d}'.format(int(chrom)), # hacky
        software_rscript = config.software_rscript,
        MAF_threshold_exposure = float(config.maf_threshold_exposure),
        MAF_threshold_outcome = float(config.maf_threshold_outcome),
        excl_cov_exposure_file = config.excl_cov_exposure_file,
        excl_cov_outcome_file = config.excl_cov_outcome_file,
        gsmr_r2 = config.gsmr_r2,
        gsmr_p = config.gsmr_p,
        arg_qtltools_mode = "--nominal " + str(config.qtl_nom_pvalue)
            if config.qtl_mode == 'nominal' else "--permute " + str(config.qtl_permutations),
        qtl_jobs=intge(config.qtl_jobs, 1),
        gsmr_max_job_idx=intge(config.gsmr_jobs, 1)-1,
        )

def jobs_for_region(region):
    paths = []

    qtltools_exposure = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 35:00:00
    #SBATCH --mem 5G
    #SBATCH --array=1-{qtl_jobs}

    if [ ! -f "{exposure_qtl}.${{SLURM_ARRAY_TASK_ID}}" ]; then
        vendor/qtltools_v1.2-stderr cis \
                {arg_qtltools_mode} \
                --vcf       "{vcf}" \
                --bed       "{exposure_bed}" \
                --cov       "{covariance}" \
                --out       "{exposure_qtl}".${{SLURM_ARRAY_TASK_ID}} \
                --window    "{qtl_window}" \
                --exclude-covariates "{excl_cov_outcome_file}" \
                --seed      "{qtl_seed}" \
                --std-err \
                --chunk ${{SLURM_ARRAY_TASK_ID}} ${{SLURM_ARRAY_TASK_COUNT}}
    fi
    ''')

    qtltools_exposure_collect = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 4:00:00
    #SBATCH --mem 20G

    if [ -z "$(ls -A {exposure_cojo_dir} | head -n 1)" ]; then
        python3 _scripts/split_qtl_to_cojo.py \
                "{exposure_qtl}".'*' \
                "{sumstats}" \
                "{exposure_cojo_dir}" \
                "{MAF_threshold_exposure}"
    fi

    find "{exposure_cojo_dir}" -type f \
            | awk -F / '{{print $NF " " $0}}' \
            > "{gsmr_exposure}"
    ''')

    qtltools_outcome = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 35:00:00
    #SBATCH --mem 5G
    #SBATCH --array=1-{qtl_jobs}

    if [ ! -f "{outcome_qtl}.${{SLURM_ARRAY_TASK_ID}}" ]; then
        vendor/qtltools_v1.2-stderr cis \
                {arg_qtltools_mode} \
                --vcf       "{vcf}" \
                --bed       "{outcome_bed}" \
                --cov       "{covariance}" \
                --out       "{outcome_qtl}".${{SLURM_ARRAY_TASK_ID}} \
                --window    "{qtl_window}" \
                --exclude-covariates "{excl_cov_exposure_file}" \
                --seed      "{qtl_seed}" \
                --std-err \
                --chunk ${{SLURM_ARRAY_TASK_ID}} ${{SLURM_ARRAY_TASK_COUNT}}
    fi
    ''')

    qtltools_outcome_collect = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 2:00:00
    #SBATCH --mem 20G

    if [ -z "$(ls -A {outcome_cojo_dir} | head -n 1)" ]; then
        python3 _scripts/split_qtl_to_cojo.py \
                "{outcome_qtl}".'*' \
                "{sumstats}" \
                "{outcome_cojo_dir}" \
                "{MAF_threshold_outcome}"
    fi

    find "{outcome_cojo_dir}" -type f \
            | awk -F / '{{print $NF " " $0}}' \
            > "{gsmr_outcome}"
    ''')

    create_bed = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 1:00:00
    #SBATCH --mem 40G

    #!/usr/bin/env bash
    if [ ! -f "{gen_bed}.bed" ]; then
        plink2 --make-bed \
                --vcf "{vcf}" \
                --out "{gen_bed}"
    fi
    ''')

    gsmr_pairs = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 15:00:00
    #SBATCH --mem 60G

    if [ ! -f "{gsmr_combinations}" ]; then
        python3 _scripts/gsmr_pairs.py \
                "{gsmr_exposure}" \
                "{gsmr_outcome}" \
                > "{gsmr_combinations}"
    fi
    ''')

    run_gsmr = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 2:00:00
    #SBATCH --mem 40G
    #SBATCH --array 0-{gsmr_max_job_idx}

    cat "{gsmr_combinations}" \
        | awk "(NR-1)%${{SLURM_ARRAY_TASK_COUNT}} == ${{SLURM_ARRAY_TASK_ID}}" \
        | while read line; do
        if [ ! -f "{gsmr_out_dir}"/"${{combination}}.log" ]; then
            combination=$(echo "$line" | awk '{{print $1 "-" $2}}')
            exposure=$(echo "$line" | awk '{{print $1}}')
            outcome=$(echo "$line" | awk '{{print $2}}')
            exposure_file=$(echo "$line" | awk '{{print $3}}')
            outcome_file=$(echo "$line" | awk '{{print $4}}')
            echo "${{exposure}}" "${{exposure_file}}" > "{gsmr_out_dir}/${{combination}}.gsmr_exposure.txt"
            echo "${{outcome}}" "${{outcome_file}}" > "{gsmr_out_dir}/${{combination}}.gsmr_outcome.txt"
            gcta_1.92.1b6 \
                --bfile "{gen_bed}" \
                --gsmr-file \
                    "{gsmr_out_dir}/${{combination}}.gsmr_exposure.txt" \
                    "{gsmr_out_dir}/${{combination}}.gsmr_outcome.txt" \
                --gsmr-direction 2 \
                --out "{gsmr_out_dir}"/"${{combination}}" \
                --gwas-thresh {gsmr_p} \
                --effect-plot \
                --clump-r2 {gsmr_r2} || true
        fi
    done
    ''')

    gsmr_summarize = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 8:00:00
    #SBATCH --mem 40G

    (
    cat "{gsmr_out_dir}"/*.gsmr | head -n 1
    cat "{gsmr_out_dir}"/*.gsmr  | grep -vE 'Exposure|nan'
    ) | column -t > "{gsmr_out_filtered}"

    python3 _scripts/summarize_gsmr.py \
            "{gsmr_out_filtered}" \
            "{exposure_qtl}.*" \
            "{outcome_qtl}.*"\
            > "{gsmr_out_filtered}.summary"
    ''')

    plot_gsmr = inject(region, paths, r'''
    #!/usr/bin/env bash
    #SBATCH --time 8:00:00
    #SBATCH --mem 40G

    for file in "{gsmr_out_dir}*.eff_plot.gz"
    do
        base={gsmr_out_dir}$(basename $file .eff_plot.gz)
        echo $base
        "/hpc/local/CentOS7/dhl_ec/software/R-3.4.0/bin/Rscript" _scripts/run_gsmr_plot.r \
            "$file" \
            "$base.gsmr" \
            "{gsmr_plot_dir}"
    done
    ''')

    basedirs = sorted(list(set(map(os.path.dirname, paths))))

    create_dirs = '\n'.join('mkdir -p {0}'.format(basedir) for basedir in basedirs) + '\n'
    create_dirs = textwrap.dedent(r'''
    #!/usr/bin/env bash
    {0}
    '''.lstrip('\n')).format(create_dirs)

    return [
        dict(jobname='create_dirs', job=create_dirs, depends=[]),
        dict(jobname='exp', job=qtltools_exposure, depends=['create_dirs']),
        dict(jobname='expc', job=qtltools_exposure_collect, depends=['exp']),
        dict(jobname='out', job=qtltools_outcome, depends=['create_dirs']),
        dict(jobname='outc', job=qtltools_outcome_collect, depends=['out']),
        dict(jobname='bed', job=create_bed, depends=['create_dirs']),
        dict(jobname='pairs', job=gsmr_pairs, depends=['expc', 'outc', 'bed']),
        dict(jobname='gsmr', job=run_gsmr, depends=['pairs']),
        dict(jobname='summary', job=gsmr_summarize, depends=['gsmr']),
        dict(jobname='plot', job=plot_gsmr, depends=['summary']),
    ]
