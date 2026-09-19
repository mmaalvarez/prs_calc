process PLOT_PRS {

    tag "${meta.id}"

    label 'process_low'

    conda "${projectDir}/envs/prs_calc.yml"

    input:
    tuple val(meta), path(summary), path(details)

    output:
    tuple val(meta),
          path("${meta.id}.contributions.png"),
          path("${meta.id}.variant_status.png"),
          emit: plots

    script:
    def topVariants = params.plot_top_variants as Integer

    """
    Rscript "${projectDir}/bin/plot_prs.R" \
        --summary "${summary}" \
        --details "${details}" \
        --sample-id "${meta.id}" \
        --top-n "${topVariants}" \
        --contribution-output "${meta.id}.contributions.png" \
        --status-output "${meta.id}.variant_status.png"
    """
}


process PLOT_PRS_COHORT {

    tag 'all_samples'

    label 'process_low'

    conda "${projectDir}/envs/prs_calc.yml"

    input:
    path summary,    stageAs: 'cohort_prs_scores.tsv'
    path phenotypes, stageAs: 'cohort_phenotypes.tsv'

    output:
    tuple path('prs_roc.png'),
          path('prs_or_deciles.png'),
          emit: plots

    script:
    """
    Rscript "${projectDir}/bin/plot_prs.R" \
        --summary "${summary}" \
        --phenotypes "${phenotypes}" \
        --roc-output "prs_roc.png" \
        --or-output "prs_or_deciles.png"
    """
}
