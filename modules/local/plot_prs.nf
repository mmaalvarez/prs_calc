process PLOT_PRS {

    tag "${meta.id}"

    label 'process_low'
    label 'prs_calc_env'

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
    plot_prs.R \
        --summary "${summary}" \
        --details "${details}" \
        --sample_id "${meta.id}" \
        --top_n "${topVariants}" \
        --contribution_output "${meta.id}.contributions.png" \
        --status_output "${meta.id}.variant_status.png"
    """
}


process PLOT_PRS_COHORT {

    tag 'all_samples'

    label 'process_low'
    label 'prs_calc_env'

    input:
    path summary,    stageAs: 'cohort_prs_scores.tsv'
    path phenotypes, stageAs: 'cohort_phenotypes.tsv'

    output:
    tuple path('prs_roc.png'),
          path('prs_or_deciles.png'),
          emit: plots

    script:
    """
    plot_prs.R \
        --summary "${summary}" \
        --phenotypes "${phenotypes}" \
        --roc_output "prs_roc.png" \
        --or_output "prs_or_deciles.png"
    """
}
