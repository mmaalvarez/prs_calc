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
    def q = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'"}
    def contributionFile = "${meta.id}.contributions.png"
    def statusFile = "${meta.id}.variant_status.png"
    """
    plot_prs.R \
        --summary ${q.call(summary)} \
        --details ${q.call(details)} \
        --sample_id ${q.call(meta.id)} \
        --top_n "${topVariants}" \
        --contribution_output ${q.call(contributionFile)} \
        --status_output ${q.call(statusFile)}
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
        --summary cohort_prs_scores.tsv \
        --phenotypes cohort_phenotypes.tsv \
        --roc_output prs_roc.png \
        --or_output prs_or_deciles.png
    """
}
