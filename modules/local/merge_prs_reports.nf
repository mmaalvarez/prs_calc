process MERGE_PRS_REPORTS {

    tag 'all_samples'

    label 'process_low'

    conda "${projectDir}/envs/prs_calc.yml"

    input:
    path summaries, stageAs: 'summaries/*'
    path details,   stageAs: 'details/*'

    output:
    path 'prs_scores.tsv',          emit: report
    path 'prs_variant_details.tsv', emit: variant_report

    script:
    """
    Rscript "${projectDir}/bin/merge_prs_reports.R" \
        --summary-dir "summaries" \
        --details-dir "details" \
        --summary-output "prs_scores.tsv" \
        --details-output "prs_variant_details.tsv"
    """
}
