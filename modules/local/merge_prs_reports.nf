process MERGE_PRS_REPORTS {

    tag 'all_samples'
    
    label 'process_medium'
    label 'prs_calc_env'

    input:
    path summaries, stageAs: 'summaries/*'
    path details,   stageAs: 'details/*'

    output:
    path 'prs_scores.tsv',          emit: report
    path 'prs_variant_details.tsv', emit: variant_report

    script:
    """
    merge_prs_reports.R \
        --summary_dir "summaries" \
        --details_dir "details" \
        --summary_output "prs_scores.tsv" \
        --details_output "prs_variant_details.tsv"
    """
}
