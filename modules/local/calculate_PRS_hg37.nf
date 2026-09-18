process calculate_PRS_hg37 {

    label 'short_higher'
    
    conda '/home/malvarez/.conda/envs/R/'

    input:
    val(sampleId)
    path(metadata_table)
    path(prs_model)

    output:
    path('predictions.tsv'), emit: predictions

    script:
    """
    calculate_PRS_hg37.R ${sampleId} ${metadata_table} ${prs_model}
    """
}
