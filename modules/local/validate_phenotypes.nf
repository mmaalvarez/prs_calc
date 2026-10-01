process VALIDATE_PHENOTYPES {

    label 'process_low'
    label 'prs_calc_env'

    input:
    path phenotypes, stageAs: 'phenotypes.tsv'
    val sample_ids

    output:
    path 'checked_phenotypes.tsv', emit: checked

    script:
    def quotedIds = sample_ids.collect { id ->
        "'" + id.toString().replace("'", "'\"'\"'") + "'"
    }.join(' ')

    """
    printf '%s\\n' ${quotedIds} > sample_ids.txt

    validate_phenotypes.py \
        --phenotypes phenotypes.tsv \
        --sample-ids sample_ids.txt \
        --output checked_phenotypes.tsv
    """
}
