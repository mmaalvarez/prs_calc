process PREPARE_SCOREFILE {

    tag "${scorefile.name}"

    label 'process_low'
    label 'prs_calc_env'

    input:
    path scorefile
    val target_build

    output:
    path 'normalized_scorefile.tsv', emit: scorefile
    path 'scorefile_qc.tsv',         emit: qc

    script:
    def q = { value ->
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }

    def nonAdditive = params.non_additive
        .toString().trim().toLowerCase()

    """
    prepare_scorefile.R \
        --scorefile ${q.call(scorefile)} \
        --target_build ${q.call(target_build)} \
        --non_additive ${q.call(nonAdditive)} \
        --output normalized_scorefile.tsv \
        --qc_output scorefile_qc.tsv
    """
}
