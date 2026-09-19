process PREPARE_SCOREFILE {

    tag "${scorefile.name}"

    label 'process_low'

    conda "${projectDir}/envs/prs_calc.yml"

    input:
    path scorefile
    val target_build

    output:
    path 'normalized_scorefile.tsv', emit: scorefile
    path 'scorefile_qc.tsv',         emit: qc

    script:
    """
    Rscript "${projectDir}/bin/prepare_scorefile.R" \
        --scorefile "${scorefile}" \
        --target-build "${target_build}" \
        --output "normalized_scorefile.tsv" \
        --qc-output "scorefile_qc.tsv"
    """
}
