#!/usr/bin/env nextflow

nextflow.enable.dsl=2

println "Project directory: $projectDir"
println "Launch directory: $launchDir"
println "Working directory: $workDir"

include { calculate_PRS_hg37 } from './modules/calculate_PRS_hg37'


workflow {

    // channel by vcf chunks

    sample_IDs_list = Channel
        .fromPath(params.metadata_table)
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.IID } // Extracts only the value under the "IID" column

    calculate_PRS_hg37(
        sample_IDs_list,
        file(params.metadata_table),
        file(params.prs_model)
    )

    predictions = calculate_PRS_hg37.out.predictions
        .flatten()
        .collectFile(name: 'predictions', keepHeader: true,
                     storeDir: 'res/')
}


workflow.onError {
    println "Pipeline execution stopped with error: ${workflow.errorMessage}"
}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    println "Duration: $workflow.duration"
}
