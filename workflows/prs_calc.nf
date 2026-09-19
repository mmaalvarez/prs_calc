include { PREPARE_SCOREFILE  } from '../modules/local/prepare_scorefile'
include { CALCULATE_PRS      } from '../modules/local/calculate_prs'
include { PLOT_PRS; PLOT_PRS_COHORT } from '../modules/local/plot_prs'
include { MERGE_PRS_REPORTS  } from '../modules/local/merge_prs_reports'


workflow PRSCALC {

    take:
    ch_samples
    ch_scorefile
    ch_target_build
    ch_phenotypes

    main:

    /*
     * Normalise and validate the score file once.
     */
    PREPARE_SCOREFILE(
        ch_scorefile,
        ch_target_build
    )

    /*
     * Cross each sample with the one prepared score file.
     */
    ch_calculation_inputs = ch_samples.combine(
        PREPARE_SCOREFILE.out.scorefile
    )

    /*
     * One independent task is created per sample.
     */
    CALCULATE_PRS(
        ch_calculation_inputs,
        ch_target_build
    )

    /*
     * Per-sample contribution and variant-status plots.
     */
    PLOT_PRS(
        CALCULATE_PRS.out.results
    )

    /*
     * Merge all per-sample outputs after every sample has completed.
     */
    ch_summary_files = CALCULATE_PRS.out.results.map {
        meta, summary, details -> summary
    }

    ch_detail_files = CALCULATE_PRS.out.results.map {
        meta, summary, details -> details
    }

    MERGE_PRS_REPORTS(
        ch_summary_files.collect(),
        ch_detail_files.collect()
    )

    /*
     * Cohort-level plots. When ch_phenotypes is empty because
     * --phenotypes was omitted, this process does not execute.
     */
    PLOT_PRS_COHORT(
        MERGE_PRS_REPORTS.out.report,
        ch_phenotypes
    )

    emit:
    report         = MERGE_PRS_REPORTS.out.report
    variant_report = MERGE_PRS_REPORTS.out.variant_report
    plots          = PLOT_PRS.out.plots
    cohort_plots   = PLOT_PRS_COHORT.out.plots
}
