#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(readr)
    library(dplyr)
})

stopf <- function(...) {
    stop(sprintf(...), call. = FALSE)
}

option_list <- list(
    make_option(
        "--summary-dir",
        dest = "summary_dir",
        type = "character",
        default = "summaries"
    ),
    make_option(
        "--details-dir",
        dest = "details_dir",
        type = "character",
        default = "details"
    ),
    make_option(
        "--summary-output",
        dest = "summary_output",
        type = "character",
        default = "prs_scores.tsv"
    ),
    make_option(
        "--details-output",
        dest = "details_output",
        type = "character",
        default = "prs_variant_details.tsv"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)

summary_files <- list.files(
    opt$summary_dir,
    pattern = "\\.prs\\.tsv$",
    full.names = TRUE
)

detail_files <- list.files(
    opt$details_dir,
    pattern = "\\.variants\\.tsv$",
    full.names = TRUE
)

if (length(summary_files) == 0L) {
    stopf(
        "No per-sample PRS summary files were found in '%s'",
        opt$summary_dir
    )
}

if (length(detail_files) == 0L) {
    stopf(
        "No per-sample variant detail files were found in '%s'",
        opt$details_dir
    )
}

summary_tables <- lapply(
    summary_files,
    function(filename) {
        table <- read_tsv(
            filename,
            show_col_types = FALSE,
            progress = FALSE
        )

        if (nrow(table) != 1L) {
            stopf(
                "Expected one row in '%s', found %d",
                basename(filename),
                nrow(table)
            )
        }

        table
    }
)

combined_summary <- bind_rows(
    summary_tables
)

duplicated_samples <- unique(
    combined_summary$sample_id[
        duplicated(combined_summary$sample_id)
    ]
)

if (length(duplicated_samples) > 0L) {
    stopf(
        "Duplicate sample IDs in summary files: %s",
        paste(
            duplicated_samples,
            collapse = ", "
        )
    )
}

combined_summary <- combined_summary %>%
    arrange(sample_id)

detail_tables <- lapply(
    detail_files,
    function(filename) {
        read_tsv(
            filename,
            show_col_types = FALSE,
            progress = FALSE
        )
    }
)

combined_details <- bind_rows(
    detail_tables
) %>%
    arrange(
        sample_id,
        score_row
    )

write_tsv(
    combined_summary,
    opt$summary_output,
    na = "NA"
)

write_tsv(
    combined_details,
    opt$details_output,
    na = "NA"
)
