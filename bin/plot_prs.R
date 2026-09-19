#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(ragg)
})

stopf <- function(...) {
    stop(sprintf(...), call. = FALSE)
}

has_value <- function(value) {
    !is.null(value) &&
        length(value) == 1L &&
        !is.na(value) &&
        nzchar(trimws(as.character(value)))
}

format_values <- function(values, limit = 20L) {
    values <- as.character(values)
    shown <- head(values, limit)

    suffix <- if (length(values) > limit) {
        sprintf(" ... and %d more", length(values) - limit)
    } else {
        ""
    }

    paste0(
        paste(shown, collapse = ", "),
        suffix
    )
}

option_list <- list(
    make_option(
        "--summary",
        dest = "summary",
        type = "character"
    ),
    make_option(
        "--details",
        dest = "details",
        type = "character"
    ),
    make_option(
        "--sample-id",
        dest = "sample_id",
        type = "character"
    ),
    make_option(
        "--top-n",
        dest = "top_n",
        type = "integer",
        default = 25L
    ),
    make_option(
        "--contribution-output",
        dest = "contribution_output",
        type = "character"
    ),
    make_option(
        "--status-output",
        dest = "status_output",
        type = "character"
    ),
    make_option(
        "--phenotypes",
        dest = "phenotypes",
        type = "character",
        default = NULL,
        help = paste0(
            "Optional two-column TSV containing sample ID and ",
            "case/control phenotype"
        )
    ),
    make_option(
        "--roc-output",
        dest = "roc_output",
        type = "character",
        default = NULL
    ),
    make_option(
        "--or-output",
        dest = "or_output",
        type = "character",
        default = NULL
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)


write_cohort_plots <- function() {
    required_options <- c(
        "summary",
        "phenotypes",
        "roc_output",
        "or_output"
    )

    missing_options <- required_options[
        !vapply(
            required_options,
            function(option_name) {
                has_value(opt[[option_name]])
            },
            logical(1)
        )
    ]

    if (length(missing_options) > 0L) {
        option_labels <- paste0(
            "--",
            gsub(
                "_",
                "-",
                missing_options,
                fixed = TRUE
            )
        )

        stopf(
            "Missing required cohort plotting option(s): %s",
            paste(option_labels, collapse = ", ")
        )
    }

    if (!requireNamespace("pROC", quietly = TRUE)) {
        stopf(
            paste0(
                "The pROC R package is required for cohort ROC plots. ",
                "Install it with the Conda package 'r-proc'."
            )
        )
    }

    phenotype_raw <- tryCatch(
        {
            read_tsv(
                opt$phenotypes,
                col_names = FALSE,
                comment = "#",
                na = c("", "NA"),
                trim_ws = TRUE,
                col_types = cols(
                    .default = col_character()
                ),
                show_col_types = FALSE,
                progress = FALSE
            )
        },
        error = function(error) {
            stopf(
                "Could not read phenotype file '%s': %s",
                opt$phenotypes,
                conditionMessage(error)
            )
        }
    )

    if (nrow(phenotype_raw) == 0L) {
        stopf(
            "The phenotype file '%s' is empty",
            opt$phenotypes
        )
    }

    if (ncol(phenotype_raw) != 2L) {
        stopf(
            paste0(
                "The phenotype file must contain exactly two ",
                "tab-separated columns: sample ID and phenotype. ",
                "Found %d columns."
            ),
            ncol(phenotype_raw)
        )
    }

    #
    # A header is optional. The documented format is headerless, but
    # accepting a simple header avoids treating it as a sample.
    #
    first_id <- tolower(
        trimws(
            as.character(
                phenotype_raw[[1]][[1]]
            )
        )
    )

    first_phenotype <- tolower(
        trimws(
            as.character(
                phenotype_raw[[2]][[1]]
            )
        )
    )

    has_header <- (
        first_id %in% c(
            "sample",
            "sample_id",
            "id"
        ) &&
        first_phenotype %in% c(
            "phenotype",
            "status",
            "label",
            "case_control",
            "case-control"
        )
    )

    if (has_header) {
        phenotype_raw <- phenotype_raw[-1L, ]
    }

    if (nrow(phenotype_raw) == 0L) {
        stopf(
            "The phenotype file does not contain any phenotype records"
        )
    }

    phenotypes <- tibble::tibble(
        sample_id = trimws(
            as.character(
                phenotype_raw[[1]]
            )
        ),
        phenotype = tolower(
            trimws(
                as.character(
                    phenotype_raw[[2]]
                )
            )
        )
    )

    missing_phenotype_values <- (
        is.na(phenotypes$sample_id) |
        !nzchar(phenotypes$sample_id) |
        is.na(phenotypes$phenotype) |
        !nzchar(phenotypes$phenotype)
    )

    if (any(missing_phenotype_values)) {
        stopf(
            "Missing sample ID or phenotype at phenotype row(s): %s",
            paste(
                which(missing_phenotype_values),
                collapse = ", "
            )
        )
    }

    duplicated_ids <- unique(
        phenotypes$sample_id[
            duplicated(phenotypes$sample_id)
        ]
    )

    if (length(duplicated_ids) > 0L) {
        stopf(
            "Duplicate sample IDs in the phenotype file: %s",
            format_values(duplicated_ids)
        )
    }

    valid_labels <- c(
        "control",
        "case"
    )

    invalid_labels <- sort(
        unique(
            phenotypes$phenotype[
                !(phenotypes$phenotype %in% valid_labels)
            ]
        )
    )

    if (length(invalid_labels) > 0L) {
        stopf(
            paste0(
                "Unsupported phenotype label(s): %s. ",
                "Phenotypes must be 'case' or 'control'."
            ),
            paste(invalid_labels, collapse = ", ")
        )
    }

    score_table <- tryCatch(
        {
            read_tsv(
                opt$summary,
                na = c("", "NA"),
                col_types = cols(
                    .default = col_guess(),
                    sample_id = col_character(),
                    prs = col_double()
                ),
                show_col_types = FALSE,
                progress = FALSE
            )
        },
        error = function(error) {
            stopf(
                "Could not read PRS summary '%s': %s",
                opt$summary,
                conditionMessage(error)
            )
        }
    )

    required_score_columns <- c(
        "sample_id",
        "prs"
    )

    missing_score_columns <- setdiff(
        required_score_columns,
        names(score_table)
    )

    if (length(missing_score_columns) > 0L) {
        stopf(
            "PRS summary is missing required column(s): %s",
            paste(
                missing_score_columns,
                collapse = ", "
            )
        )
    }

    if (nrow(score_table) == 0L) {
        stopf("The merged PRS summary is empty")
    }

    sample_ids <- trimws(
        as.character(
            score_table$sample_id
        )
    )

    missing_sample_ids <- (
        is.na(sample_ids) |
        !nzchar(sample_ids)
    )

    if (any(missing_sample_ids)) {
        stopf(
            "Missing sample IDs in the merged PRS summary"
        )
    }

    duplicated_score_ids <- unique(
        sample_ids[
            duplicated(sample_ids)
        ]
    )

    if (length(duplicated_score_ids) > 0L) {
        stopf(
            "Duplicate sample IDs in the PRS summary: %s",
            format_values(duplicated_score_ids)
        )
    }

    prs_values <- suppressWarnings(
        as.numeric(
            score_table$prs
        )
    )

    invalid_prs <- (
        is.na(prs_values) |
        !is.finite(prs_values)
    )

    if (any(invalid_prs)) {
        stopf(
            "Missing or non-finite PRS values for sample(s): %s",
            format_values(
                sample_ids[invalid_prs]
            )
        )
    }

    scores <- tibble::tibble(
        sample_id = sample_ids,
        prs = prs_values
    )

    samples_without_phenotypes <- setdiff(
        scores$sample_id,
        phenotypes$sample_id
    )

    if (length(samples_without_phenotypes) > 0L) {
        stopf(
            paste0(
                "No phenotype was provided for PRS sample(s): %s. ",
                "Phenotype IDs must match the pipeline sample_id values."
            ),
            format_values(samples_without_phenotypes)
        )
    }

    unused_phenotypes <- setdiff(
        phenotypes$sample_id,
        scores$sample_id
    )

    if (length(unused_phenotypes) > 0L) {
        warning(
            sprintf(
                "Ignoring phenotype IDs with no PRS result: %s",
                format_values(unused_phenotypes)
            ),
            call. = FALSE
        )
    }

    cohort <- scores %>%
        left_join(
            phenotypes,
            by = "sample_id"
        )

    n_samples <- nrow(cohort)
    n_cases <- sum(
        cohort$phenotype == "case"
    )
    n_controls <- sum(
        cohort$phenotype == "control"
    )

    if (n_cases == 0L || n_controls == 0L) {
        stopf(
            paste0(
                "ROC and odds-ratio plots require both cases and ",
                "controls. Found %d cases and %d controls."
            ),
            n_cases,
            n_controls
        )
    }

    if (n_samples < 10L) {
        stopf(
            paste0(
                "At least 10 samples are required to form ten PRS ",
                "deciles. Found %d samples."
            ),
            n_samples
        )
    }

    if (length(unique(cohort$prs)) < 2L) {
        warning(
            paste0(
                "All samples have the same PRS. The ROC curve and ",
                "decile assignments will not be informative."
            ),
            call. = FALSE
        )
    }

    #
    # Force the biologically expected direction: increasing PRS is
    # interpreted as increasing probability of being a case. This
    # allows an AUC below 0.5 rather than automatically reversing it.
    #
    roc_object <- tryCatch(
        {
            pROC::roc(
                response = factor(
                    cohort$phenotype,
                    levels = c(
                        "control",
                        "case"
                    )
                ),
                predictor = cohort$prs,
                levels = c(
                    "control",
                    "case"
                ),
                direction = "<",
                quiet = TRUE
            )
        },
        error = function(error) {
            stopf(
                "Could not calculate the ROC curve: %s",
                conditionMessage(error)
            )
        }
    )

    auc_value <- as.numeric(
        pROC::auc(roc_object)
    )

    roc_data <- tibble::tibble(
        false_positive_rate = 1 -
            as.numeric(
                roc_object$specificities
            ),
        sensitivity = as.numeric(
            roc_object$sensitivities
        )
    ) %>%
        arrange(
            false_positive_rate,
            sensitivity
        ) %>%
        distinct(
            false_positive_rate,
            sensitivity,
            .keep_all = TRUE
        )

    roc_plot <- ggplot(
        roc_data,
        aes(
            x = false_positive_rate,
            y = sensitivity
        )
    ) +
        geom_abline(
            intercept = 0,
            slope = 1,
            linetype = "dashed",
            colour = "grey55",
            linewidth = 0.7
        ) +
        geom_step(
            direction = "vh",
            colour = "#2166AC",
            linewidth = 1.2
        ) +
        scale_x_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.2),
            expand = expansion(mult = 0)
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.2),
            expand = expansion(mult = 0)
        ) +
        coord_equal() +
        labs(
            title = "Polygenic risk score ROC curve",
            subtitle = sprintf(
                "AUC = %.3f; n = %d (%d cases, %d controls)",
                auc_value,
                n_samples,
                n_cases,
                n_controls
            ),
            x = "False-positive rate (1 - specificity)",
            y = "Sensitivity"
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title.position = "plot"
        )

    ggsave(
        filename = opt$roc_output,
        plot = roc_plot,
        device = ragg::agg_png,
        width = 12.5,
        height = 7,
        units = "in",
        dpi = 300,
        bg = "white"
    )

    #
    # Sort by PRS and then sample ID so that tied PRS values are split
    # deterministically between deciles.
    #
    cohort <- cohort %>%
        arrange(
            prs,
            sample_id
        )

    cohort$decile <- dplyr::ntile(
        seq_len(nrow(cohort)),
        10L
    )

    or_table <- cohort %>%
        group_by(decile) %>%
        summarise(
            cases = sum(
                phenotype == "case"
            ),
            controls = sum(
                phenotype == "control"
            ),
            n_samples = n(),
            .groups = "drop"
        ) %>%
        arrange(decile)

    expected_deciles <- seq_len(10L)

    if (!all(expected_deciles %in% or_table$decile)) {
        stopf(
            "Could not assign samples to all ten PRS deciles"
        )
    }

    reference_decile <- 1L

    reference_cases <- or_table$cases[
        or_table$decile == reference_decile
    ][[1]]

    reference_controls <- or_table$controls[
        or_table$decile == reference_decile
    ][[1]]

    or_table$odds_ratio <- NA_real_
    or_table$lower <- NA_real_
    or_table$upper <- NA_real_
    or_table$continuity_corrected <- FALSE

    for (i in seq_len(nrow(or_table))) {
        decile <- or_table$decile[[i]]

        if (decile == reference_decile) {
            or_table$odds_ratio[[i]] <- 1
            or_table$lower[[i]] <- 1
            or_table$upper[[i]] <- 1
            next
        }

        cell_counts <- c(
            or_table$cases[[i]],
            or_table$controls[[i]],
            reference_cases,
            reference_controls
        )

        #
        # Use a Haldane-Anscombe correction for comparisons containing
        # a zero cell, avoiding infinite ORs and confidence intervals.
        #
        correction <- if (any(cell_counts == 0L)) {
            0.5
        } else {
            0
        }

        or_table$continuity_corrected[[i]] <- (
            correction > 0
        )

        cases <- or_table$cases[[i]] +
            correction

        controls <- or_table$controls[[i]] +
            correction

        ref_cases <- reference_cases +
            correction

        ref_controls <- reference_controls +
            correction

        log_or <- (
            log(cases) -
            log(controls) -
            log(ref_cases) +
            log(ref_controls)
        )

        standard_error <- sqrt(
            1 / cases +
            1 / controls +
            1 / ref_cases +
            1 / ref_controls
        )

        or_table$odds_ratio[[i]] <- exp(log_or)

        or_table$lower[[i]] <- exp(
            log_or -
            1.96 * standard_error
        )

        or_table$upper[[i]] <- exp(
            log_or +
            1.96 * standard_error
        )
    }

    or_table$decile_factor <- factor(
        or_table$decile,
        levels = expected_deciles
    )

    decile_labels <- sprintf(
        "%d\n(%d/%d)",
        or_table$decile,
        or_table$cases,
        or_table$controls
    )

    names(decile_labels) <- as.character(
        or_table$decile
    )

    correction_note <- if (
        any(or_table$continuity_corrected)
    ) {
        paste0(
            "A 0.5 continuity correction was used for ",
            "comparisons containing a zero cell."
        )
    } else {
        "No zero-cell continuity correction was needed."
    }

    or_plot <- ggplot(
        or_table,
        aes(
            x = decile_factor,
            y = odds_ratio
        )
    ) +
        geom_hline(
            yintercept = 1,
            linetype = "dashed",
            colour = "#B2182B",
            linewidth = 0.7
        ) +
        geom_errorbar(
            aes(
                ymin = lower,
                ymax = upper
            ),
            width = 0.22,
            colour = "grey30",
            linewidth = 0.7
        ) +
        geom_point(
            size = 3.2,
            colour = "#2166AC"
        ) +
        scale_x_discrete(
            labels = decile_labels
        ) +
        scale_y_log10() +
        labs(
            title = "PRS odds ratios by decile",
            subtitle = sprintf(
                paste0(
                    "Reference = decile %d; %d cases and ",
                    "%d controls. %s"
                ),
                reference_decile,
                n_cases,
                n_controls,
                correction_note
            ),
            x = "PRS decile (cases/controls)",
            y = paste0(
                "Odds ratio relative to decile ",
                reference_decile,
                " (log scale)"
            ),
            caption = paste0(
                "Decile 1 contains the lowest PRS values and decile 10 ",
                "the highest. Error bars are 95% Wald confidence intervals."
            )
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title.position = "plot",
            axis.text.x = element_text(
                size = 10
            )
        )

    ggsave(
        filename = opt$or_output,
        plot = or_plot,
        device = ragg::agg_png,
        width = 12.5,
        height = 7,
        units = "in",
        dpi = 300,
        bg = "white"
    )
}


if (has_value(opt$phenotypes)) {
    write_cohort_plots()

    quit(
        save = "no",
        status = 0L
    )
}


#
# Per-sample plotting mode.
#
required_options <- c(
    "summary",
    "details",
    "sample_id",
    "contribution_output",
    "status_output"
)

for (option_name in required_options) {
    value <- opt[[option_name]]

    if (!has_value(value)) {
        stopf(
            "Missing required option: --%s",
            gsub(
                "_",
                "-",
                option_name,
                fixed = TRUE
            )
        )
    }
}

top_n <- as.integer(opt$top_n)

if (is.na(top_n) || top_n < 1L) {
    stopf("--top-n must be a positive integer")
}

summary_table <- read_tsv(
    opt$summary,
    show_col_types = FALSE,
    progress = FALSE
)

details <- read_tsv(
    opt$details,
    show_col_types = FALSE,
    progress = FALSE
)

if (nrow(summary_table) != 1L) {
    stopf(
        "Expected one summary row, found %d",
        nrow(summary_table)
    )
}

if (!"prs" %in% names(summary_table)) {
    stopf(
        "The summary table is missing the 'prs' column"
    )
}

required_detail_columns <- c(
    "score_row",
    "variant_id",
    "chrom",
    "position",
    "contribution",
    "status"
)

missing_detail_columns <- setdiff(
    required_detail_columns,
    names(details)
)

if (length(missing_detail_columns) > 0L) {
    stopf(
        "The details table is missing required column(s): %s",
        paste(
            missing_detail_columns,
            collapse = ", "
        )
    )
}

if (nrow(details) == 0L) {
    stopf("The variant details table is empty")
}

prs_value <- summary_table$prs[[1]]

prs_label <- if (is.na(prs_value)) {
    "NA"
} else {
    sprintf("%.6f", prs_value)
}

plot_data <- details %>%
    mutate(
        absolute_contribution = ifelse(
            is.na(contribution),
            -Inf,
            abs(contribution)
        ),
        plot_contribution = ifelse(
            is.na(contribution),
            0,
            contribution
        ),
        direction = case_when(
            is.na(contribution) ~ "Unscored",
            contribution > 0 ~ "Positive",
            contribution < 0 ~ "Negative",
            TRUE ~ "Zero"
        ),
        label = paste0(
            variant_id,
            " (",
            chrom,
            ":",
            position,
            ")"
        )
    ) %>%
    arrange(
        desc(absolute_contribution),
        score_row
    ) %>%
    slice_head(
        n = min(
            top_n,
            nrow(details)
        )
    ) %>%
    arrange(
        absolute_contribution,
        score_row
    )

plot_data$label <- make.unique(
    plot_data$label
)

plot_data$label <- factor(
    plot_data$label,
    levels = plot_data$label
)

contribution_plot <- ggplot(
    plot_data,
    aes(
        x = plot_contribution,
        y = label,
        fill = direction
    )
) +
    geom_col(width = 0.75) +
    geom_vline(
        xintercept = 0,
        colour = "grey30",
        linewidth = 0.4
    ) +
    scale_fill_manual(
        values = c(
            Positive = "#B2182B",
            Negative = "#2166AC",
            Zero = "#999999",
            Unscored = "#F4A582"
        ),
        drop = FALSE
    ) +
    labs(
        title = paste0(
            opt$sample_id,
            ": top PRS variant contributions"
        ),
        subtitle = paste0(
            "PRS = ",
            prs_label,
            "; displaying up to ",
            top_n,
            " variants"
        ),
        x = "Effect weight \u00d7 effect-allele dosage",
        y = NULL,
        fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "bottom",
        plot.title.position = "plot"
    )

contribution_height <- max(
    5,
    min(
        14,
        2.5 + 0.28 * nrow(plot_data)
    )
)

ggsave(
    filename = opt$contribution_output,
    plot = contribution_plot,
    device = ragg::agg_png,
    width = 11,
    height = contribution_height,
    units = "in",
    dpi = 180,
    bg = "white"
)

status_data <- details %>%
    count(
        status,
        name = "n_variants"
    ) %>%
    mutate(
        status_label = ifelse(
            is.na(status),
            "missing status",
            gsub(
                "_",
                " ",
                status,
                fixed = TRUE
            )
        )
    ) %>%
    arrange(n_variants)

status_data$status_label <- factor(
    status_data$status_label,
    levels = status_data$status_label
)

status_plot <- ggplot(
    status_data,
    aes(
        x = n_variants,
        y = status_label,
        fill = status_label
    )
) +
    geom_col(width = 0.75) +
    geom_text(
        aes(label = n_variants),
        hjust = -0.2,
        size = 4
    ) +
    scale_x_continuous(
        expand = expansion(
            mult = c(0, 0.15)
        )
    ) +
    guides(fill = "none") +
    labs(
        title = paste0(
            opt$sample_id,
            ": score-variant QC"
        ),
        subtitle = paste0(
            "PRS = ",
            prs_label
        ),
        x = "Number of score variants",
        y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title.position = "plot"
    )

status_height <- max(
    4.5,
    2.5 + 0.5 * nrow(status_data)
)

ggsave(
    filename = opt$status_output,
    plot = status_plot,
    device = ragg::agg_png,
    width = 10,
    height = status_height,
    units = "in",
    dpi = 180,
    bg = "white"
)
