#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(readr)
    library(tibble)
})

stopf <- function(...) {
    stop(sprintf(...), call. = FALSE)
}

option_list <- list(
    make_option(
        "--scorefile",
        dest = "scorefile",
        type = "character"
    ),
    make_option(
        "--target-build",
        dest = "target_build",
        type = "character"
    ),
    make_option(
        "--output",
        dest = "output",
        type = "character",
        default = "normalized_scorefile.tsv"
    ),
    make_option(
        "--qc-output",
        dest = "qc_output",
        type = "character",
        default = "scorefile_qc.tsv"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)

if (is.null(opt$scorefile) || !nzchar(opt$scorefile)) {
    stopf("--scorefile is required")
}

if (is.null(opt$target_build) || !nzchar(opt$target_build)) {
    stopf("--target-build is required")
}

target_build <- tolower(trimws(opt$target_build))

if (!target_build %in% c("hg38", "hg37", "hg19")) {
    stopf(
        "Unsupported target build '%s'. Expected hg38, hg37 or hg19.",
        target_build
    )
}

position_column <- paste0("chr_position_", target_build)

score <- read_tsv(
    opt$scorefile,
    comment = "#",
    na = c("", "NA", "."),
    trim_ws = TRUE,
    name_repair = "check_unique",
    show_col_types = FALSE,
    progress = FALSE
)

if (nrow(score) == 0L) {
    stopf("The score file does not contain any variants")
}

required_columns <- c(
    "effect_allele",
    "effect_weight",
    position_column
)

missing_columns <- setdiff(required_columns, names(score))

if (length(missing_columns) > 0L) {
    stopf(
        "Score file is missing required column(s): %s. Available columns: %s",
        paste(missing_columns, collapse = ", "),
        paste(names(score), collapse = ", ")
    )
}

if (!any(c("chr_name", "hm_chr") %in% names(score))) {
    stopf(
        "The score file must contain either 'chr_name' or 'hm_chr'"
    )
}

is_missing_text <- function(x) {
    is.na(x) | trimws(as.character(x)) %in% c("", ".", "NA")
}

normalise_chromosome <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^chr", "", x, ignore.case = TRUE)
    x <- toupper(x)

    x[x %in% c("23")] <- "X"
    x[x %in% c("24")] <- "Y"
    x[x %in% c("M", "MTDNA", "25")] <- "MT"

    x
}

clean_allele <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x[x %in% c("", ".", "NA")] <- NA_character_
    x
}

chromosome <- rep(NA_character_, nrow(score))

if ("chr_name" %in% names(score)) {
    chromosome <- as.character(score$chr_name)
}

if ("hm_chr" %in% names(score)) {
    use_hm_chr <- is_missing_text(chromosome)
    chromosome[use_hm_chr] <- as.character(score$hm_chr[use_hm_chr])
}

chromosome <- normalise_chromosome(chromosome)

position_raw <- score[[position_column]]
position_numeric <- suppressWarnings(as.numeric(position_raw))

bad_position <- (
    is.na(position_numeric) |
    !is.finite(position_numeric) |
    position_numeric < 1 |
    abs(position_numeric - round(position_numeric)) > 0
)

if (any(bad_position)) {
    bad_rows <- which(bad_position)

    stopf(
        "Invalid or missing values in '%s' at row(s): %s",
        position_column,
        paste(head(bad_rows, 20L), collapse = ", ")
    )
}

position <- as.integer(round(position_numeric))

if (any(is_missing_text(chromosome))) {
    bad_rows <- which(is_missing_text(chromosome))

    stopf(
        "Missing chromosome values at row(s): %s",
        paste(head(bad_rows, 20L), collapse = ", ")
    )
}

effect_allele <- clean_allele(score$effect_allele)

if (any(is.na(effect_allele))) {
    bad_rows <- which(is.na(effect_allele))

    stopf(
        "Missing effect alleles at row(s): %s",
        paste(head(bad_rows, 20L), collapse = ", ")
    )
}

other_allele <- rep(NA_character_, nrow(score))

if ("other_allele" %in% names(score)) {
    other_allele <- clean_allele(score$other_allele)
}

if ("hm_inferOtherAllele" %in% names(score)) {
    inferred_other <- clean_allele(score$hm_inferOtherAllele)

    # PGS Catalog harmonized files can contain slash-separated candidate
    # alleles at multiallelic SNV positions, for example "A/T" or
    # "A/C/T". These are sets of possible other alleles, not one literal
    # allele. The normalized schema has a scalar other_allele field, so
    # only copy an inferred value when it contains exactly one allele.
    # Ambiguous candidate sets remain NA rather than selecting an allele
    # arbitrarily or duplicating the score row.
    infer_rows <- which(
        is.na(other_allele) & !is.na(inferred_other)
    )

    if (length(infer_rows) > 0L) {
        valid_inferred <- grepl(
            "^[ACGT](/[ACGT])*$",
            inferred_other[infer_rows]
        )

        if (any(!valid_inferred)) {
            bad_rows <- infer_rows[!valid_inferred]

            stopf(
                paste0(
                    "Unsupported value(s) in 'hm_inferOtherAllele' ",
                    "at row(s): %s. Expected a single A/C/G/T allele ",
                    "or a slash-separated list of A/C/G/T alleles."
                ),
                paste(head(bad_rows, 20L), collapse = ", ")
            )
        }

        single_inferred_rows <- infer_rows[
            grepl(
                "^[ACGT]$",
                inferred_other[infer_rows]
            )
        ]

        other_allele[single_inferred_rows] <- inferred_other[single_inferred_rows]
    }
}

# This implementation intentionally supports single-nucleotide variants.
# Rejecting indels is safer than silently applying the wrong reference
# allele to an absent indel.

bad_effect_allele <- !grepl("^[ACGT]$", effect_allele)
bad_other_allele <- !is.na(other_allele) &
    !grepl("^[ACGT]$", other_allele)

if (any(bad_effect_allele | bad_other_allele)) {
    bad_rows <- which(bad_effect_allele | bad_other_allele)

    stopf(
        paste0(
            "Only single-nucleotide A/C/G/T alleles are supported. ",
            "Non-SNV allele(s) were found at row(s): %s"
        ),
        paste(head(bad_rows, 20L), collapse = ", ")
    )
}

same_allele <- !is.na(other_allele) &
    effect_allele == other_allele

if (any(same_allele)) {
    bad_rows <- which(same_allele)

    stopf(
        "Effect and other allele are identical at row(s): %s",
        paste(head(bad_rows, 20L), collapse = ", ")
    )
}

effect_weight <- suppressWarnings(
    as.numeric(score$effect_weight)
)

if (
    any(is.na(effect_weight)) ||
    any(!is.finite(effect_weight))
) {
    bad_rows <- which(
        is.na(effect_weight) |
        !is.finite(effect_weight)
    )

    stopf(
        "Invalid effect weights at row(s): %s",
        paste(head(bad_rows, 20L), collapse = ", ")
    )
}

variant_id <- rep(NA_character_, nrow(score))

for (candidate in c("rsID", "hm_rsID", "variant_id")) {
    if (candidate %in% names(score)) {
        candidate_values <- trimws(as.character(score[[candidate]]))
        candidate_values[
            candidate_values %in% c("", ".", "NA")
        ] <- NA_character_

        fill_rows <- is.na(variant_id) &
            !is.na(candidate_values)

        variant_id[fill_rows] <- candidate_values[fill_rows]
    }
}

fallback_id <- paste0(
    chromosome,
    ":",
    position,
    ":",
    effect_allele,
    ":",
    ifelse(is.na(other_allele), "NA", other_allele)
)

variant_id[is.na(variant_id)] <- fallback_id[is.na(variant_id)]

normalized <- tibble(
    score_row = seq_len(nrow(score)),
    variant_id = variant_id,
    chrom = chromosome,
    position = position,
    effect_allele = effect_allele,
    other_allele = other_allele,
    effect_weight = effect_weight
)

coordinate <- paste(
    normalized$chrom,
    normalized$position,
    sep = ":"
)

coordinate_counts <- table(coordinate)

qc <- tibble(
    source_scorefile = basename(opt$scorefile),
    target_build = target_build,
    position_column = position_column,
    n_input_rows = nrow(score),
    n_output_rows = nrow(normalized),
    n_unique_positions = length(unique(coordinate)),
    n_positions_with_multiple_score_rows = sum(
        coordinate_counts > 1L
    ),
    n_missing_other_allele = sum(
        is.na(normalized$other_allele)
    )
)

write_tsv(
    normalized,
    opt$output,
    na = "NA"
)

write_tsv(
    qc,
    opt$qc_output,
    na = "NA"
)
