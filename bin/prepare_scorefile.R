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
        "--target_build",
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
        "--qc_output",
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
    stopf("--target_build is required")
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
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    progress = FALSE
)

if (nrow(score) == 0L) {
    stopf("The score file does not contain any variants")
}

n_input_rows <- nrow(score)
original_score_rows <- seq_len(n_input_rows)

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
    x <- toupper(
        gsub(
            "[[:space:]]+",
            "",
            as.character(x)
        )
    )
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

present_chrom_rows <- which(!is_missing_text(chromosome))

if (length(present_chrom_rows) > 0L) {
    has_chr_prefix <- grepl(
        "^chr",
        chromosome[present_chrom_rows],
        ignore.case = TRUE
    )

    if (any(has_chr_prefix) && any(!has_chr_prefix)) {
        warning(
            paste0(
                "Score file mixes chromosome names with and without ",
                "'chr'; equivalent names will be treated as the ",
                "same chromosome."
            ),
            call. = FALSE
        )
    }
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

invalid_inferred_other <- rep(FALSE, n_input_rows)

if ("hm_inferOtherAllele" %in% names(score)) {
    inferred_other <- clean_allele(score$hm_inferOtherAllele)

# Under this policy a slash-separated inferred value is a set of
    # allowed non-effect alleles, not one literal allele.
    infer_rows <- which(
        is.na(other_allele) & !is.na(inferred_other)
    )

    if (length(infer_rows) > 0L) {
        valid_inferred <- grepl(
            "^[ACGT](/[ACGT])*$",
            inferred_other[infer_rows]
        )

        invalid_inferred_other[
            infer_rows[!valid_inferred]
        ] <- TRUE

        other_allele[
            infer_rows[valid_inferred]
        ] <- inferred_other[
            infer_rows[valid_inferred]
        ]
    }
}

# Keep only single-base A/C/G/T score variants. The missing
# other_allele allowed by this schema is not itself a non-SNV.
bad_effect_allele <- !grepl(
    "^[ACGT](/[ACGT])*$",
    effect_allele
)

bad_other_allele <- !is.na(other_allele) &
    !grepl(
        "^[ACGT](/[ACGT])*$",
        other_allele
    )

skip_non_snv <- bad_effect_allele |
    bad_other_allele |
    invalid_inferred_other

skipped_non_snv_rows <- original_score_rows[skip_non_snv]
n_skipped_non_snv <- length(skipped_non_snv_rows)

if (n_skipped_non_snv > 0L) {
    warning(
        sprintf(
            paste0(
                "Skipped %d non-SNV/unsupported-allele score row(s); ",
                "first original row number(s): %s"
            ),
            n_skipped_non_snv,
            paste(head(skipped_non_snv_rows, 20L), collapse = ", ")
        ),
        call. = FALSE
    )
}

keep <- !skip_non_snv

if (!any(keep)) {
    stopf(
        "No SNV score variants remain after skipping %d unsupported row(s)",
        n_skipped_non_snv
    )
}

score <- score[keep, , drop = FALSE]
chromosome <- chromosome[keep]
position <- position[keep]
effect_allele <- effect_allele[keep]
other_allele <- other_allele[keep]
original_score_rows <- original_score_rows[keep]

normalise_allele_set <- function(value) {
    if (is.na(value)) {
        return(NA_character_)
    }

    paste(
        sort(
            unique(
                strsplit(
                    value,
                    "/",
                    fixed = TRUE
                )[[1]]
            )
        ),
        collapse = "/"
    )
}

effect_allele <- vapply(
    effect_allele,
    normalise_allele_set,
    character(1),
    USE.NAMES = FALSE
)

other_allele <- vapply(
    other_allele,
    normalise_allele_set,
    character(1),
    USE.NAMES = FALSE
)

overlapping_alleles <- vapply(
    seq_along(effect_allele),
    function(i) {
        if (is.na(other_allele[[i]])) {
            return(FALSE)
        }

        effect_set <- strsplit(
            effect_allele[[i]],
            "/",
            fixed = TRUE
        )[[1]]

        other_set <- strsplit(
            other_allele[[i]],
            "/",
            fixed = TRUE
        )[[1]]

        length(intersect(effect_set, other_set)) > 0L
    },
    logical(1)
)

if (any(overlapping_alleles)) {
    stopf(
        "Effect and other allele sets overlap at original score row(s): %s",
        paste(
            head(
                original_score_rows[overlapping_alleles],
                20L
            ),
            collapse = ", "
        )
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
    score_row = original_score_rows,
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
    n_input_rows = n_input_rows,
    n_output_rows = nrow(normalized),
    n_skipped_non_snv = n_skipped_non_snv,
    skipped_non_snv_score_rows_first_20 = if (
        n_skipped_non_snv > 0L
    ) {
        paste(
            head(skipped_non_snv_rows, 20L),
            collapse = ","
        )
    } else {
        NA_character_
    },
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
