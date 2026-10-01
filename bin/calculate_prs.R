#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(readr)
    library(dplyr)
    library(tibble)
    library(Biostrings)
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
        "--genotypes",
        dest = "genotypes",
        type = "character"
    ),
    make_option(
        "--sample_id",
        dest = "sample_id",
        type = "character"
    ),
    make_option(
        "--vcf_sample_file",
        dest = "vcf_sample_file",
        type = "character"
    ),
    make_option(
        "--source_vcf",
        dest = "source_vcf",
        type = "character"
    ),
    make_option(
        "--target_build",
        dest = "target_build",
        type = "character"
    ),
    make_option(
        "--missing_genotype",
        dest = "missing_genotype",
        type = "character",
        default = NA_character_,
        help = paste0(
            "Always required. How to score an absent SNP (i.e. no row or hom-ref block covering it): ",
            "'reference' (assume homozygous for the reference-genome base; VCF-only), ",
            "'zero' (dosage 0), or 'error' (abort)."
        )
    ),
    make_option(
        "--no_calls",
        dest = "no_calls",
        type = "character",
        default = NA_character_,
        help = paste0(
            "Always required. For an explicit score-SNP GT=./.: ",
            "'zero' assigns dosage 0; 'error' aborts."
        )
    ),
    make_option(
        "--gvcf_mode",
        dest = "gvcf_mode",
        type = "character",
        default = "auto",
        help = paste0(
            "Asserts how gVCF status was decided upstream: 'auto', ",
            "'gvcf' or 'plain'. Cross-checked against --input_is_gvcf ",
            "and recorded in the summary; 'gvcf' or 'plain' disagreeing ",
            "with --input_is_gvcf is an error. [default: %default]"
        )
    ),
    make_option(
        "--default_ploidy",
        dest = "default_ploidy",
        type = "integer",
        default = 2L
    ),
    make_option(
        "--strict_alleles",
        dest = "strict_alleles",
        type = "character",
        default = "true"
    ),
    make_option(
        "--input_is_gvcf",
        dest = "input_is_gvcf",
        type = "character",
        default = "false"
    ),
    make_option(
        "--min_covered_fraction",
        dest = "min_covered_fraction",
        type = "double",
        default = 0,
        help = paste0(
            "Abort if the fraction of score variants with a usable ",
            "called SNP genotype or validated hom-ref reference block ",
            "falls below this value. No-calls, allele mismatches and ",
            "assumed missing genotypes are not covered. ",
            "0 disables this threshold; it does not turn unscored ",
            "variants into a complete PRS. [default: %default]"
        )
    ),
    make_option(
        "--output_summary",
        dest = "output_summary",
        type = "character"
    ),
    make_option(
        "--output_details",
        dest = "output_details",
        type = "character"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)

required_options <- c(
    "scorefile",
    "genotypes",
    "sample_id",
    "vcf_sample_file",
    "source_vcf",
    "target_build",
    "output_summary",
    "output_details"
)

missing_options <- required_options[
    vapply(
        required_options,
        function(option_name) {
            value <- opt[[option_name]]
            is.null(value) || !nzchar(value)
        },
        logical(1)
    )
]

if (length(missing_options) > 0L) {
    stopf(
        "Missing required option(s): %s",
        paste(missing_options, collapse = ", ")
    )
}

if (identical(opt$sample_id, "NA")) {
    stopf("Literal pipeline sample ID 'NA' is reserved")
}

parse_boolean <- function(value, option_name) {
    value <- tolower(trimws(as.character(value)))

    if (value == "true") {
        return(TRUE)
    }

    if (value == "false") {
        return(FALSE)
    }

    stopf("%s must be true or false", option_name)
}

strict_alleles <- parse_boolean(
    opt$strict_alleles,
    "--strict_alleles"
)

if (!strict_alleles) {
    stopf(
        "--strict_alleles false is unsupported: allele mismatches must abort, not yield NA dosages"
    )
}

input_is_gvcf <- parse_boolean(
    opt$input_is_gvcf,
    "--input_is_gvcf"
)


gvcf_mode_declared <- tolower(trimws(as.character(opt$gvcf_mode)))

if (!gvcf_mode_declared %in% c("auto", "gvcf", "plain")) {
    stopf(
        "--gvcf_mode must be auto, gvcf or plain, not '%s'",
        gvcf_mode_declared
    )
}

min_covered_fraction <- suppressWarnings(
    as.numeric(opt$min_covered_fraction)
)

if (
    length(min_covered_fraction) != 1L ||
    is.na(min_covered_fraction) ||
    !is.finite(min_covered_fraction) ||
    min_covered_fraction < 0 ||
    min_covered_fraction > 1
) {
    stopf("--min_covered_fraction must be a number between 0 and 1")
}

#
# --gvcf_mode records what the caller declared; --input_is_gvcf records
# what was resolved upstream. If they disagree, one of the two is wrong
# and the consequences are silent (a declared 'plain' input that is
# handled as a gVCF discards --missing_genotype without notice), so this
# is an error rather than a warning.
#
if (gvcf_mode_declared == "gvcf" && !input_is_gvcf) {
    stopf(
        paste0(
            "--gvcf_mode gvcf was declared but --input_is_gvcf is ",
            "false. Pass --input_is_gvcf true, or declare ",
            "--gvcf_mode auto and let the upstream detection decide."
        )
    )
}

if (gvcf_mode_declared == "plain" && input_is_gvcf) {
    stopf(
        paste0(
            "--gvcf_mode plain was declared but --input_is_gvcf is ",
            "true. Pass --input_is_gvcf false to force plain handling ",
            "(--missing_genotype then becomes mandatory), or declare ",
            "--gvcf_mode auto."
        )
    )
}

missing_supplied <- !is.null(opt$missing_genotype) &&
    !is.na(opt$missing_genotype) &&
    nzchar(trimws(as.character(opt$missing_genotype)))

requested_missing_mode <- if (missing_supplied) {
    tolower(trimws(as.character(opt$missing_genotype)))
} else {
    NA_character_
}

if (
    missing_supplied &&
    !requested_missing_mode %in% c("reference", "zero", "error")
) {
    stopf(
        "--missing_genotype must be reference, zero or error, not '%s'",
        requested_missing_mode
    )
}

if (!missing_supplied) {
    stopf(
        "--missing_genotype is required for both plain VCF and gVCF input; choose reference, zero or error"
    )
}

missing_mode <- requested_missing_mode

if (input_is_gvcf && missing_mode == "reference") {
    stopf(
        "--missing_genotype reference is allowed only for plain VCF input, not gVCF input"
    )
}

no_calls_supplied <- !is.null(opt$no_calls) &&
    length(opt$no_calls) == 1L &&
    !is.na(opt$no_calls) &&
    nzchar(trimws(opt$no_calls))

if (!no_calls_supplied) {
    stopf(
        "Missing required option: --no_calls. Choose zero or error."
    )
}

no_calls_mode <- tolower(trimws(opt$no_calls))

if (!(no_calls_mode %in% c("zero", "error"))) {
    stopf(
        "--no_calls must be zero or error, not '%s'",
        no_calls_mode
    )
}

default_ploidy <- as.integer(opt$default_ploidy)

if (
    is.na(default_ploidy) ||
    default_ploidy < 1L
) {
    stopf("--default_ploidy must be a positive integer")
}

target_build <- tolower(
    trimws(opt$target_build)
)

genome_package <- switch(
    target_build,
    hg38 = "BSgenome.Hsapiens.UCSC.hg38",
    hg37 = "BSgenome.Hsapiens.UCSC.hg19",
    hg19 = "BSgenome.Hsapiens.UCSC.hg19",
    stopf(
        "Unsupported target build '%s'",
        target_build
    )
)

suppressPackageStartupMessages(
    library(
        genome_package,
        character.only = TRUE
    )
)

genome_environment <- as.environment(
    paste0("package:", genome_package)
)

if (
    !exists(
        "Hsapiens",
        envir = genome_environment,
        inherits = FALSE
    )
) {
    stopf(
        "Could not load the Hsapiens object from %s",
        genome_package
    )
}

reference_genome <- get(
    "Hsapiens",
    envir = genome_environment,
    inherits = FALSE
)

available_seqnames <- names(reference_genome)

normalise_chromosome <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^chr", "", x, ignore.case = TRUE)
    x <- toupper(x)

    x[x %in% c("23")] <- "X"
    x[x %in% c("24")] <- "Y"
    x[x %in% c("M", "MTDNA", "25")] <- "MT"

    x
}

genome_seqname <- function(chromosome) {
    chromosome <- normalise_chromosome(chromosome)

    candidates <- if (chromosome == "MT") {
        c("chrM", "MT", "M")
    } else {
        c(
            paste0("chr", chromosome),
            chromosome
        )
    }

    matches <- candidates[
        candidates %in% available_seqnames
    ]

    if (length(matches) == 0L) {
        return(NA_character_)
    }

    matches[[1]]
}

reference_cache <- new.env(parent = emptyenv())

get_reference_sequence <- function(chromosome, position, width = 1L) {
    cache_key <- paste(chromosome, position, width, sep = ":")
    cached <- reference_cache[[cache_key]]

    if (!is.null(cached)) {
        return(cached)
    }

    seqname <- genome_seqname(chromosome)

    value <- if (is.na(seqname)) {
        NA_character_
    } else {
        tryCatch(
            {
                sequence <- Biostrings::getSeq(
                    reference_genome,
                    names = seqname,
                    start = as.integer(position),
                    width = as.integer(width)
                )

                toupper(as.character(sequence))
            },
            error = function(error) {
                NA_character_
            }
        )
    }

    reference_cache[[cache_key]] <- value
    value
}

normalise_allele <- function(x) {
    x <- toupper(trimws(as.character(x)))

    x[
        is.na(x) |
            x %in% c("", ".", "NA")
    ] <- NA_character_

    x
}

score_allele_set <- function(value) {
    value <- normalise_allele(value)

    if (length(value) != 1L || is.na(value)) {
        return(character(0))
    }

    if (!grepl("^[ACGT](/[ACGT])*$", value)) {
        stopf(
            "Invalid normalised score allele set '%s'",
            value
        )
    }

    unique(
        strsplit(
            value,
            "/",
            fixed = TRUE
        )[[1]]
    )
}

record_alleles <- function(ref, alt) {
    ref <- normalise_allele(ref)
    alt <- normalise_allele(alt)

    alternate_alleles <- character(0)

    if (!is.na(alt)) {
        alternate_alleles <- normalise_allele(
            strsplit(
                alt,
                ",",
                fixed = TRUE
            )[[1]]
        )
    }

    c(ref, alternate_alleles)
}

#
# A record carries no concrete alternate allele when ALT is missing or is
# a symbolic non-reference placeholder. Such records assert reference
# evidence over an interval and must never enter allele matching, because
# the effect allele cannot be present in their allele list.
#
reference_only_alt <- function(alt) {
    alt <- normalise_allele(alt)

    is.na(alt) | alt %in% c("<NON_REF>", "<*>")
}

symbolic_only_alt <- function(alt) {
    alt <- normalise_allele(alt)

    if (is.na(alt)) {
        return(TRUE)
    }

    parts <- normalise_allele(
        strsplit(alt, ",", fixed = TRUE)[[1]]
    )

    all(!grepl("^[ACGTN]+$", parts))
}

validate_snp_alleles <- function(
    ref,
    alt,
    effect_set,
    other_set,
    chromosome,
    position,
    variant_id
) {
    alleles <- record_alleles(ref, alt)

    # Symbolic gVCF placeholders do not assert a particular base.
    # A GT selecting one will be rejected separately below.
    concrete <- unique(
        alleles[
            !is.na(alleles) &
                grepl("^[ACGT]$", alleles)
        ]
    )

    if (length(other_set) > 0L) {
        unexpected <- setdiff(
            concrete,
            c(effect_set, other_set)
        )

        if (length(unexpected) > 0L) {
            stopf(
                paste0(
                    "Score/VCF allele mismatch for sample '%s' at ",
                    "%s:%d (variant '%s'): VCF allele(s) %s are ",
                    "absent from score allele set %s"
                ),
                opt$sample_id,
                chromosome,
                position,
                variant_id,
                paste(unexpected, collapse = "/"),
                paste(
                    c(effect_set, other_set),
                    collapse = "/"
                )
            )
        }
    } else if (
        length(intersect(effect_set, concrete)) == 0L &&
        !symbolic_only_alt(alt)
    ) {
        # Retain the existing effect-only model behaviour when the
        # score file provides no other_allele.
        stopf(
            paste0(
                "Score/VCF allele mismatch for sample '%s' at ",
                "%s:%d (variant '%s'): no effect allele %s occurs ",
                "in the VCF REF/ALT alleles"
            ),
            opt$sample_id,
            chromosome,
            position,
            variant_id,
            paste(effect_set, collapse = "/")
        )
    }

    alleles
}


# Only single-base A/C/G/T alleles and the two gVCF placeholders
# are SNP-compatible. Bare * represents a spanning deletion;
# other symbolic ALTs, such as <DEL>, are not SNP alleles.
is_snv_or_symbolic_record <- function(ref, alt, position, end) {
    ref <- normalise_allele(ref)
    alt <- normalise_allele(alt)

    if (is.na(ref) || !grepl("^[ACGT]$", ref)) {
        return(FALSE)
    }

    if (
        !is.na(end) &&
        end > position &&
        !reference_only_alt(alt)
    ) {
        return(FALSE)
    }

    if (is.na(alt)) {
        return(TRUE)
    }

    parts <- strsplit(alt, ",", fixed = TRUE)[[1]]

    all(
        grepl("^[ACGT]$", parts) |
            parts %in% c("<NON_REF>", "<*>")
    )
}


is_full_no_call <- function(gt) {
    if (is.na(gt) || !nzchar(gt)) {
        return(FALSE)
    }

    primary <- sub(":.*$", "", gt)

    # ., ./., .|., and wholly missing polyploid GTs.
    grepl("^\\.([/|]\\.)*$", primary)
}

check_nonoverlapping <- function(intervals, label) {
    if (nrow(intervals) < 2L) {
        return(invisible(NULL))
    }

    for (chromosome in unique(intervals$chrom)) {
        rows <- which(intervals$chrom == chromosome)

        if (length(rows) < 2L) {
            next
        }

        previous <- rows[-length(rows)]
        following <- rows[-1L]

        overlaps <- which(
            intervals$start[following] <=
                intervals$end[previous]
        )

        if (length(overlaps) > 0L) {
            first <- previous[[overlaps[[1]]]]
            second <- following[[overlaps[[1]]]]

            stopf(
                "Overlapping %s on %s: %d-%d and %d-%d",
                label,
                chromosome,
                intervals$start[[first]],
                intervals$end[[first]],
                intervals$start[[second]],
                intervals$end[[second]]
            )
        }
    }

    invisible(NULL)
}

fallback_genotype <- function(
    reason,
    variant_id,
    effect_is_reference,
    model_mismatch,
    reference_base,
    policy,
    ploidy = default_ploidy
) {
    if (
        length(policy) != 1L ||
        is.na(policy) ||
        !policy %in% c("reference", "zero", "error")
    ) {
        stopf(
            "Internal error: unsupported missing-genotype policy '%s'",
            as.character(policy)
        )
    }
    if (policy == "error") {
        stopf(
            paste0(
                "Variant '%s' has an unusable genotype/site (%s) and ",
                "--missing_genotype=error"
            ),
            variant_id,
            reason
        )
    }

    if (policy == "zero") {
        return(list(
            dosage = 0,
            status = paste0(reason, "_set_to_zero"),
            dosage_from_reference = FALSE
        ))
    }

    if (is.na(reference_base)) {
        stopf(
            "Cannot score variant '%s': its BSgenome reference base is unavailable",
            variant_id
        )
    }

    if (model_mismatch) {
        stopf(
            paste0(
                "Score/reference allele mismatch for variant '%s': ",
                "the BSgenome reference allele %s is in neither ",
                "score allele set"
            ),
            variant_id,
            reference_base
        )
    }

    dosage <- if (effect_is_reference) {
        ploidy
    } else {
        0
    }

    list(
        dosage = as.numeric(dosage),
        status = paste0(
            reason,
            "_dosage_from_reference"
        ),
        dosage_from_reference = TRUE
    )
}

genotype_header <- names(
    read_tsv(
        opt$genotypes,
        n_max = 0L,
        show_col_types = FALSE,
        progress = FALSE
    )
)

genotype_cols <- list(
    chrom = col_character(),
    position = col_integer(),
    ref = col_character(),
    alt = col_character(),
    genotype = col_character()
)

if ("end" %in% genotype_header) {
    genotype_cols$end <- col_character()
}

genotypes <- read_tsv(
    opt$genotypes,
    na = c("", "NA"),
    show_col_types = FALSE,
    progress = FALSE,
    col_types = do.call(cols, genotype_cols)
)

if ("end" %in% names(genotypes)) {
    raw_end <- trimws(genotypes$end)
    has_end <- !is.na(raw_end) & raw_end != "."

    bad_end <- has_end & !grepl("^[0-9]+$", raw_end)
    if (any(bad_end)) {
        stopf(
            "Invalid INFO/END at queried genotype row(s): %s",
            paste(head(which(bad_end), 20L), collapse = ", ")
        )
    }

    parsed_end <- suppressWarnings(as.numeric(raw_end[has_end]))
    if (any(
        !is.finite(parsed_end) |
            parsed_end < 1 |
            parsed_end > .Machine$integer.max
    )) {
        stopf("INFO/END is outside the supported integer range")
    }

    end_values <- rep(NA_integer_, nrow(genotypes))
    end_values[has_end] <- as.integer(parsed_end)
    genotypes$end <- end_values
} else {
    genotypes$end <- NA_integer_
}

scores <- read_tsv(
    opt$scorefile,
    na = c("", "NA"),
    show_col_types = FALSE,
    progress = FALSE,
    col_types = cols(
        score_row = col_integer(),
        variant_id = col_character(),
        chrom = col_character(),
        position = col_integer(),
        effect_allele = col_character(),
        other_allele = col_character(),
        effect_weight = col_double()
    )
)

if (nrow(scores) == 0L) {
    stopf("The normalised score file contains no variants")
}


vcf_sample_ids <- readLines(
    opt$vcf_sample_file,
    warn = FALSE
)

vcf_sample_ids <- vcf_sample_ids[
    nzchar(trimws(vcf_sample_ids))
]

if (length(vcf_sample_ids) != 1L) {
    stopf(
        "Expected exactly one VCF sample ID, found %d",
        length(vcf_sample_ids)
    )
}

vcf_sample_id <- vcf_sample_ids[[1]]

if (identical(vcf_sample_id, "NA")) {
    stopf("Literal VCF sample ID 'NA' is reserved")
}

scores <- scores %>%
    mutate(
        chrom = normalise_chromosome(chrom),
        effect_allele = normalise_allele(effect_allele),
        other_allele = normalise_allele(other_allele),
        key = paste(chrom, position, sep = ":")
    )

non_snv_at_score <- rep(NA_integer_, nrow(scores))

reference_blocks <- tibble(
    chrom = character(0),
    start = integer(0),
    end = integer(0),
    ref = character(0),
    alt = character(0),
    genotype = character(0),
    ploidy = integer(0)
)

no_call_intervals <- tibble(
    chrom = character(0),
    start = integer(0),
    end = integer(0),
    ref = character(0),
    genotype = character(0),
    source_row = integer(0)
)

block_index <- list()
no_call_index <- list()
genotype_index <- list()

if (nrow(genotypes) > 0L) {
    genotypes <- genotypes %>%
        mutate(
            chrom = normalise_chromosome(chrom),
            ref = normalise_allele(ref),
            alt = normalise_allele(alt),
            key = paste(chrom, position, sep = ":"),
            is_reference_only = reference_only_alt(alt),
            is_snv_compatible = mapply(
                is_snv_or_symbolic_record,
                ref,
                alt,
                position,
                end,
                USE.NAMES = FALSE
            ),
            full_no_call = vapply(
                genotype,
                is_full_no_call,
                logical(1)
            )
        )

    bad_end_rows <- which(
        !is.na(genotypes$end) &
            genotypes$end < genotypes$position
    )

    if (length(bad_end_rows) > 0L) {
        stopf(
            "INFO/END precedes POS at queried genotype row(s): %s",
            paste(head(bad_end_rows, 20L), collapse = ", ")
        )
    }

    n_ignored_non_snv <- sum(
        !genotypes$is_snv_compatible
    )

    if (n_ignored_non_snv > 0L) {

        non_snv_rows <- which(
            !genotypes$is_snv_compatible
        )

        ref_width <- nchar(genotypes$ref[non_snv_rows])
        ref_width[is.na(ref_width) | ref_width < 1L] <- 1L

        vcf_end <- genotypes$end[non_snv_rows]
        missing_end <- is.na(vcf_end)
        vcf_end[missing_end] <-
            genotypes$position[non_snv_rows][missing_end]

        span_end <- pmax(
            as.double(genotypes$position[non_snv_rows]) +
                ref_width - 1,
            as.double(vcf_end)
        )

        if (
            any(!is.finite(span_end)) ||
            any(span_end > .Machine$integer.max)
        ) {
            stopf("Invalid or out-of-range non-SNV VCF record span")
        }

        non_snv_intervals <- tibble(
            chrom = genotypes$chrom[non_snv_rows],
            start = genotypes$position[non_snv_rows],
            end = as.integer(span_end),
            source_row = as.integer(non_snv_rows)
        )

        for (
            chromosome in intersect(
                unique(scores$chrom),
                unique(non_snv_intervals$chrom)
            )
        ) {
            score_rows <- which(scores$chrom == chromosome)
            interval_rows <- which(
                non_snv_intervals$chrom == chromosome
            )

            hits <- IRanges::findOverlaps(
                IRanges::IRanges(
                    start = scores$position[score_rows],
                    width = 1L
                ),
                IRanges::IRanges(
                    start = non_snv_intervals$start[interval_rows],
                    end = non_snv_intervals$end[interval_rows]
                ),
                select = "first"
            )

            matched <- which(!is.na(hits))

            non_snv_at_score[score_rows[matched]] <-
                non_snv_intervals$source_row[
                    interval_rows[hits[matched]]
                ]
        }

        warning(
            sprintf(
                paste0(
                    "Sample '%s': ignored %d non-SNV VCF record(s) ",
                    "when selecting SNP genotypes."
                ),
                opt$sample_id,
                n_ignored_non_snv
            ),
            call. = FALSE
        )
    }

    genotype_primary <- sub(
        ":.*$",
        "",
        ifelse(
            is.na(genotypes$genotype),
            "",
            genotypes$genotype
        )
    )

    is_hom_ref_call <- grepl(
        "^0([/|]0)*$",
        genotype_primary
    )

    bad_interval_rows <- which(
        genotypes$is_snv_compatible &
            genotypes$is_reference_only &
            !is.na(genotypes$end) &
            genotypes$end > genotypes$position &
            !(is_hom_ref_call | genotypes$full_no_call)
    )

    if (length(bad_interval_rows) > 0L) {
        row <- bad_interval_rows[[1]]

        stopf(
            paste0(
                "Reference-only interval at %s:%d has unsupported ",
                "GT '%s'; an interval must be hom-ref or a full no-call"
            ),
            genotypes$chrom[[row]],
            genotypes$position[[row]],
            ifelse(
                is.na(genotypes$genotype[[row]]),
                "NA",
                genotypes$genotype[[row]]
            )
        )
    }

    block_rows <- which(
        genotypes$is_snv_compatible &
            genotypes$is_reference_only &
            is_hom_ref_call
    )

    if (length(block_rows) > 0L) {
        block_end <- ifelse(
            is.na(genotypes$end[block_rows]),
            genotypes$position[block_rows],
            genotypes$end[block_rows]
        )

        if (
            any(
                block_end <
                    genotypes$position[block_rows]
            )
        ) {
            stopf(
                "A hom-ref reference block has END before POS"
            )
        }

        block_ploidy <- vapply(
            genotype_primary[block_rows],
            function(gt) {
                length(
                    strsplit(
                        gt,
                        "[/|]",
                        perl = TRUE
                    )[[1]]
                )
            },
            integer(1)
        )

        reference_blocks <- tibble(
            chrom = genotypes$chrom[block_rows],
            start = genotypes$position[block_rows],
            end = as.integer(block_end),
            ref = genotypes$ref[block_rows],
            alt = genotypes$alt[block_rows],
            genotype = genotypes$genotype[block_rows],
            ploidy = as.integer(block_ploidy)
        ) %>%
            distinct(
                chrom, start, end, ref, alt, genotype, ploidy,
                .keep_all = TRUE
            ) %>%
            arrange(chrom, start, end)

        check_nonoverlapping(
            reference_blocks,
            "hom-ref reference blocks"
        )

        # Validate REF at the ACTUAL BLOCK START, even if no score
        # variant is located at that start.
        for (block in seq_len(nrow(reference_blocks))) {
            expected <- get_reference_sequence(
                reference_blocks$chrom[[block]],
                reference_blocks$start[[block]],
                nchar(reference_blocks$ref[[block]])
            )

            if (is.na(expected)) {
                stopf(
                    "Cannot verify reference block REF at %s:%d",
                    reference_blocks$chrom[[block]],
                    reference_blocks$start[[block]]
                )
            }

            if (
                !identical(
                    reference_blocks$ref[[block]],
                    expected
                )
            ) {
                stopf(
                    paste0(
                        "Reference block REF mismatch at %s:%d: ",
                        "VCF REF=%s; %s REF=%s"
                    ),
                    reference_blocks$chrom[[block]],
                    reference_blocks$start[[block]],
                    reference_blocks$ref[[block]],
                    target_build,
                    expected
                )
            }
        }

        block_index <- split(
            seq_len(nrow(reference_blocks)),
            reference_blocks$chrom
        )
    }

    # Reference-only ./., including ALT=., is NOT a hom-ref block.
    # Index it as a no-call at its start. When END spans an interval,
    # also recognise its internal positions as no-calls.
    no_call_rows <- which(
        genotypes$is_snv_compatible &
            genotypes$is_reference_only &
            genotypes$full_no_call &
            !is.na(genotypes$end)
    )

    if (length(no_call_rows) > 0L) {
        if (
            any(
                genotypes$end[no_call_rows] <
                    genotypes$position[no_call_rows]
            )
        ) {
            stopf(
                "A no-call interval has END before POS"
            )
        }

        interval_rows <- no_call_rows[
            genotypes$end[no_call_rows] >
                genotypes$position[no_call_rows]
        ]

        if (length(interval_rows) > 0L) {
            no_call_intervals <- tibble(
                chrom = genotypes$chrom[interval_rows],
                start = genotypes$position[interval_rows],
                end = genotypes$end[interval_rows],
                ref = genotypes$ref[interval_rows],
                genotype = genotypes$genotype[interval_rows],
                source_row = as.integer(interval_rows)
            ) %>%
                distinct(
                    chrom,
                    start,
                    end,
                    ref,
                    genotype,
                    .keep_all = TRUE
                ) %>%
                arrange(chrom, start, end)

            check_nonoverlapping(
                no_call_intervals,
                "no-call intervals"
            )

            for (interval in seq_len(nrow(no_call_intervals))) {
                expected <- get_reference_sequence(
                    no_call_intervals$chrom[[interval]],
                    no_call_intervals$start[[interval]]
                )

                if (
                    is.na(expected) ||
                    !identical(
                        expected,
                        no_call_intervals$ref[[interval]]
                    )
                ) {
                    stopf(
                        "No-call interval REF mismatch at %s:%d",
                        no_call_intervals$chrom[[interval]],
                        no_call_intervals$start[[interval]]
                    )
                }
            }

            no_call_index <- split(
                seq_len(nrow(no_call_intervals)),
                no_call_intervals$chrom
            )
        }
    }

    # Includes explicit no-calls and called symbolic ALT records;
    # excludes validated hom-ref blocks and concrete indels.
    variant_rows <- which(
        genotypes$is_snv_compatible &
            !(
                genotypes$is_reference_only &
                    is_hom_ref_call
            )
    )

    if (length(variant_rows) > 0L) {
        genotype_index <- split(
            variant_rows,
            genotypes$key[variant_rows]
        )
    }

    # A called variant or explicit no-call contradicts an overlapping
    # interval that asserts homozygous reference.
    if (nrow(reference_blocks) > 0L) {
        for (row in variant_rows) {
            overlaps <- (
                reference_blocks$chrom ==
                    genotypes$chrom[[row]] &
                    reference_blocks$start <=
                    genotypes$position[[row]] &
                    reference_blocks$end >=
                    genotypes$position[[row]]
            )

            if (any(overlaps)) {
                stopf(
                    paste0(
                        "VCF record at %s:%d overlaps a hom-ref ",
                        "reference block"
                    ),
                    genotypes$chrom[[row]],
                    genotypes$position[[row]]
                )
            }
        }
    }

    for (interval in seq_len(nrow(no_call_intervals))) {
        chromosome <- no_call_intervals$chrom[[interval]]
        first <- no_call_intervals$start[[interval]]
        last <- no_call_intervals$end[[interval]]

        overlapping_blocks <- (
            reference_blocks$chrom == chromosome &
                reference_blocks$start <= last &
                reference_blocks$end >= first
        )

        if (any(overlapping_blocks)) {
            stopf(
                "No-call interval %s:%d-%d overlaps a hom-ref block",
                chromosome,
                first,
                last
            )
        }

        other_rows <- variant_rows[
            variant_rows !=
                no_call_intervals$source_row[[interval]]
        ]

        contradictory_records <- (
            genotypes$chrom[other_rows] == chromosome &
                genotypes$position[other_rows] >= first &
                genotypes$position[other_rows] <= last
        )

        if (any(contradictory_records)) {
            stopf(
                paste0(
                    "No-call interval %s:%d-%d overlaps ",
                    "another VCF record"
                ),
                chromosome,
                first,
                last
            )
        }
    }
}

n_homozygous_reference_blocks <- sum(
    reference_blocks$end > reference_blocks$start
)
n_hom_ref_site_records <- sum(
    reference_blocks$end == reference_blocks$start
)

#
# Blocks within a chromosome are sorted by start and, in a well-formed
# gVCF, non-overlapping, so a binary search on the start coordinates
# locates the only candidate block.
#
find_reference_block <- function(chromosome, position) {
    rows <- block_index[[chromosome]]

    if (is.null(rows) || length(rows) == 0L) {
        return(NA_integer_)
    }

    candidate <- findInterval(
        position,
        reference_blocks$start[rows]
    )

    if (candidate < 1L) {
        return(NA_integer_)
    }

    row <- rows[[candidate]]

    if (reference_blocks$end[[row]] >= position) {
        return(as.integer(row))
    }

    NA_integer_
}

find_no_call_interval <- function(chromosome, position) {
    rows <- no_call_index[[chromosome]]

    if (is.null(rows) || length(rows) == 0L) {
        return(NA_integer_)
    }

    candidate <- findInterval(
        position,
        no_call_intervals$start[rows]
    )

    if (candidate < 1L) {
        return(NA_integer_)
    }

    row <- rows[[candidate]]

    if (no_call_intervals$end[[row]] >= position) {
        return(as.integer(row))
    }

    NA_integer_
}

describe_reference_block <- function(block_row, position) {
    if (is.na(block_row)) {
        return(list(
            found = FALSE,
            text = NA_character_,
            ref = NA_character_
        ))
    }

    list(
        found = TRUE,
        text = sprintf(
            "%s:%d-%d",
            reference_blocks$chrom[[block_row]],
            reference_blocks$start[[block_row]],
            reference_blocks$end[[block_row]]
        ),
        # REF belongs to the block START, not necessarily to the
        # score position reported on this detail row.
        ref = reference_blocks$ref[[block_row]]
    )
}

resolve_missing <- function(
    reason,
    variant_id,
    effect_is_reference,
    model_mismatch,
    reference_base,
    block_row
) {
    # A validated block is observed genotype evidence, not missingness.
    if (!is.na(block_row)) {
        is_spanning_block <- (
            reference_blocks$end[[block_row]] >
                reference_blocks$start[[block_row]]
        )

        result <- fallback_genotype(
            reason = "observed_hom_ref",
            variant_id = variant_id,
            effect_is_reference = effect_is_reference,
            model_mismatch = model_mismatch,
            reference_base = reference_base,
            policy = "reference",
            ploidy = reference_blocks$ploidy[[block_row]]
        )

        result$status <- if (is_spanning_block) {
            "observed_hom_ref_block"
        } else {
            "observed_hom_ref_site"
        }

        # A VCF GT supplied this dosage. It was not an assumed
        # reference genotype at an absent position.
        result$dosage_from_reference <- FALSE
        return(result)
    }

    if (input_is_gvcf && missing_mode == "error") {
        stopf(
            "Variant '%s' is absent/uncovered in this gVCF (%s) and --missing_genotype=error",
            variant_id,
            reason
        )
    }

    fallback_genotype(
        reason = if (input_is_gvcf) {
            paste0(reason, "_uncovered")
        } else {
            reason
        },
        variant_id = variant_id,
        effect_is_reference = effect_is_reference,
        model_mismatch = model_mismatch,
        reference_base = reference_base,
        policy = missing_mode
    )
}

detail_rows <- vector(
    "list",
    nrow(scores)
)

for (i in seq_len(nrow(scores))) {
    score_row <- scores[i, ]

    chromosome <- score_row$chrom[[1]]
    position <- score_row$position[[1]]
    variant_id <- score_row$variant_id[[1]]
    effect_allele <- score_row$effect_allele[[1]]
    other_allele <- score_row$other_allele[[1]]
    effect_weight <- score_row$effect_weight[[1]]
    key <- score_row$key[[1]]
    partial_no_call <- FALSE
    gt_known_alleles <- NA_integer_
    gt_missing_alleles <- NA_integer_
    vcf_ref_relation <- NA_character_

    if (
        is.na(other_allele) ||
        !nzchar(other_allele)
    ) {
        other_allele <- NA_character_
    }

    effect_set <- score_allele_set(effect_allele)
    other_set <- score_allele_set(other_allele)

    if (
        length(effect_set) == 0L ||
        length(intersect(effect_set, other_set)) > 0L
    ) {
        stopf(
            "Invalid or overlapping score allele sets for variant '%s'",
            variant_id
        )
    }

    if (
        is.na(effect_weight) ||
        !is.finite(effect_weight)
    ) {
        stopf(
            "Invalid effect weight for variant '%s'",
            variant_id
        )
    }

    reference_base <- get_reference_sequence(
        chromosome,
        position,
        1L
    )

    if (
        is.na(reference_base) ||
        !(reference_base %in% c("A", "C", "G", "T"))
    ) {
        stopf(
            "Cannot determine a concrete BSgenome reference base for score variant '%s' at %s:%d",
            variant_id,
            chromosome,
            position
        )
    }

    effect_is_reference <- reference_base %in% effect_set
    other_is_reference <- reference_base %in% other_set

    model_mismatch <- (
        length(other_set) > 0L &&
            !effect_is_reference &&
            !other_is_reference
    )

    non_snv_row <- non_snv_at_score[[i]]

    if (!is.na(non_snv_row)) {
        offending <- genotypes[non_snv_row, , drop = FALSE]

        stopf(
            paste0(
                "Score SNP '%s' at %s:%d overlaps a non-SNV VCF ",
                "record starting at %s:%d (REF=%s, ALT=%s, ",
                "INFO/END=%s); it cannot be scored as a SNP or ",
                "treated as absent"
            ),
            variant_id,
            chromosome,
            position,
            offending$chrom[[1]],
            offending$position[[1]],
            offending$ref[[1]],
            ifelse(is.na(offending$alt[[1]]),
                   ".", offending$alt[[1]]),
            ifelse(is.na(offending$end[[1]]),
                   ".", as.character(offending$end[[1]]))
        )
    }

    record_indices <- genotype_index[[key]]

    interval_row <- find_no_call_interval(
        chromosome, position
    )
    block_row <- find_reference_block(
        chromosome, position
    )

    exact_hom_ref_record <- (
        !is.na(block_row) &&
        reference_blocks$start[[block_row]] == position
    )

    # An interval record synthesized for an interior position is
    # evidence, but it is not a record beginning at that position.
    vcf_record_at_position <- (
        length(record_indices) > 0L ||
        exact_hom_ref_record
    )

    no_call_interval_found <- !is.na(interval_row)

    if (
        length(record_indices) == 0L &&
        no_call_interval_found
    ) {
        record_indices <-
            no_call_intervals$source_row[[interval_row]]
    }

    # Retain this internal variable for choosing the scoring branch.
    record_found <- length(record_indices) > 0L

    record_ambiguous <- FALSE
    vcf_ref <- NA_character_
    vcf_alt <- NA_character_
    genotype <- NA_character_
    called_alleles_text <- NA_character_
    dosage <- NA_real_
    contribution <- NA_real_
    status <- NA_character_
    no_call <- FALSE
    non_model_call <- FALSE
    dosage_from_reference <- FALSE
    allele_mismatch <- FALSE

    reference_block_found <- (
        !is.na(block_row) &&
        reference_blocks$end[[block_row]] >
            reference_blocks$start[[block_row]]
    )
    hom_ref_site_found <- (
        !is.na(block_row) &&
        reference_blocks$end[[block_row]] ==
            reference_blocks$start[[block_row]]
    )
    reference_block_text <- NA_character_
    reference_block_ref <- NA_character_

    if (record_found && !is.na(block_row)) {
        stopf(
            paste0(
                "VCF record at %s:%d contradicts an ",
                "overlapping hom-ref block"
            ),
            chromosome,
            position
        )
    }

    if (record_found) {
        records <- genotypes[
            record_indices,
            ,
            drop = FALSE
        ]

        # In particular, two different SNP GTs at the same canonical
        # chr22/22 position must NOT be decided by input order or by
        # the previous compatibility-score tie breaker.
        signatures <- paste(
            records$ref,
            ifelse(
                is.na(records$alt),
                ".",
                records$alt
            ),
            ifelse(
                is.na(records$genotype),
                ".",
                records$genotype
            ),
            sep = ":"
        )

        if (length(unique(signatures)) > 1L) {
            stopf(
                paste0(
                    "Conflicting VCF records at score position ",
                    "%s:%d (variant '%s')"
                ),
                chromosome,
                position,
                variant_id
            )
        }

        # Identical repeated records are unambiguous. Conflicting
        # ones have already failed, regardless of their order.
        selected_record <- records[1L, , drop = FALSE]
        genotype <- selected_record$genotype[[1]]

        if (
            selected_record$position[[1]] ==
                position
        ) {
            vcf_ref <- selected_record$ref[[1]]
            vcf_alt <- selected_record$alt[[1]]
        }

        # An interior score position in a no-call interval has no VCF REF
        # of its own; the interval's start REF was validated separately.
        if (!is.na(vcf_ref)) {
            if (vcf_ref == reference_base) {
                vcf_ref_relation <- "reference_match"
            } else {
                vcf_alts <- record_alleles(
                    vcf_ref,
                    vcf_alt
                )[-1L]

                concrete_alts <- vcf_alts[
                    !is.na(vcf_alts) &
                        grepl("^[ACGT]$", vcf_alts)
                ]

                if (reference_base %in% concrete_alts) {
                    vcf_ref_relation <- "swap_compatible"
                } else {
                    stopf(
                        paste0(
                            "Unresolved VCF/reference mismatch for sample ",
                            "'%s', score row %d at %s:%d: BSgenome REF=%s, ",
                            "VCF REF=%s, VCF ALT=%s. The BSgenome base is ",
                            "absent from the concrete VCF alleles; this ",
                            "cannot be explained by a REF/ALT swap."
                        ),
                        opt$sample_id,
                        score_row$score_row[[1]],
                        chromosome,
                        position,
                        reference_base,
                        vcf_ref,
                        ifelse(is.na(vcf_alt), ".", vcf_alt)
                    )
                }
            }
        }

        # For the INSIDE of a no-call interval, its start REF is not
        # the REF base at this scored position. Leave vcf_ref as NA.

        if (is_full_no_call(genotype)) {
            if (no_calls_mode == "error") {
                stopf(
                    "Score SNP '%s' at %s:%d has GT=%s and --no_calls=error",
                    variant_id,
                    chromosome,
                    position,
                    genotype
                )
            }

            if (!is.na(vcf_ref)) {
                validate_snp_alleles(
                    vcf_ref,
                    vcf_alt,
                    effect_set,
                    other_set,
                    chromosome,
                    position,
                    variant_id
                )
            }

            no_call <- TRUE
            dosage <- 0
            status <- "no_call_set_to_zero"
        } else {
            if (is.na(genotype) || !nzchar(genotype)) {
                stopf(
                    "Score SNP '%s' at %s:%d has no usable GT field",
                    variant_id,
                    chromosome,
                    position
                )
            }

            if (is.na(vcf_ref)) {
                stopf(
                    "Score SNP '%s' at %s:%d has an unsupported interval genotype '%s'",
                    variant_id,
                    chromosome,
                    position,
                    genotype
                )
            }

            alleles <- validate_snp_alleles(
                vcf_ref,
                vcf_alt,
                effect_set,
                other_set,
                chromosome,
                position,
                variant_id
            )

            primary_gt <- sub(":.*$", "", genotype)

            # Accept numeric allele indexes and '.', in any ploidy.
            # A wholly missing GT was handled above.
            valid_gt <- grepl(
                "^(\\.|[0-9]+)([/|](\\.|[0-9]+))*$",
                primary_gt
            )

            if (!valid_gt) {
                stopf(
                    "Invalid SNP genotype '%s' for variant '%s' at %s:%d",
                    genotype,
                    variant_id,
                    chromosome,
                    position
                )
            }

            gt_tokens <- strsplit(
                primary_gt,
                "[/|]",
                perl = TRUE
            )[[1]]

            missing_token <- gt_tokens == "."

            if (all(missing_token)) {
                stopf(
                    "Internal error: full no-call GT '%s' was not recognised",
                    genotype
                )
            }

            partial_no_call <- any(missing_token)

            if (partial_no_call && no_calls_mode == "error") {
                stopf(
                    paste0(
                        "Score SNP '%s' at %s:%d has partially missing ",
                        "GT=%s and --no_calls=error"
                    ),
                    variant_id,
                    chromosome,
                    position,
                    genotype
                )
            }

            # Only the called indexes are converted to alleles. Missing entries
            # have no effect-allele dosage under --no_calls=zero.
            known_indices <- suppressWarnings(
                as.integer(gt_tokens[!missing_token])
            )

            if (
                anyNA(known_indices) ||
                any(known_indices < 0L) ||
                any(known_indices >= length(alleles))
            ) {
                stopf(
                    "Invalid SNP genotype '%s' for variant '%s' at %s:%d",
                    genotype,
                    variant_id,
                    chromosome,
                    position
                )
            }

            called_alleles <- alleles[known_indices + 1L]

            if (
                any(
                    is.na(called_alleles) |
                        !grepl("^[ACGT]$", called_alleles)
                )
            ) {
                stopf(
                    paste0(
                        "Cannot score symbolic/unknown ALT selected ",
                        "by GT '%s' for variant '%s' at %s:%d"
                    ),
                    genotype,
                    variant_id,
                    chromosome,
                    position
                )
            }

            displayed_alleles <- rep(".", length(gt_tokens))
            displayed_alleles[!missing_token] <- called_alleles
            called_alleles_text <- paste(
                displayed_alleles,
                collapse = "/"
            )

            gt_known_alleles <- length(known_indices)
            gt_missing_alleles <- sum(missing_token)

            non_model_call <- (
                length(other_set) == 0L &&
                    any(!(called_alleles %in% effect_set))
            )

            dosage <- as.numeric(
                sum(called_alleles %in% effect_set)
            )

            status <- if (partial_no_call) {
                "partial_no_call_counted_known_alleles"
            } else if (non_model_call) {
                "observed_nonmodel_allele"
            } else {
                "observed"
            }
        }
    } else {
        fallback <- resolve_missing(
            reason = "absent",
            variant_id = variant_id,
            effect_is_reference = effect_is_reference,
            model_mismatch = model_mismatch,
            reference_base = reference_base,
            block_row = block_row
        )

        dosage <- fallback$dosage
        status <- fallback$status
        dosage_from_reference <-
            fallback$dosage_from_reference

        block_details <- describe_reference_block(
            block_row, position
        )

        if (reference_block_found) {
            reference_block_text <- block_details$text
            reference_block_ref <- block_details$ref
        }

        if (!is.na(block_row)) {
            genotype <- reference_blocks$genotype[[block_row]]
            gt_known_alleles <-
                reference_blocks$ploidy[[block_row]]
            gt_missing_alleles <- 0L
            called_alleles_text <- paste(
                rep(reference_base, gt_known_alleles),
                collapse = "/"
            )

            if (exact_hom_ref_record) {
                vcf_ref <- reference_blocks$ref[[block_row]]
                vcf_alt <- reference_blocks$alt[[block_row]]
                vcf_ref_relation <- "reference_match"
            }
        }
    }
    
    if (!is.na(dosage)) {
        contribution <- effect_weight * dosage
    }

    detail_rows[[i]] <- tibble(
        sample_id = opt$sample_id,
        vcf_sample_id = vcf_sample_id,
        source_vcf = opt$source_vcf,
        target_build = target_build,
        score_row = score_row$score_row[[1]],
        variant_id = variant_id,
        chrom = chromosome,
        position = position,
        effect_allele = effect_allele,
        other_allele = other_allele,
        effect_weight = effect_weight,
        reference_allele = reference_base,
        effect_allele_is_reference =
            effect_is_reference,
        other_allele_is_reference =
            other_is_reference,
        vcf_record_found = vcf_record_at_position,
        reference_block_found = reference_block_found,
        reference_block = reference_block_text,
        reference_block_ref = reference_block_ref,
        hom_ref_site_found = hom_ref_site_found,
        no_call_interval_found = no_call_interval_found,
        vcf_ref = vcf_ref,
        vcf_alt = vcf_alt,
        genotype = genotype,
        called_alleles = called_alleles_text,
        vcf_no_call = no_call,
        vcf_partial_no_call = partial_no_call,
        gt_known_alleles = gt_known_alleles,
        gt_missing_alleles = gt_missing_alleles,
        vcf_ref_relation = vcf_ref_relation,
        non_model_allele_call = non_model_call,
        dosage_from_reference =
            dosage_from_reference,
        effect_allele_dosage = dosage,
        contribution = contribution,
        status = status
    )
}

details <- bind_rows(detail_rows) %>%
    arrange(score_row)

if (
    anyNA(details$effect_allele_dosage) ||
    anyNA(details$contribution) ||
    any(!is.finite(details$effect_allele_dosage)) ||
    any(!is.finite(details$contribution))
) {
    stopf(
        "Internal error: an SNP has a missing or non-finite dosage/contribution"
    )
}

write_tsv(
    details,
    opt$output_details,
    na = "NA"
)

#
# Shared counts, computed once so the warnings and the summary cannot
# drift apart.
#

n_reference_lookup_failed <- sum(is.na(details$reference_allele))

# A bad block REF now fails before scoring, rather than being counted
# and merely warned about after publishing a PRS.
n_block_ref_not_genome_ref <- 0L

if (n_reference_lookup_failed != 0L) {
    stopf(
        "Internal error: a score SNP passed scoring without a BSgenome reference base"
    )
}

n_vcf_records_found <- sum(details$vcf_record_found)

n_positions_from_reference_block <- sum(
    details$reference_block_found
)

n_positions_from_hom_ref_site <- sum(
    details$hom_ref_site_found
)

n_uncovered_positions <- sum(
    !details$vcf_record_found &
        !details$reference_block_found &
        !details$no_call_interval_found
)

n_assumed_reference <- sum(details$dosage_from_reference)


usable_hom_ref_evidence <- details$status %in% c(
    "observed_hom_ref_block",
    "observed_hom_ref_site"
)

# An exact hom-ref record is also an observed GT at that position;
# a block covering only an interior position is not an exact record.
usable_observed_call <- (
    details$status %in% c(
        "observed",
        "observed_nonmodel_allele"
    ) |
    (
        usable_hom_ref_evidence &
            details$vcf_record_found
    )
)

usable_reference_block <- (
    details$reference_block_found &
        details$status ==
        "absent_reference_block_dosage_from_reference"
)

usable_coverage <- (
    usable_observed_call |
        usable_hom_ref_evidence
)

n_unusable_positions <- sum(!usable_coverage)

covered_fraction <- mean(
    usable_coverage
)

partial_rows <- which(details$vcf_partial_no_call)
n_vcf_partial_no_calls <- length(partial_rows)

n_known_alleles_in_partial_calls <- sum(
    details$gt_known_alleles[partial_rows]
)

n_missing_alleles_in_partial_calls <- sum(
    details$gt_missing_alleles[partial_rows]
)

partial_no_call_warning <- if (n_vcf_partial_no_calls > 0L) {
    sprintf(
        paste0(
            "Sample '%s': --no_calls=zero partially scored %d score ",
            "row(s) (first score_row IDs: %s). Only known genotype ",
            "alleles were used to count effect copies (%d known ",
            "copy/copies); %d missing copy/copies contributed zero. ",
            "These SNP contributions are incomplete and are excluded ",
            "from covered_fraction."
        ),
        opt$sample_id,
        n_vcf_partial_no_calls,
        paste(
            head(details$score_row[partial_rows], 20L),
            collapse = ","
        ),
        n_known_alleles_in_partial_calls,
        n_missing_alleles_in_partial_calls
    )
} else {
    NA_character_
}

if (!is.na(partial_no_call_warning)) {
    warning(partial_no_call_warning, call. = FALSE)
}


#
# gVCF declaration versus what the records actually look like. Advisory
# only: a gVCF whose blocks all happen to miss the score positions is
# legal; 'plain' mode is rejected if hom-ref blocks are detected.
#
if (
    input_is_gvcf &&
    n_positions_from_reference_block == 0L &&
    n_uncovered_positions > 0L
) {
    warning(
        sprintf(
            paste0(
                "Sample '%s' is treated as a gVCF, but no hom-ref ",
                "reference block overlapped any score position, and ",
                "%d position(s) are reported as uncovered. If this is ",
                "really a sites-only VCF, use --input_is_gvcf false ",
                "with an explicit --missing_genotype instead."
            ),
            opt$sample_id,
            n_uncovered_positions
        ),
        call. = FALSE
    )
}

if (!input_is_gvcf && n_positions_from_reference_block > 0L) {
    warning(
        sprintf(
            paste0(
                "Sample '%s' is treated as a plain VCF, but %d score ",
                "position(s) were resolved from a hom-ref reference ",
                "block. Those blocks are honoured as observed hom-ref ",
                "evidence; positions in no block are resolved with ",
                "--missing_genotype=%s. If this is an all-sites or ",
                "gVCF input, use --input_is_gvcf true."
            ),
            opt$sample_id,
            n_positions_from_reference_block,
            missing_mode
        ),
        call. = FALSE
    )
}

# Successful rows with a VCF REF/BSgenome disagreement have already
# been checked: the BSgenome base must occur among concrete VCF ALTs.
n_vcf_ref_swap_compatible <- sum(
    details$vcf_ref_relation == "swap_compatible",
    na.rm = TRUE
)

# Retain the existing QC column for compatibility. On successful
# samples it equals the swap-compatible count: unresolved cases abort.
n_vcf_ref_not_genome_ref <- n_vcf_ref_swap_compatible

# Only call a score row biallelic-palindromic when each score allele
# is one base. Slash-separated sets require separate consideration.
palindromic_biallelic <- (
    (details$effect_allele == "A" &
         details$other_allele == "T") |
    (details$effect_allele == "T" &
         details$other_allele == "A") |
    (details$effect_allele == "C" &
         details$other_allele == "G") |
    (details$effect_allele == "G" &
         details$other_allele == "C")
)

n_palindromic_biallelic <- sum(
    palindromic_biallelic,
    na.rm = TRUE
)

n_palindromic_swap_compatible <- sum(
    details$vcf_ref_relation == "swap_compatible" &
        palindromic_biallelic,
    na.rm = TRUE
)

if (n_vcf_ref_swap_compatible > 0L) {
    orientation_message <- sprintf(
        paste0(
            "Sample '%s': %d score row(s) have VCF REF differing ",
            "from the %s BSgenome base, but that base occurs among ",
            "their concrete VCF ALTs (swap-compatible, not proof ",
            "of correct build or strand)."
        ),
        opt$sample_id,
        n_vcf_ref_swap_compatible,
        target_build
    )

    if (input_is_gvcf) {
        orientation_message <- paste(
            orientation_message,
            "A gVCF REF disagreement is unexpected: investigate",
            "the calling FASTA, target build and coordinates."
        )
    }

    if (
        missing_mode == "reference" &&
        n_assumed_reference > 0L
    ) {
        orientation_message <- paste(
            orientation_message,
            sprintf(
                paste0(
                    "%d other score row(s) used assumed reference ",
                    "genotypes; those assumptions may be unreliable ",
                    "until reference orientation is resolved."
                ),
                n_assumed_reference
            )
        )
    }

    if (n_palindromic_swap_compatible > 0L) {
        orientation_message <- paste(
            orientation_message,
            sprintf(
                paste0(
                    "%d swap-compatible row(s) are biallelic A/T ",
                    "or C/G; their strand cannot be established ",
                    "from these allele labels alone."
                ),
                n_palindromic_swap_compatible
            )
        )
    }

    warning(orientation_message, call. = FALSE)
}

# This identifies a model/reference disagreement. It does NOT
# establish that the score and VCF have been jointly strand-flipped.
n_model_ref_allele_discrepancies <- sum(
    !is.na(details$other_allele) &
        !details$effect_allele_is_reference &
        !details$other_allele_is_reference,
    na.rm = TRUE
)

if (n_model_ref_allele_discrepancies > 0L) {
    warning(
        sprintf(
            paste0(
                "Sample '%s': %d score row(s) specify neither ",
                "score allele set as the %s BSgenome reference ",
                "base. This is an unresolved score/reference ",
                "discrepancy, not proof of a strand flip or of ",
                "correct called dosages. Verify build, position ",
                "and allele harmonization."
            ),
            opt$sample_id,
            n_model_ref_allele_discrepancies,
            target_build
        ),
        call. = FALSE
    )
}

n_scored <- nrow(details)
prs <- sum(details$contribution)
total_effect_alleles <- sum(details$effect_allele_dosage)

if (
    !is.finite(prs) ||
    !is.finite(total_effect_alleles)
) {
    stopf(
        "The PRS or total effect-allele dosage is non-finite"
    )
}

summary <- tibble(
    sample_id = opt$sample_id,
    vcf_sample_id = vcf_sample_id,
    source_vcf = opt$source_vcf,
    target_build = target_build,
    n_score_variants = nrow(details),
    n_scored_variants = n_scored,
    n_vcf_records_found = n_vcf_records_found,
    n_absent_from_vcf = sum(!details$vcf_record_found),
    n_vcf_no_calls = sum(details$vcf_no_call),
    n_vcf_partial_no_calls = n_vcf_partial_no_calls,
    n_known_alleles_in_partial_calls = n_known_alleles_in_partial_calls,
    n_missing_alleles_in_partial_calls = n_missing_alleles_in_partial_calls,
    partial_no_call_warning = partial_no_call_warning,
    n_dosage_from_reference = n_assumed_reference,
    n_non_model_allele_calls = sum(details$non_model_allele_call),
    n_vcf_ref_not_genome_ref = n_vcf_ref_not_genome_ref,
    n_vcf_ref_swap_compatible = n_vcf_ref_swap_compatible,
    n_model_ref_allele_discrepancies = n_model_ref_allele_discrepancies,
    n_palindromic_biallelic = n_palindromic_biallelic,
    n_palindromic_swap_compatible = n_palindromic_swap_compatible,
    n_block_ref_not_genome_ref = n_block_ref_not_genome_ref,
    input_is_gvcf = input_is_gvcf,
    gvcf_mode_declared = gvcf_mode_declared,
    missing_genotype_requested = requested_missing_mode,
    missing_genotype_mode = missing_mode,
    no_calls_mode = no_calls_mode,
    min_covered_fraction = min_covered_fraction,
    n_homozygous_reference_blocks = n_homozygous_reference_blocks,
    n_positions_from_reference_block = n_positions_from_reference_block,
    n_hom_ref_site_records = n_hom_ref_site_records,
    n_positions_from_hom_ref_site = n_positions_from_hom_ref_site,
    n_reference_lookup_failed = n_reference_lookup_failed,
    n_uncovered_positions = n_uncovered_positions,
    n_unusable_positions = n_unusable_positions,
    observed_fraction = mean(usable_observed_call),
    covered_fraction = covered_fraction,
    total_effect_alleles = total_effect_alleles,
    prs = prs
)

write_tsv(
    summary,
    opt$output_summary,
    na = "NA"
)

# Checked after both outputs are written so the QC record survives the
# failure and can be inspected.

if (
    min_covered_fraction > 0 &&
    covered_fraction < min_covered_fraction
) {
    stopf(
        paste0(
            "Sample '%s' has usable coverage %.4f across %d score ",
            "variants, below --min_covered_fraction=%.4f; ",
            "%d position(s) lack a usable SNP call or validated ",
            "hom-ref block (%d have no VCF record or block)."
        ),
        opt$sample_id,
        covered_fraction,
        nrow(details),
        min_covered_fraction,
        n_unusable_positions,
        n_uncovered_positions
    )
}
