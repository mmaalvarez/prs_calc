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

# This recognises SNP records, ALT=. reference/no-call records, and
# symbolic gVCF records. A concrete multi-base REF or ALT makes the
# entire record unsuitable as an SNP genotype.
is_snv_or_symbolic_record <- function(ref, alt) {
    ref <- normalise_allele(ref)
    alt <- normalise_allele(alt)

    if (is.na(ref) || !grepl("^[ACGT]$", ref)) {
        return(FALSE)
    }

    if (is.na(alt)) {
        return(TRUE)
    }

    parts <- strsplit(alt, ",", fixed = TRUE)[[1]]

    all(
        grepl("^[ACGT]$", parts) |
            grepl("^<[^<>]+>$", parts) |
            parts == "*"
    )
}

# A partial GT such as 0/. is invalid.
# Only a wholly missing GT is handled by --no_calls.
genotype_has_missing <- function(gt) {
    if (is.na(gt) || !nzchar(gt)) {
        return(TRUE)
    }

    primary <- sub(":.*$", "", gt)
    tokens <- strsplit(primary, "[/|]", perl = TRUE)[[1]]

    length(tokens) == 0L || any(tokens == ".")
}

is_full_no_call <- function(gt) {
    if (is.na(gt) || !nzchar(gt)) {
        return(FALSE)
    }

    primary <- sub(":.*$", "", gt)

    # Includes ./., .|., and wholly missing polyploid GTs.
    grepl(
        "^\\.([/|]\\.)+$",
        primary
    )
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
    genotype_cols$end <- col_integer()
}

genotypes <- read_tsv(
    opt$genotypes,
    na = c("", "NA", "."),
    show_col_types = FALSE,
    progress = FALSE,
    col_types = do.call(cols, genotype_cols)
)

if (!"end" %in% names(genotypes)) {
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

scores <- scores %>%
    mutate(
        chrom = normalise_chromosome(chrom),
        effect_allele = normalise_allele(effect_allele),
        other_allele = normalise_allele(other_allele),
        key = paste(chrom, position, sep = ":")
    )

reference_blocks <- tibble(
    chrom = character(0),
    start = integer(0),
    end = integer(0),
    ref = character(0),
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
non_snv_index <- list()

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
                USE.NAMES = FALSE
            ),
            full_no_call = vapply(
                genotype,
                is_full_no_call,
                logical(1)
            )
        )

    n_ignored_non_snv <- sum(
        !genotypes$is_snv_compatible
    )

    if (n_ignored_non_snv > 0L) {

        non_snv_rows <- which(
            !genotypes$is_snv_compatible
        )

        non_snv_index <- split(
            non_snv_rows,
            genotypes$key[non_snv_rows]
        )

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
            ploidy = as.integer(block_ploidy)
        ) %>%
            distinct(
                chrom,
                start,
                end,
                ref,
                ploidy,
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

n_homozygous_reference_blocks <- nrow(reference_blocks)


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
        return(
            fallback_genotype(
                reason = paste0(reason, "_reference_block"),
                variant_id = variant_id,
                effect_is_reference = effect_is_reference,
                model_mismatch = model_mismatch,
                reference_base = reference_base,
                policy = "reference",
                ploidy = reference_blocks$ploidy[[block_row]]
            )
        )
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

    if (length(non_snv_index[[key]]) > 0L) {
        stopf(
            paste0(
                "Score SNP '%s' at %s:%d has a non-SNV VCF record ",
                "at that position; it cannot be treated as an SNP ",
                "no-call or as an absent site"
            ),
            variant_id,
            chromosome,
            position
        )
    }

    record_indices <- genotype_index[[key]]

    interval_row <- find_no_call_interval(
        chromosome,
        position
    )

    if (
        length(record_indices) == 0L &&
        !is.na(interval_row)
    ) {
        record_indices <- no_call_intervals$source_row[[interval_row]]
    }

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

    reference_block_found <- FALSE
    reference_block_text <- NA_character_
    reference_block_ref <- NA_character_

    block_row <- find_reference_block(
        chromosome,
        position
    )

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

            # An interval's REF is at its start, not necessarily at
            # this score position; only validate alleles when this
            # record starts at the score position.
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
            if (genotype_has_missing(genotype)) {
                stopf(
                    paste0(
                        "Score SNP '%s' at %s:%d has missing or ",
                        "partially missing GT '%s'; only a full ",
                        "GT=./. no-call is handled by --no_calls"
                    ),
                    variant_id,
                    chromosome,
                    position,
                    ifelse(is.na(genotype), "NA", genotype)
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

            primary_gt <- sub(
                ":.*$",
                "",
                genotype
            )

            valid_gt <- grepl(
                "^[0-9]+([/|][0-9]+)*$",
                primary_gt
            )

            genotype_indices <- if (valid_gt) {
                suppressWarnings(
                    as.integer(
                        strsplit(
                            primary_gt,
                            "[/|]",
                            perl = TRUE
                        )[[1]]
                    )
                )
            } else {
                NA_integer_
            }

            callable <- (
                valid_gt &&
                    all(!is.na(genotype_indices)) &&
                    all(genotype_indices >= 0L) &&
                    all(genotype_indices < length(alleles))
            )

            if (!callable) {
                stopf(
                    "Invalid SNP genotype '%s' for variant '%s' at %s:%d",
                    genotype,
                    variant_id,
                    chromosome,
                    position
                )
            }

            called_alleles <- alleles[
                genotype_indices + 1L
            ]

            if (
                any(
                    is.na(called_alleles) |
                        !grepl(
                            "^[ACGT]$",
                            called_alleles
                        )
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

            called_alleles_text <- paste(
                called_alleles,
                collapse = "/"
            )

            non_model_call <- (
                length(other_set) == 0L &&
                    any(
                        !(called_alleles %in% effect_set)
                    )
            )

            dosage <- as.numeric(
                sum(
                    called_alleles %in% effect_set
                )
            )

            status <- if (non_model_call) {
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
            block_row,
            position
        )

        reference_block_found <- block_details$found
        reference_block_text <- block_details$text
        reference_block_ref <- block_details$ref
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
        vcf_record_found = record_found,
        vcf_record_ambiguous = record_ambiguous,
        reference_block_found =
            reference_block_found,
        reference_block = reference_block_text,
        reference_block_ref = reference_block_ref,
        vcf_ref = vcf_ref,
        vcf_alt = vcf_alt,
        genotype = genotype,
        called_alleles = called_alleles_text,
        vcf_no_call = no_call,
        non_model_allele_call = non_model_call,
        dosage_from_reference =
            dosage_from_reference,
        allele_mismatch = allele_mismatch,
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

n_positions_from_reference_block <- sum(details$reference_block_found)

n_uncovered_positions <- sum(
    !details$vcf_record_found &
        !details$reference_block_found
)

n_assumed_reference <- sum(
    details$dosage_from_reference &
        !details$reference_block_found
)

usable_observed_call <- details$status %in% c(
    "observed",
    "observed_nonmodel_allele"
)

usable_reference_block <- (
    details$reference_block_found &
        details$status ==
        "absent_reference_block_dosage_from_reference"
)

usable_coverage <- (
    usable_observed_call |
        usable_reference_block
)

n_unusable_positions <- sum(!usable_coverage)

covered_fraction <- mean(
    usable_coverage
)


#
# gVCF declaration versus what the records actually look like. Advisory
# only: a gVCF whose blocks all happen to miss the score positions is
# legal, as is a plain all-sites VCF that carries hom-ref blocks.
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

#
# Reference orientation check. Dosage for called genotypes is counted by
# allele index and is therefore unaffected by a REF/ALT swap, but the
# 'reference' fallback compares the effect allele against the reference
# genome, so a non-reference-oriented VCF makes assumed genotypes
# unreliable.
#
n_vcf_ref_not_genome_ref <- sum(
    details$vcf_record_found &
        !is.na(details$vcf_ref) &
        !is.na(details$reference_allele) &
        substr(details$vcf_ref, 1L, 1L) !=
            details$reference_allele,
    na.rm = TRUE
)

n_palindromic_not_oriented <- sum(
    details$vcf_record_found &
        !is.na(details$vcf_ref) &
        !is.na(details$reference_allele) &
        substr(details$vcf_ref, 1L, 1L) != details$reference_allele &
        paste0(details$effect_allele, details$other_allele) %in%
            c("AT", "TA", "CG", "GC"),
    na.rm = TRUE
)

palindromic_note <- if (n_palindromic_not_oriented > 0L) {
    sprintf(
        paste0(
            " %d of them are palindromic (A/T or C/G), where a swap ",
            "cannot be distinguished from a strand flip, so those ",
            "dosages may be inverted."
        ),
        n_palindromic_not_oriented
    )
} else {
    paste0(
        " None are palindromic, so all are distinguishable REF/ALT ",
        "swaps and called dosages are correct."
    )
}

if (n_vcf_ref_not_genome_ref > 0L) {
    orientation_message <- sprintf(
        paste0(
            "Sample '%s': %d of %d VCF records have a REF allele that ",
            "disagrees with the %s reference base."
        ),
        opt$sample_id,
        n_vcf_ref_not_genome_ref,
        n_vcf_records_found,
        target_build
    )

    orientation_message <- if (input_is_gvcf) {
        paste0(
            orientation_message,
            " A gVCF's REF fields come from the FASTA it was called ",
            "against, so this indicates the wrong --target_build ",
            "rather than a REF/ALT swap, and all dosages are suspect."
        )
    } else if (
        missing_mode == "reference" &&
        n_assumed_reference > 0L
    ) {
        paste0(
            orientation_message,
            sprintf(
                paste0(
                    " This VCF is not reference-oriented, and %d ",
                    "variant(s) were assumed homozygous reference, so ",
                    "those dosages are unreliable. Consider ",
                    "--missing_genotype=zero, or re-orient the VCF ",
                    "with 'bcftools norm --check_ref s'."
                ),
                n_assumed_reference
            ),
            palindromic_note
        )
    } else {
        paste0(
            orientation_message,
            " This VCF is not reference-oriented. Called genotypes are ",
            "scored by allele index, so simple REF/ALT swaps are ",
            "handled correctly.",
            palindromic_note
        )
    }

    warning(orientation_message, call. = FALSE)
}


# Neither model allele matching the reference strand at a biallelic
# non-palindromic site means the score file and the VCF are both
# flipped relative to the genome. Called genotypes still score
# correctly, but the same file's palindromic sites will be flipped too
# and there no detector exists.
#
n_model_strand_mismatch <- sum(
    details$vcf_record_found &
        !details$allele_mismatch &
        !is.na(details$other_allele) &
        !is.na(details$reference_allele) &
        !details$effect_allele_is_reference &
        !details$other_allele_is_reference,
    na.rm = TRUE
)

if (n_model_strand_mismatch > 0L) {
    warning(
        sprintf(
            paste0(
                "Sample '%s': %d variant(s) have neither model allele ",
                "matching the %s reference base, which is the ",
                "signature of a score file and VCF that are jointly ",
                "strand-flipped. Called genotypes still score ",
                "correctly, but palindromic sites in the same file ",
                "cannot be checked and may be inverted."
            ),
            opt$sample_id,
            n_model_strand_mismatch,
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
    n_dosage_from_reference = n_assumed_reference,
    n_allele_mismatches = sum(
        details$allele_mismatch,
        na.rm = TRUE
    ),
    n_non_model_allele_calls = sum(
        details$non_model_allele_call
    ),
    n_vcf_ref_not_genome_ref = n_vcf_ref_not_genome_ref,
    n_model_strand_mismatch = n_model_strand_mismatch,
    n_palindromic_not_oriented = n_palindromic_not_oriented,
    n_block_ref_not_genome_ref = n_block_ref_not_genome_ref,
    input_is_gvcf = input_is_gvcf,
    gvcf_mode_declared = gvcf_mode_declared,
    missing_genotype_requested = requested_missing_mode,
    missing_genotype_mode = missing_mode,
    no_calls_mode = no_calls_mode,
    min_covered_fraction = min_covered_fraction,
    n_homozygous_reference_blocks = n_homozygous_reference_blocks,
    n_positions_from_reference_block = n_positions_from_reference_block,
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
