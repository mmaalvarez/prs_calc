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
        "--sample-id",
        dest = "sample_id",
        type = "character"
    ),
    make_option(
        "--vcf-sample-file",
        dest = "vcf_sample_file",
        type = "character"
    ),
    make_option(
        "--source-vcf",
        dest = "source_vcf",
        type = "character"
    ),
    make_option(
        "--target-build",
        dest = "target_build",
        type = "character"
    ),
    make_option(
        "--missing-genotype",
        dest = "missing_genotype",
        type = "character",
        default = NA_character_,
        help = paste0(
            "Required unless the input is treated as a gVCF. How to ",
            "score variants with no usable genotype: 'reference' ",
            "(assume homozygous for the reference-genome base), 'zero' ",
            "(dosage 0) or 'error' (abort)."
        )
    ),
    make_option(
        "--gvcf-mode",
        dest = "gvcf_mode",
        type = "character",
        default = "auto",
        help = paste0(
            "Asserts how gVCF status was decided upstream: 'auto', ",
            "'gvcf' or 'plain'. Cross-checked against --input-is-gvcf ",
            "and recorded in the summary; 'gvcf' or 'plain' disagreeing ",
            "with --input-is-gvcf is an error. [default: %default]"
        )
    ),
    make_option(
        "--default-ploidy",
        dest = "default_ploidy",
        type = "integer",
        default = 2L
    ),
    make_option(
        "--strict-alleles",
        dest = "strict_alleles",
        type = "character",
        default = "true"
    ),
    make_option(
        "--input-is-gvcf",
        dest = "input_is_gvcf",
        type = "character",
        default = "false"
    ),
    make_option(
        "--min-covered-fraction",
        dest = "min_covered_fraction",
        type = "double",
        default = 0,
        help = paste0(
            "Abort if the fraction of score variants with either a VCF ",
            "record or an overlapping hom-ref reference block falls ",
            "below this value. 0 disables the check. [default: %default]"
        )
    ),
    make_option(
        "--output-summary",
        dest = "output_summary",
        type = "character"
    ),
    make_option(
        "--output-details",
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
    "--strict-alleles"
)

input_is_gvcf <- parse_boolean(
    opt$input_is_gvcf,
    "--input-is-gvcf"
)


gvcf_mode_declared <- tolower(trimws(as.character(opt$gvcf_mode)))

if (!gvcf_mode_declared %in% c("auto", "gvcf", "plain")) {
    stopf(
        "--gvcf-mode must be auto, gvcf or plain, not '%s'",
        gvcf_mode_declared
    )
}

min_covered_fraction <- suppressWarnings(
    as.numeric(opt$min_covered_fraction)
)

if (
    length(min_covered_fraction) != 1L ||
    is.na(min_covered_fraction) ||
    min_covered_fraction < 0 ||
    min_covered_fraction > 1
) {
    stopf("--min-covered-fraction must be a number between 0 and 1")
}

#
# --gvcf-mode records what the caller declared; --input-is-gvcf records
# what was resolved upstream. If they disagree, one of the two is wrong
# and the consequences are silent (a declared 'plain' input that is
# handled as a gVCF discards --missing-genotype without notice), so this
# is an error rather than a warning.
#
if (gvcf_mode_declared == "gvcf" && !input_is_gvcf) {
    stopf(
        paste0(
            "--gvcf-mode gvcf was declared but --input-is-gvcf is ",
            "false. Pass --input-is-gvcf true, or declare ",
            "--gvcf-mode auto and let the upstream detection decide."
        )
    )
}

if (gvcf_mode_declared == "plain" && input_is_gvcf) {
    stopf(
        paste0(
            "--gvcf-mode plain was declared but --input-is-gvcf is ",
            "true. Pass --input-is-gvcf false to force plain handling ",
            "(--missing-genotype then becomes mandatory), or declare ",
            "--gvcf-mode auto."
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
        "--missing-genotype must be reference, zero or error, not '%s'",
        requested_missing_mode
    )
}

if (input_is_gvcf) {
    #
    # Reference blocks assert hom-ref explicitly, and absence outside
    # every block means "not covered", so neither 'reference' nor 'zero'
    # is a decision the caller has to make. 'error' stays honoured
    # because aborting on uncovered positions is still meaningful.
    #
    missing_mode <- if (
        missing_supplied &&
        requested_missing_mode == "error"
    ) {
        "error"
    } else {
        "gvcf"
    }

    if (
        missing_supplied &&
        requested_missing_mode != "error"
    ) {
        ignored_notice <- sprintf(
            paste0(
                "--missing-genotype='%s' is ignored because the input ",
                "is treated as a gVCF: hom-ref blocks are honoured as ",
                "observed evidence, and positions outside every block ",
                "are scored as uncovered (dosage 0) rather than ",
                "assumed homozygous reference."
            ),
            requested_missing_mode
        )

        #
        # Always a warning, including under --gvcf-mode auto: a caller
        # that passes a flag which cannot apply should hear about it
        # even when the gVCF decision was made upstream.
        #
        warning(ignored_notice, call. = FALSE)
    }
} else {
    if (!missing_supplied) {
        stop(
            paste0(
                "--missing-genotype is required for a plain VCF/BCF. ",
                "Absence of a record in a sites-only file is ",
                "ambiguous, so the handling of unobserved variants ",
                "must be stated explicitly:\n",
                "  --missing-genotype=reference  assume homozygous for ",
                "the reference-genome base\n",
                "  --missing-genotype=zero       assign dosage 0 ",
                "(biased low, but never inverted)\n",
                "  --missing-genotype=error      abort on any missing ",
                "genotype\n",
                "'reference' is only defensible for a ",
                "reference-oriented file; see the orientation warning ",
                "in the run log."
            ),
            call. = FALSE
        )
    }

    missing_mode <- requested_missing_mode
}

default_ploidy <- as.integer(opt$default_ploidy)

if (
    is.na(default_ploidy) ||
    default_ploidy < 1L
) {
    stopf("--default-ploidy must be a positive integer")
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

allele_is_reference <- function(
    allele,
    chromosome,
    position
) {
    allele <- normalise_allele(allele)

    if (
        length(allele) != 1L ||
        is.na(allele) ||
        !nzchar(allele) ||
        !grepl("^[ACGTN]+$", allele)
    ) {
        return(FALSE)
    }

    observed_reference <- get_reference_sequence(
        chromosome,
        position,
        nchar(allele)
    )

    !is.na(observed_reference) &&
        identical(
            toupper(allele),
            observed_reference
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

find_allele_index <- function(ref, alt, allele) {
    allele <- normalise_allele(allele)

    if (
        length(allele) != 1L ||
        is.na(allele)
    ) {
        return(NA_integer_)
    }

    alleles <- record_alleles(ref, alt)
    index <- match(allele, alleles)

    if (is.na(index)) {
        return(NA_integer_)
    }

    # VCF genotype indices are zero-based:
    # REF = 0, first ALT = 1, second ALT = 2, ...
    as.integer(index - 1L)
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
                "--missing-genotype=error"
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
        return(list(
            dosage = NA_real_,
            status = paste0(
                reason,
                "_unscored_no_reference"
            ),
            dosage_from_reference = FALSE
        ))
    }

    if (model_mismatch) {
        return(list(
            dosage = NA_real_,
            status = paste0(
                reason,
                "_unscored_reference_mismatch"
            ),
            dosage_from_reference = FALSE
        ))
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

block_index <- list()
genotype_index <- list()

if (nrow(genotypes) > 0L) {
    genotypes <- genotypes %>%
        mutate(
            chrom = normalise_chromosome(chrom),
            ref = normalise_allele(ref),
            alt = normalise_allele(alt),
            key = paste(chrom, position, sep = ":"),
            is_reference_only = reference_only_alt(alt)
        )

    #
    # Reference blocks must assert a homozygous reference call. A
    # reference-only record with a no-call genotype carries no
    # information, so the position stays uncovered.
    #
    genotype_primary <- sub(
        ":.*$",
        "",
        genotypes$genotype
    )

    is_hom_ref_call <- !is.na(genotype_primary) &
        grepl("^0([/|]0)*$", genotype_primary)

    block_ploidy <- ifelse(
        is_hom_ref_call,
        lengths(
            strsplit(
                ifelse(
                    is.na(genotype_primary),
                    "",
                    genotype_primary
                ),
                "[/|]"
            )
        ),
        NA_integer_
    )

    is_block <- genotypes$is_reference_only & is_hom_ref_call

    if (any(is_block)) {
        block_ref_width <- ifelse(
            is.na(genotypes$ref),
            1L,
            nchar(genotypes$ref)
        )

        block_end <- ifelse(
            !is.na(genotypes$end),
            genotypes$end,
            genotypes$position + block_ref_width - 1L
        )

        reference_blocks <- tibble(
            chrom = genotypes$chrom[is_block],
            start = genotypes$position[is_block],
            end = as.integer(
                pmax(
                    block_end[is_block],
                    genotypes$position[is_block]
                )
            ),
            ref = genotypes$ref[is_block],
            ploidy = as.integer(block_ploidy[is_block])
        ) %>%
            arrange(chrom, start)

        block_index <- split(
            seq_len(nrow(reference_blocks)),
            reference_blocks$chrom
        )
    }

    variant_rows <- which(!genotypes$is_reference_only)

    if (length(variant_rows) > 0L) {
        genotype_index <- split(
            variant_rows,
            genotypes$key[variant_rows]
        )
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
        # Only comparable to the genome when the block starts here.
        ref = if (
            reference_blocks$start[[block_row]] == position
        ) {
            reference_blocks$ref[[block_row]]
        } else {
            NA_character_
        }
    )
}

#
# Coverage policy. A hom-ref block is observed evidence and is honoured
# regardless of --missing-genotype. Only genuinely uncovered positions
# consult the policy, and in a gVCF absence means "not assessed", so
# those are set to zero rather than assumed homozygous reference.
#
resolve_missing <- function(
    reason,
    variant_id,
    effect_is_reference,
    model_mismatch,
    reference_base,
    block_row
) {
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

    if (input_is_gvcf) {
        return(
            fallback_genotype(
                reason = paste0(reason, "_uncovered"),
                variant_id = variant_id,
                effect_is_reference = effect_is_reference,
                model_mismatch = model_mismatch,
                reference_base = reference_base,
                policy = if (missing_mode == "error") {
                    "error"
                } else {
                    "zero"
                }
            )
        )
    }

    fallback_genotype(
        reason = reason,
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

    reference_base <- get_reference_sequence(
        chromosome,
        position,
        1L
    )

    effect_is_reference <- allele_is_reference(
        effect_allele,
        chromosome,
        position
    )

    other_is_reference <- allele_is_reference(
        other_allele,
        chromosome,
        position
    )

    #
    # If the other allele is unknown, a non-reference effect allele is
    # still valid. If both alleles are known, one of them should agree
    # with the reference genome.
    #
    model_mismatch <- (
        !is.na(other_allele) &&
        !effect_is_reference &&
        !other_is_reference
    )

    record_indices <- genotype_index[[key]]

    record_found <- (
        !is.null(record_indices) &&
        length(record_indices) > 0L
    )

    selected_record <- NULL
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
    effect_allele_index <- NA_integer_
    other_allele_index <- NA_integer_

    reference_block_found <- FALSE
    reference_block_text <- NA_character_
    reference_block_ref <- NA_character_

    block_row <- find_reference_block(
        chromosome,
        position
    )

    if (record_found) {
        records <- genotypes[record_indices, ]

        compatibility_scores <- vapply(
            seq_len(nrow(records)),
            function(record_number) {
                record <- records[record_number, ]
                alleles <- record_alleles(
                    record$ref[[1]],
                    record$alt[[1]]
                )

                score <- 0L

                if (effect_allele %in% alleles) {
                    score <- score + 4L
                }

                if (
                    !is.na(other_allele) &&
                    other_allele %in% alleles
                ) {
                    score <- score + 2L
                }

                if (
                    allele_is_reference(
                        record$ref[[1]],
                        chromosome,
                        position
                    )
                ) {
                    score <- score + 1L
                }

                score
            },
            integer(1)
        )

        best_score <- max(compatibility_scores)
        best_records <- which(
            compatibility_scores == best_score
        )

        record_ambiguous <- length(best_records) > 1L
        selected_record <- records[best_records[[1]], ]

        vcf_ref <- normalise_allele(
            selected_record$ref[[1]]
        )

        vcf_alt <- normalise_allele(
            selected_record$alt[[1]]
        )

        genotype <- selected_record$genotype[[1]]

        alleles <- record_alleles(
            vcf_ref,
            vcf_alt
        )

        effect_allele_index <- find_allele_index(
            vcf_ref,
            vcf_alt,
            effect_allele
        )

        if (!is.na(other_allele)) {
            other_allele_index <- find_allele_index(
                vcf_ref,
                vcf_alt,
                other_allele
            )
        }

        # The score alleles may appear in either VCF orientation:
        #
        #   REF=other,  ALT=effect
        #   REF=effect, ALT=other
        #
        # The effect allele must occur in the VCF allele list. When an
        # other allele is supplied, it must also occur and must be
        # distinct from the effect allele.
        allele_mismatch <- (
            is.na(effect_allele_index) ||
            (
                !is.na(other_allele) &&
                (
                    is.na(other_allele_index) ||
                    other_allele_index ==
                        effect_allele_index
                )
            )
        )

        genotype_tokens <- if (
            is.na(genotype) ||
            !nzchar(genotype)
        ) {
            character(0)
        } else {
            strsplit(
                sub(":.*$", "", genotype),
                "[/|]",
                perl = TRUE
            )[[1]]
        }

        callable <- (
            length(genotype_tokens) > 0L &&
            all(grepl("^[0-9]+$", genotype_tokens))
        )

        genotype_indices <- suppressWarnings(
            as.integer(genotype_tokens)
        )

        callable <- (
            callable &&
            all(!is.na(genotype_indices)) &&
            all(genotype_indices >= 0L) &&
            all(genotype_indices < length(alleles))
        )

        if (callable) {
            called_alleles <- alleles[
                genotype_indices + 1L
            ]

            called_alleles_text <- paste(
                called_alleles,
                collapse = "/"
            )

            #
            # A record whose ALT is symbolic only (e.g. <DEL>) cannot
            # contain the effect allele, but a hom-ref call on it is
            # still positive evidence of the reference genotype, so the
            # dosage comes from the genome rather than the allele list.
            # <NON_REF>/<*> records never reach here: they are
            # classified as reference-only and become blocks instead.
            #
            hom_ref_symbolic <- allele_mismatch &&
                is.na(effect_allele_index) &&
                symbolic_only_alt(vcf_alt) &&
                all(genotype_indices == 0L)

            if (hom_ref_symbolic) {
                fallback <- fallback_genotype(
                    reason = "hom_ref_symbolic_alt",
                    variant_id = variant_id,
                    effect_is_reference = effect_is_reference,
                    model_mismatch = model_mismatch,
                    reference_base = reference_base,
                    policy = "reference",
                    ploidy = length(genotype_indices)
                )

                dosage <- fallback$dosage
                status <- fallback$status
                dosage_from_reference <- fallback$dosage_from_reference
                allele_mismatch <- FALSE
            } else {
                # Count copies of the effect allele by its VCF allele
                # index. This works whether the effect allele is REF,
                # ALT, or one allele in a multiallelic record.
                dosage <- if (
                    is.na(effect_allele_index)
                ) {
                    NA_real_
                } else {
                    as.numeric(
                        sum(
                            genotype_indices ==
                                effect_allele_index
                        )
                    )
                }

                model_alleles <- unique(
                    na.omit(
                        c(
                            effect_allele,
                            other_allele
                        )
                    )
                )

                non_model_call <- any(
                    !called_alleles %in% model_alleles
                )

                status <- if (allele_mismatch) {
                    "observed_allele_mismatch"
                } else if (record_ambiguous) {
                    "observed_ambiguous_record"
                } else if (non_model_call) {
                    "observed_nonmodel_allele"
                } else {
                    "observed"
                }
            }
        } else {
            no_call <- TRUE

            fallback <- resolve_missing(
                reason = "no_call",
                variant_id = variant_id,
                effect_is_reference = effect_is_reference,
                model_mismatch = model_mismatch,
                reference_base = reference_base,
                block_row = block_row
            )

            dosage <- fallback$dosage
            status <- fallback$status
            dosage_from_reference <- fallback$dosage_from_reference

            block_details <- describe_reference_block(
                block_row,
                position
            )

            reference_block_found <- block_details$found
            reference_block_text <- block_details$text
            reference_block_ref <- block_details$ref
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
        dosage_from_reference <- fallback$dosage_from_reference

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
        effect_allele_is_reference = effect_is_reference,
        other_allele_is_reference = other_is_reference,
        vcf_record_found = record_found,
        vcf_record_ambiguous = record_ambiguous,
        reference_block_found = reference_block_found,
        reference_block = reference_block_text,
        reference_block_ref = reference_block_ref,
        vcf_ref = vcf_ref,
        vcf_alt = vcf_alt,
        genotype = genotype,
        called_alleles = called_alleles_text,
        vcf_no_call = no_call,
        non_model_allele_call = non_model_call,
        dosage_from_reference = dosage_from_reference,
        allele_mismatch = allele_mismatch,
        effect_allele_dosage = dosage,
        contribution = contribution,
        status = status
    )
}

details <- bind_rows(detail_rows) %>%
    arrange(score_row)

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

if (n_reference_lookup_failed == nrow(details)) {
    stopf(
        paste0(
            "No score position could be resolved in %s. Check that ",
            "--target-build matches the score file and that chromosome ",
            "names are recognised (seen: %s)."
        ),
        genome_package,
        paste(head(unique(details$chrom), 5L), collapse = ", ")
    )
} else if (n_reference_lookup_failed > 0L) {
    warning(
        sprintf(
            "Sample '%s': %d of %d score positions had no %s reference base; orientation checks skip those rows.",
            opt$sample_id, n_reference_lookup_failed, nrow(details), target_build
        ),
        call. = FALSE
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
        !details$reference_block_found &
        !grepl("^hom_ref_symbolic_alt", details$status)
)

covered_fraction <- mean(
    details$vcf_record_found |
        details$reference_block_found
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
                "really a sites-only VCF, use --input-is-gvcf false ",
                "with an explicit --missing-genotype instead."
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
                "--missing-genotype=%s. If this is an all-sites or ",
                "gVCF input, use --input-is-gvcf true."
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
            "against, so this indicates the wrong --target-build ",
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
                    "--missing-genotype=zero, or re-orient the VCF ",
                    "with 'bcftools norm --check-ref s'."
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

#
# A gVCF's REF fields come from the FASTA it was called against, so a
# block REF that disagrees with BSgenome indicates the wrong
# --target-build rather than a REF/ALT swap. Only checkable when a block
# happens to start at a score position.
#
n_block_ref_not_genome_ref <- sum(
    !is.na(details$reference_block_ref) &
        !is.na(details$reference_allele) &
        substr(details$reference_block_ref, 1L, 1L) !=
            details$reference_allele,
    na.rm = TRUE
)

if (n_block_ref_not_genome_ref > 0L) {
    warning(
        sprintf(
            paste0(
                "Sample '%s': %d reference block(s) begin at a score ",
                "position with a REF allele that disagrees with the ",
                "%s reference base. The input was probably called ",
                "against a different reference build, so ",
                "block-derived dosages are unreliable."
            ),
            opt$sample_id,
            n_block_ref_not_genome_ref,
            target_build
        ),
        call. = FALSE
    )
}

#
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

n_scored <- sum(
    !is.na(details$contribution)
)

prs <- if (n_scored == 0L) {
    NA_real_
} else {
    sum(
        details$contribution,
        na.rm = TRUE
    )
}

total_effect_alleles <- if (
    all(is.na(details$effect_allele_dosage))
) {
    NA_real_
} else {
    sum(
        details$effect_allele_dosage,
        na.rm = TRUE
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
    min_covered_fraction = min_covered_fraction,
    n_homozygous_reference_blocks = n_homozygous_reference_blocks,
    n_positions_from_reference_block = n_positions_from_reference_block,
    n_reference_lookup_failed = n_reference_lookup_failed,
    n_uncovered_positions = n_uncovered_positions,
    observed_fraction = mean(details$vcf_record_found),
    covered_fraction = covered_fraction,
    total_effect_alleles = total_effect_alleles,
    prs = prs
)

write_tsv(
    summary,
    opt$output_summary,
    na = "NA"
)

if (
    strict_alleles &&
    any(
        details$allele_mismatch,
        na.rm = TRUE
    )
) {
    bad_variants <- unique(
        details$variant_id[
            details$allele_mismatch %in% TRUE
        ]
    )

    stopf(
        paste0(
            "Score/VCF allele mismatch for sample '%s': %s. ",
            "Check --target_build, the score position and allele ",
            "columns, and the VCF allele definitions. To retain ",
            "mismatched rows as QC results, use ",
            "--strict-alleles false."
        ),
        opt$sample_id,
        paste(head(bad_variants, 10L), collapse = ", ")
    )
}

if (
    strict_alleles &&
    any(is.na(details$contribution))
) {
    unscored_variants <- details$variant_id[
        is.na(details$contribution)
    ]

    stopf(
        "Sample '%s' contains unscored variant(s): %s",
        opt$sample_id,
        paste(
            head(unique(unscored_variants), 10L),
            collapse = ", "
        )
    )
}

#
# Checked after both outputs are written so the QC record survives the
# failure and can be inspected.
#
if (
    min_covered_fraction > 0 &&
    covered_fraction < min_covered_fraction
) {
    stopf(
        paste0(
            "Sample '%s' covers %.4f of %d score variants, below ",
            "--min-covered-fraction=%.4f. %d position(s) had neither a ",
            "VCF record nor an overlapping hom-ref reference block."
        ),
        opt$sample_id,
        covered_fraction,
        nrow(details),
        min_covered_fraction,
        n_uncovered_positions
    )
}
