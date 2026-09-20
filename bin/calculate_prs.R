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
        default = "reference"
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

missing_mode <- tolower(
    trimws(opt$missing_genotype)
)

if (!missing_mode %in% c("reference", "zero", "error")) {
    stopf(
        "--missing-genotype must be reference, zero or error"
    )
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

get_reference_sequence <- function(
    chromosome,
    position,
    width = 1L
) {
    seqname <- genome_seqname(chromosome)

    if (is.na(seqname)) {
        return(NA_character_)
    }

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
    reference_base
) {
    if (missing_mode == "error") {
        stopf(
            paste0(
                "Variant '%s' has a %s genotype/site and ",
                "--missing-genotype=error"
            ),
            variant_id,
            reason
        )
    }

    if (missing_mode == "zero") {
        return(list(
            dosage = 0,
            status = paste0(reason, "_set_to_zero"),
            assumed_reference = FALSE
        ))
    }

    if (is.na(reference_base)) {
        return(list(
            dosage = NA_real_,
            status = paste0(
                reason,
                "_unscored_no_reference"
            ),
            assumed_reference = FALSE
        ))
    }

    if (model_mismatch) {
        return(list(
            dosage = NA_real_,
            status = paste0(
                reason,
                "_unscored_reference_mismatch"
            ),
            assumed_reference = FALSE
        ))
    }

    dosage <- if (effect_is_reference) {
        default_ploidy
    } else {
        0
    }

    list(
        dosage = as.numeric(dosage),
        status = paste0(
            reason,
            "_assumed_reference"
        ),
        assumed_reference = TRUE
    )
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

genotypes <- read_tsv(
    opt$genotypes,
    na = c("", "NA"),
    show_col_types = FALSE,
    progress = FALSE,
    col_types = cols(
        chrom = col_character(),
        position = col_integer(),
        ref = col_character(),
        alt = col_character(),
        genotype = col_character()
    )
)

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

if (nrow(genotypes) > 0L) {
    genotypes <- genotypes %>%
        mutate(
            chrom = normalise_chromosome(chrom),
            ref = normalise_allele(ref),
            alt = normalise_allele(alt),
            key = paste(chrom, position, sep = ":")
        )

    genotype_index <- split(
        seq_len(nrow(genotypes)),
        genotypes$key
    )
} else {
    genotype_index <- list()
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

    effect_is_reference <- (
        !is.na(reference_base) &&
        effect_allele == reference_base
    )

    other_is_reference <- (
        !is.na(other_allele) &&
        !is.na(reference_base) &&
        other_allele == reference_base
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
    assumed_reference <- FALSE
    allele_mismatch <- FALSE
    effect_allele_index <- NA_integer_
    other_allele_index <- NA_integer_

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

            if (allele_mismatch) {
                status <- "observed_allele_mismatch"
            } else if (record_ambiguous) {
                status <- "observed_ambiguous_record"
            } else if (non_model_call) {
                status <- "observed_nonmodel_allele"
            } else {
                status <- "observed"
            }
        } else {
            no_call <- TRUE

            fallback <- fallback_genotype(
                reason = "no_call",
                variant_id = variant_id,
                effect_is_reference = effect_is_reference,
                model_mismatch = model_mismatch,
                reference_base = reference_base
            )

            dosage <- fallback$dosage
            status <- fallback$status
            assumed_reference <- fallback$assumed_reference
        }
    } else {
        fallback <- fallback_genotype(
            reason = "absent",
            variant_id = variant_id,
            effect_is_reference = effect_is_reference,
            model_mismatch = model_mismatch,
            reference_base = reference_base
        )

        dosage <- fallback$dosage
        status <- fallback$status
        assumed_reference <- fallback$assumed_reference
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
        vcf_ref = vcf_ref,
        vcf_alt = vcf_alt,
        genotype = genotype,
        called_alleles = called_alleles_text,
        vcf_no_call = no_call,
        non_model_allele_call = non_model_call,
        assumed_reference = assumed_reference,
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
    n_vcf_records_found = sum(
        details$vcf_record_found
    ),
    n_absent_from_vcf = sum(
        !details$vcf_record_found
    ),
    n_vcf_no_calls = sum(
        details$vcf_no_call
    ),
    n_assumed_reference = sum(
        details$assumed_reference
    ),
    n_allele_mismatches = sum(
        details$allele_mismatch,
        na.rm = TRUE
    ),
    n_non_model_allele_calls = sum(
        details$non_model_allele_call
    ),
    observed_fraction = mean(
        details$vcf_record_found
    ),
    total_effect_alleles = total_effect_alleles,
    prs = prs
)

write_tsv(
    summary,
    opt$output_summary,
    na = "NA"
)
