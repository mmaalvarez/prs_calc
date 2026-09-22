#!/usr/bin/env bash

nextflow run "$PWD/../main.nf" \
    --input "$PWD/samplesheet.tsv" \
    --scorefile "$PWD/prs_calc_toy/PGS000392_hg19.txt" \
    --phenotypes "$PWD/prs_calc_toy/phenotypes.tsv" \
    --target_build hg19 \
    --missing_genotype reference \
    -profile conda \
    -resume
