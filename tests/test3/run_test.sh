#!/usr/bin/env bash

nextflow run "$PWD/../../main.nf" \
    --input "$PWD/prs_calc_toy/samplesheet" \
    --scorefile "$PWD/prs_calc_toy/PGS001229_22.txt" \
    --phenotypes "$PWD/prs_calc_toy/fake_phenotypes" \
    --target_build hg37 \
    --missing_genotype reference \
    --no_calls zero \
    --gvcf_mode plain \
    -profile conda \
    -resume
