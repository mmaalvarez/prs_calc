#!/usr/bin/env bash

set -euo pipefail

unset R_HOME
export NXF_HOME="$HOME/.nextflow"

nextflow run "$PWD/../main.nf" \
    --input "$PWD/prs_calc_toy/samplesheet.tsv" \
    --scorefile "$PWD/prs_calc_toy/PGS000392_hg19.txt" \
    --phenotypes "$PWD/prs_calc_toy/phenotypes.tsv" \
    --target_build hg19 \
    -profile conda \
    -resume
