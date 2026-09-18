#!/bin/bash
#SBATCH --mem=2G
#SBATCH --time=5:00
#SBATCH -c 1
#SBATCH --partition=normal_prio

conda activate R

vcf=$1
prs_model=$2

Rscript 2_calculate_prs.R $vcf $prs_model
