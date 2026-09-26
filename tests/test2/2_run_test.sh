# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_reference.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference --strict_alleles true \
#   --default_ploidy 2 --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_reference

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_zero.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype zero --strict_alleles true \
#   --default_ploidy 2 --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_zero

# nextflow run ../../main.nf \
#   --input prs_calc_toy/gvcf_full.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode gvcf \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_gvcf_full

# # This specifically tests INFO/END + ALT=. auto-detection.
# # auto requires an explicit missing-genotype policy at the main.nf level.
# nextflow run ../../main.nf \
#   --input prs_calc_toy/gvcf_star.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode auto \
#   --missing_genotype error --strict_alleles true \
#   --default_ploidy 2 --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_gvcf_star

# nextflow run ../../main.nf \
#   --input prs_calc_toy/gvcf_gap.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode gvcf \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_gvcf_gap
  
# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_single.tsv \
#   --scorefile prs_calc_toy/score_pal_effect_swapped.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_score_pal_effect_swapped
  
# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_single.tsv \
#   --scorefile prs_calc_toy/score_duplicate_position.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_score_duplicate_position

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_mismatch.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_mismatch

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_tied.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_tied

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_ambiguous_contigs.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_ambiguous_contigs

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_symbolic.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_symbolic

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_symbolic_ref_effect.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_symbolic_ref_effect

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_two_samples.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_two_samples

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_no_contigs.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_no_contigs

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_numeric_ids.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_plain_numeric_ids

# nextflow run ../../main.nf \
#   --input prs_calc_toy/gvcf_bad_block.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode gvcf \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_gvcf_bad_block

# nextflow run ../../main.nf \
#   --input prs_calc_toy/gvcf_gap.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --target_build hg38 --gvcf_mode gvcf \
#   --missing_genotype error \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_gvcf_gap_missing_error

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_zero.tsv \
#   --scorefile prs_calc_toy/score_bad_position.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_score_bad_position

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_zero.tsv \
#   --scorefile prs_calc_toy/score_missing_weight.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_score_missing_weight

# nextflow run ../../main.nf \
#   --input prs_calc_toy/plain_zero.tsv \
#   --scorefile prs_calc_toy/score_non_snv.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_score_non_snv

# nextflow run ../../main.nf \
#   --input prs_calc_toy/cohort.tsv \
#   --scorefile prs_calc_toy/score.tsv \
#   --phenotypes prs_calc_toy/phenotypes.tsv \
#   --target_build hg38 --gvcf_mode plain \
#   --missing_genotype reference \
#   --strict_alleles true --default_ploidy 2 \
#   --plot_top_variants 12 \
#   --outdir prs_calc_toy/out_cohort

# nextflow run ../../main.nf \
#   --input "prs_calc_toy/plain_reference.tsv" \
#   --scorefile "prs_calc_toy/score_s08_only.tsv" \
#   --target_build hg38 \
#   --gvcf_mode plain \
#   --missing_genotype error \
#   --outdir "prs_calc_toy/other_tests/out_reg_no_call_error"

# nextflow run ../../main.nf \
#   --input "prs_calc_toy/plain_reference.tsv" \
#   --scorefile "prs_calc_toy/score_s08_only.tsv" \
#   --target_build hg38 \
#   --gvcf_mode plain \
#   --missing_genotype zero \
#   --outdir "prs_calc_toy/other_tests/out_reg_no_call_zero"

nextflow run ../../main.nf \
  --input "prs_calc_toy/plain_single.tsv" --scorefile "prs_calc_toy/score_s06_only.tsv" --gvcf_mode plain --missing_genotype error \
  --outdir "prs_calc_toy/other_tests/out_score_s06_only"

# nextflow run ../../main.nf \
#   --input "prs_calc_toy/gvcf_gap.tsv" --scorefile "prs_calc_toy/score.tsv" --gvcf_mode gvcf --min_covered_fraction 1 \
#   --outdir "prs_calc_toy/other_tests/out_min_covered_fraction"

# nextflow run ../../main.nf \
#   --input "prs_calc_toy/plain_single.tsv" --scorefile "prs_calc_toy/score_non_snv.tsv" --gvcf_mode plain --missing_genotype reference \
#   --outdir "prs_calc_toy/other_tests/out_score_non_snv_ref"

nextflow run ../../main.nf \
  --input "prs_calc_toy/plain_mismatch.tsv" --scorefile "prs_calc_toy/score.tsv" --gvcf_mode plain --strict_alleles false \
  --missing_genotype zero \
  --outdir "prs_calc_toy/other_tests/out_plain_mismatch2"

nextflow run ../../main.nf \
  --input "prs_calc_toy/plain_tied.tsv" --scorefile "prs_calc_toy/score.tsv" --gvcf_mode plain \
  --missing_genotype zero \
  --outdir "prs_calc_toy/other_tests/out_plain_tied2"

# nextflow run ../../main.nf \
#   --input "prs_calc_toy/gvcf_bad_block.tsv" --scorefile "prs_calc_toy/score.tsv" --gvcf_mode gvcf \
#   --outdir "prs_calc_toy/other_tests/out_bad_block2"
