python3 make_prs_calc_toy.py --cohort --out prs_calc_toy
# Add --cohort only if you want the optional 10-sample plotting smoke test.

# add other tests
for id in s01 s06 s07 s08; do
  awk -F '\t' -v id="$id" 'NR == 1 || $1 == id' \
    "prs_calc_toy/score.tsv" > "prs_calc_toy/score_${id}_only.tsv"
done

## If there were multisample vcfs, compress and index the files after generating them; bgzip and tabix are required for these commands:
#for vcf in prs_calc_toy/plain/*.vcf prs_calc_toy/gvcf/*.vcf
#do
#    bgzip -f "$vcf"
#    tabix -f -p vcf "$vcf.gz"
#done
