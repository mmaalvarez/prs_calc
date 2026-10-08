python3 make_prs_calc_toy.py --out prs_calc_toy \
	  --command 'nextflow run ../../main.nf -w {work} --input {vcf} --scorefile {score} --outdir {out} --non_additive {non_additive} --genotype_calls {genotype_calls} --no_calls {no_calls}'
