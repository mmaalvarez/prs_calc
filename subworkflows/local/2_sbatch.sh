for vcf in `ls vcfs/*vcf`
do
	for prs_model in 'hung_etal_2021/PGS000740_hmPOS_GRCh37'
	do
		sbatch 2_run_calculate_prs.sh $vcf ""/g/strcombio/fsupek_data/users/malvarez/projects/lucia/data/published_PRS_models/""$prs_model"".txt""
	done
done

