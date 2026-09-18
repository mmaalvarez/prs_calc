#!/bin/bash
#SBATCH --partition="normal_prio"
#SBATCH --mem=4G
#SBATCH -c 1

snp_file="/g/strcombio/fsupek_data/users/malvarez/projects/lucia/data/published_PRS_models/hung_etal_2021/PGS000740_hmPOS_GRCh37.txt"
vcf_dir="/g/strcombio/fsupek_fisher/malvarez/gVCF/ICGC/TCGA_WGS/LUAD_LUSC_OV_PRAD_UCEC_UCS/single_sample_gVCFs"

# Build a regions file once
#grep -v "^#" "$snp_file" | awk 'NR>1 {print "chr"$2"\t"$3"\t"$3}' > regions.txt

mkdir vcfs

for vcf in "$vcf_dir"/*.gz
do
    sampleid=$(basename "$vcf" .gz)
    bcftools view -R regions.txt "$vcf" > vcfs/"${sampleid}_hung_2021_snps.vcf"
done
