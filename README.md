Miguel Martín Álvarez, PhD

miguel.m.alvarez3[--at--]gmail[--dot--]com

===========================================

- Example script to run the pipeline:

```
#!/usr/bin/env bash

set -euo pipefail

unset R_HOME
export NXF_HOME="$HOME/.nextflow"

nextflow run mmaalvarez/prs_calc -r main -latest \
	--input /path/to/samplesheet.tsv \
	--scorefile /path/to/PGSXXXXX_hg38.txt \
	--phenotypes /path/to/phenotypes.tsv \
	--target_build hg38 \
	-profile conda \
	-resume
```

To create the conda environment with the required packages, first run:
```
wget https://github.com/mmaalvarez/prs_calc/blob/main/envs/prs_calc.yml

conda env create -f prs_calc.yml -n prs_calc
```

For faster Conda environment creation, if mamba is installed:

`-profile mamba`

For Slurm plus Conda:

`-profile conda,slurm`


- Input formats

A one-column file fed to `--input` works directly:
```
/my/path/sample1.vcf.gz
/my/path/sample2.vcf.gz
/my/path/sample3.vcf.gz
/my/path/sample4.vcf.gz
```
Sample IDs become sample1, sample2, etc.

An optional explicit sample ID format is also supported:
```
patient_001	/my/path/sample1.vcf.gz
patient_002	/my/path/sample2.vcf.gz
```
Each VCF must contain exactly one sample.


An example input for `--scorefile` looks like this:
```
#HEADER
#...
#HEADER
rsID	chr_name	chr_position_hg38	effect_allele	other_allele	effect_weight	allelefrequency_effect	locus_name	variant_description	OR	hm_source	hm_rsID	hm_chr	hm_pos	hm_inferOtherAllele
rs71658797	1	77501822	A	T	0.131028262406404	0.1	AK5	35 SNP score	1.14	ENSEMBL	rs71658797	1	77501822	
rs13080835	3	189639410	T	G	-0.0618754037180875	0.49	TP63	35 SNP score	0.94	ENSEMBL	rs13080835	3	189639410	
rs7705526	5	1285859	A	C	0.113328685307003	0.34	TERT	35 SNP score	1.12	ENSEMBL	rs7705526	5	1285859	
rs112290073	5	1285917	A	G	0.3293037471426	0.01	TERT	35 SNP score	1.39	ENSEMBL	rs112290073	5	1285917	
rs2736098	5	1293971	T	C	0.131028262406404	0.28	TERT	35 SNP score	1.14	ENSEMBL	rs2736098	5	1293971	
```
The pipeline will take the effect sizes of the SNPs in this table, and as SNP ID the chromosome+position. If the flag --target_build is set to "hg38" there must be a column named 'chr_position_hg38', for "hg37" or "hg19" it would need a column "chr_position_hg37" or "chr_position_hg19".


Providing a phenotype file is optional, but without it there will be no ROC and OR deciles plots generated. 

An example `--phenotype` file looks like this:
```
patient_001	control
patient_002	case
```
i.e. it must be a headerless table with the sample ID on the first column and the phenotype on the second. Labels are case-insensitive, but must be 'case' or 'control'. IDs are matched to the pipeline sample_id column, independent of row order.


- Output layout
```
results/
├── prs_scores.tsv
├── prs_variant_details.tsv
├── per_sample/
│   ├── sample1.prs.tsv
│   ├── sample1.variants.tsv
│   ├── sample2.prs.tsv
│   └── sample2.variants.tsv
├── plots/
│   ├── sample1.contributions.png
│   ├── sample1.variant_status.png
│   ├── sample2.contributions.png
│   ├── sample2.variant_status.png
│   ├── prs_roc.png
│   └── prs_or_deciles.png
└── pipeline_info/
    ├── normalized_scorefile.tsv
    ├── scorefile_qc.tsv
    ├── execution_report_*.html
    ├── execution_timeline_*.html
    ├── execution_trace_*.txt
    └── pipeline_dag_*.html
```
The main report (`prs_scores.tsv`) contains one row per sample. `prs_variant_details.tsv` contains one row per sample and score variant, including genotype, dosage, contribution, harmonization status, and missing-site handling.
