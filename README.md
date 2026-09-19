Miguel Martín Álvarez, PhD

miguel.m.alvarez3[--at--]gmail[--dot--]com

===========================================

- Example script to run the pipeline:

```
#!/usr/bin/env bash

set -euo pipefail

unset R_HOME
export NXF_HOME="$HOME/work/.nextflow"

nextflow run mmaalvarez/prs_calc -r main -latest \
	--input /path/to/samplesheet.tsv \
	--scorefile /path/to/PGSXXXXX_hg38.txt \
	--phenotypes /path/to/phenotypes.tsv \
	--target_build hg38 \
	-profile conda \
	-resume
```

For faster Conda environment creation, if mamba is installed:

`-profile mamba`

For Slurm plus Conda:

`-profile conda,slurm`


- Input formats

A one-column file fed to --input works directly:
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


An example input for --scorefile looks like this:
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
The pipeline will take the effect sizes of the SNPs in this table, and as SNP ID the chromosome+position (so if the flag --hg says "38" then there needs to be a column named 'chr_position_hg38', if it said "37" or "19" it would need a column "chr_position_hg37" or "chr_position_hg19".


An example --phenotype file (optional, but if not provided it will not generate ROC and OR deciles plots) looks like this:
```
patient_001	control
patient_002	case
```
i.e. a table with the sample ID on the first column and the phenotype on the second. The phenotype file may be headerless. Labels are case-insensitive, but must be 'case' or 'control'. IDs are matched to the pipeline sample_id column, independent of row order.


- Output layout

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

The main report, prs_scores.tsv, contains one row per sample. prs_variant_details.tsv contains one row per sample and score variant, including genotype, dosage, contribution, harmonization status, and missing-site handling.
