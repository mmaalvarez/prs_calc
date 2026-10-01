# prs_calc

A Nextflow pipeline that calculates a per-sample polygenic risk score (PRS)
from a SNP score file and one-sample VCF, gVCF, or BCF files. It produces
variant-level scoring details, per-sample plots, and, optionally, cohort ROC
and PRS-decile odds-ratio plots.

For each retained score row, the pipeline calculates:

```text
contribution = effect_weight × effect-allele dosage
PRS          = sum of contributions
```

The dosage is the number of called genotype alleles belonging to the score
row's `effect_allele` set. The final result is a raw PRS; to standardize it, 
calculate the mean and standard deviation of raw PRSs in an appropriate 
reference population, then use:

Z_{\mathrm{PRS}}=\frac{\mathrm{PRS}_{\mathrm{raw}}-\mu_{\mathrm{reference}}} {\sigma_{\mathrm{reference}}}


## Requirements and quick start

Nextflow 23.10.0 or newer is required. The `conda` and `mamba` profiles use
`envs/prs_calc.yml` to provide bcftools, Python, R packages, and the hg19 and
hg38 BSgenome references.

Example command:
```bash
nextflow run mmaalvarez/prs_calc -r main -latest \
  --input /path/to/samplesheet.tsv \
  --gvcf_mode auto \
  --scorefile /path/to/scorefile.tsv \
  --target_build hg38 \
  --missing_genotype error \
  --no_calls error \
  -profile conda \
  -resume
```

If you downloaded the latest release to run the pipeline locally, the run 
command should instead be:
```bash
nextflow run path/to/prs_calc/main.nf \
  (...)
  (remaining parameters don't change)
```

**Both `--missing_genotype` and `--no_calls` must be supplied.** The example's
`error` policies require every score position to have usable evidence; a
variant-only VCF that omits score positions will fail. See [Missing genotypes
and no-calls](#missing-genotypes-and-no-calls) before choosing different
policies.

To request cohort plots, also provide `--phenotypes /path/to/phenotypes.tsv`.
For Slurm, use `-profile conda,slurm`; `-profile mamba` is available when
Mamba is installed.

## Input VCFs and samplesheet

Each input VCF, gVCF, or BCF must be readable by bcftools and contain
**exactly one sample**. The pipeline also requires `##contig` declarations
in its header so it can map score chromosomes to VCF contig names.

The `--input` file is a tab-separated samplesheet with either one path per
line:

```text
/path/to/sample1.vcf.gz
/path/to/sample2.gvcf.gz
```

or an explicit pipeline sample ID followed by the path:

```text
patient_001	/path/to/sample1.vcf.gz
patient_002	/path/to/sample2.gvcf.gz
```

Relative VCF paths are resolved relative to the samplesheet directory. Blank
lines, lines beginning with `#`, and simple column headers are ignored.
Without an explicit ID, the ID is inferred from the filename by removing
its `.vcf`, `.gvcf`, or `.bcf` suffix and optional `.gz`/`.bgz` suffix.

Pipeline IDs are converted to filename-safe IDs; collisions after this
conversion are errors. The ID inside the VCF can differ from the pipeline
ID. Phenotypes must use the **pipeline** IDs, not necessarily the VCF-header
sample names.

Chromosome spellings such as `1` and `chr1` are mapped to the same score
chromosome. The pipeline queries all matching contig spellings declared in
the VCF header; conflicting records at the same normalized score position
cause an error. If none of the score chromosomes map to VCF contigs, the
sample fails. Successful per-sample reports appear under 
`pipeline_info/target_mapping/`. **Only primary chromosomes** are accepted 
(e.g. chr1_KI270706v1_random is not accepted)

## Score file

`--scorefile` is a tab-separated score table. Lines beginning with `#`
before the header are permitted. Required columns are:

- `effect_allele`: an A/C/G/T base or slash-separated set, such as `A/T`.
- `effect_weight`: a finite numeric regression coefficient, **not** an
  untransformed odds ratio.
- `chr_position_<target_build>` (see [Genome build and genotype policies](#genome-build-and-genotype-policies) below).
- Either `chr_name` or `hm_chr`. When both exist, a nonmissing `chr_name`
  takes precedence.

`other_allele` is optional. If it is missing and `hm_inferOtherAllele` is
present, a valid inferred value is used; invalid, nonempty `hm_inferOtherAllele` 
values are warned about, left as unknown `other_allele`, and counted in score-file 
QC; syntactically valid values overlapping the effect set trigger an error.
Supplying `other_allele` permits stricter checking that concrete VCF REF and 
ALT alleles belong to the declared score allele sets. IDs are taken, in order, 
from `rsID`, `hm_rsID`, or `variant_id` when available; otherwise an ID is 
constructed. **Variants are matched by normalized chromosome and position, not 
by rsID.**

For example:

```text
rsID	chr_name	chr_position_hg38	effect_allele	other_allele	effect_weight
rs71658797	1	77501822	A	T	0.131028262406404
rs13080835	3	189639410	T	G	-0.0618754037180875
```

The normalizer uppercases and sorts allele sets, skips score rows with
unsupported or non-SNP allele strings, and records those skips in
`scorefile_qc.tsv`. Missing effect alleles, invalid positions, and 
overlapping effect/other allele sets cause errors. If no score rows
remain, the run fails.

**Duplicate score coordinates are not deduplicated.** Each retained score
row contributes separately to the PRS. The score-file QC report counts
positions with multiple score rows.

### Slash-separated allele sets

A **slash in the score file** denotes a *set of acceptable bases*, **not a VCF
genotype**. For a diploid genotype, each called copy belonging to the effect
set adds one to the dosage; for example:

| Effect set | Other set | VCF REF/ALT | GT | Dosage |
| --- | --- | --- | --- | ---: |
| `C` | `A/T` | `C/A` | `0/0` | 2 |
| `C` | `A/T` | `C/A` | `0/1` | 1 |
| `C` | `A/T` | `A/T` | `0/1` | 0 |
| `A/T` | `C` | `A/T` | `0/1` | 2 |
| `A/T` | `C` | `A/C` | `0/1` | 1 |
| `A/T` | `C` | `C/A,T` | `1/2` | 2 |

When `other_allele` is specified, an unexpected concrete VCF REF or ALT at
the scored SNP is an allele mismatch and aborts the sample. The pipeline
does not automatically complement alleles or resolve strand flips.

## Genome build and genotype policies

`--target_build` accepts either **hg38** or **hg19**. It defaults to `hg38`
(`GRCh38` and `38` are valid as aliases thereof); likewise, `GRCh37`, `37`, 
and `19` are all aliases of `hg19`. `--hg` is an optional alias that takes 
precedence when supplied. The selected build determines **both** the required 
score position column and the BSgenome reference used during scoring.

Choose a score file and VCF that genuinely use compatible coordinates
and reference assemblies. Merely setting `--target_build` does not
liftover either input.

### VCF versus gVCF detection

`--gvcf_mode` defaults to `auto`, and tries to infer whether the input has a gVCF 
or a VCF format. However, a VCF containing only ordinary SNP records cannot be proved to be 
plain rather than a gVCF subset with its blocks removed; `--gvcf_mode auto` will fail 
for that ambiguous input, so you'll need to specify `--gvcf_mode plain` or `--gvcf_mode gvcf` 
explicitly.

### Missing genotypes and no-calls

| Situation | Policy |
| --- | --- |
| Called SNP genotype | Count the copies of bases in the effect set. |
| Validated reference-only hom-ref record or block | Score using its hom-ref GT ploidy, regardless of the missing-site policy. |
| No matching record or hom-ref block | Apply `--missing_genotype`. |
| Explicit wholly missing SNP GT, such as `./.` | Apply `--no_calls`. |

`--missing_genotype` is always required; it takes the following options:

- `error`: abort at the first absent/uncovered score position.
- `zero`: assign effect-allele dosage 0. This is an imputation policy,
  **not evidence that the sample has zero effect alleles**.
- `reference`: **plain VCF only**. Assume the absent position is
  homozygous for the selected BSgenome reference base, using
  `--default_ploidy` (default: 2): It must be an integer. It is rejected 
  for gVCF input.

Validated reference-only hom-ref evidence is honoured even when
`--missing_genotype error` or `zero` was selected. A reference-only gVCF
record with `END` can cover score positions after its starting position;
the pipeline validates its REF at the actual block start and takes ploidy
from its GT. The genome reference base is then looked up separately at each 
scored position within the block.

`--no_calls` is always required as well; it takes the following options:

- `error`: abort on an explicit wholly missing score-SNP GT.
- `zero`: assign that score SNP dosage 0.

A no-call is **not** treated as a hom-ref block or as an absent position.
Reference-only no-call records with `END` can mark an interval of
no-calls.

For partially missing GTs, such as `0/.`, `--no_calls zero` counts effect 
alleles among the known GT entries, and the missing entry contributes zero. 
This is reported in the published per-sample QC (prs_scores.tsv). For e.g. 
REF=G, ALT=A, effect=G, weight=0.25, this gives:

| GT | called_alleles in details | Dosage | Contribution |
| --- | --- | --- | --- |
| 0/. | G/. | 1 | 0.25 | 
| 1/. | A/. | 0 | 0 |
| ./. | missing | 0 | 0 |
| 0/0 | G/G | 2 | 0.50 |

If `--no_calls error`, a partially missing GT aborts.

Other malformed GTs always abort. Allele-incompatible sites can still abort 
even if a `zero` policy was selected. Conversely, zero-imputed positions 
do not become observed genotypes.

## Coverage and quality control

`--min_covered_fraction` accepts a number from 0 to 1 and defaults to
`0`, which disables its failure threshold. Coverage is calculated across
**retained score rows**, not unique positions:

```text
covered_fraction =
  (rows with usable called SNP GTs + rows covered by validated hom-ref evidence)
  / retained score rows
```

No-calls assigned zero and absent positions assigned zero or assumed
reference are **not** covered. Partial GTs are also excluded from 
covered_fraction, even though they contribute to the raw PRS, because 
usable_observed_call explicitly lists only "observed" and 
"observed_nonmodel_allele". A threshold such as `--min_covered_fraction 0.95` 
fails a sample below 95% usable coverage. This check occurs after that 
task writes its TSVs: they can be inspected in the failed task's Nextflow 
work directory, but a failed run will not produce the merged PRS report.

In a successful summary, `n_scored_variants` equals the number of retained
score rows, **including zero-imputed rows**. Do not interpret it as the
number of observed genotypes. Inspect `covered_fraction`,
`observed_fraction`, `n_unusable_positions`, `n_uncovered_positions`,
`n_vcf_no_calls`, and `n_dosage_from_reference` instead.
`n_dosage_from_reference` counts *assumed* reference genotypes, not
positions scored from validated hom-ref blocks.

`n_absent_from_vcf` is not the uncovered count: it also includes
block-covered score positions for which there was no separately selected
SNP record. Use `n_uncovered_positions` when assessing positions with
neither a selected record nor a hom-ref block.

Genotype quality and filtering columns are ignored, so **make sure the input 
gVCF/VCF contains only good quality SNPs**.

## Optional phenotypes and cohort plots

Supply `--phenotypes` to request cohort plots. The file must contain
exactly two tab-separated columns: pipeline sample ID and phenotype.
It can be headerless:

```text
patient_001	control
patient_002	case
```

A simple `sample_id`/`phenotype` header is also accepted. Labels are
case-insensitive and must be `case` or `control`. Every scored sample
needs a phenotype; duplicate phenotype IDs cause an error. Extra
phenotype IDs are warned about and ignored. Both a case and a control
are required for cohort plotting.

The ROC curve treats **higher PRS as more case-like**; its direction is
not automatically reversed, so an AUC below 0.5 is possible. For the
decile plot, samples are ordered by PRS, with sample ID breaking ties
deterministically. Each decile's unadjusted case/control odds is compared
with the lowest-PRS decile; error bars are 95% Wald intervals. A
0.5 continuity correction is used for a comparison containing a zero
cell. Fewer than 10 samples trigger a warning and cannot populate all
ten deciles. These descriptive cohort plots do not adjust for
covariates or validate clinical predictive performance.

Without `--phenotypes`, per-sample plots are still produced, but the
cohort ROC and odds-ratio plots are not run.

## Outputs

With the default `--outdir`, successfully published outputs are:

```text
results/
├── prs_scores.tsv
├── prs_variant_details.tsv
├── plots/
│   ├── per_sample/
│   │   ├── SAMPLE.contributions.png
│   │   └── SAMPLE.variant_status.png
│   ├── prs_roc.png                 # only with --phenotypes
│   └── prs_or_deciles.png          # only with --phenotypes
└── pipeline_info/
    ├── target_mapping/
    ├── normalized_scorefile.tsv
    ├── scorefile_qc.tsv
    └── execution_*                 # Nextflow reports, trace, timeline, DAG
```

`prs_scores.tsv` contains one row per sample; `prs_variant_details.tsv`
contains one row per sample and retained score row. Details include the
allele dosage, contribution, genotype when available, scoring status,
and whether a position was covered by hom-ref evidence or required an
assumption. Per-sample contribution plots show the largest absolute
contributions; `--plot_top_variants` sets their maximum count
(default: 25). Status plots count variant-scoring statuses.

The default publishing mode is `copy`; it can be changed with
`--publish_dir_mode`.

## Important limitations

- This is a **SNP scoring** pipeline. Unsupported score alleles are
  skipped with QC accounting; a non-SNP VCF record whose REF span or declared 
  END overlaps a score position, including upstream starts, aborts. 
  **Complex SV representations without a usable linear span are not 
  comprehensively handled**.
- Called GT ploidy and validated hom-ref-block GT ploidy are read from
  their genotypes, but assumed reference genotypes use the specified 
  `--default_ploidy` everywhere: there is no chromosome-, sex-, or
  pseudoautosomal-region-specific ploidy model, so `--default_ploidy 2` 
  can mis-score sex chromosomes or other haploid contexts. Also, the assumed 
  interchangeability of GRCh37 and UCSC hg19 does not apply to mitochondrial 
  DNA. **Treat X, Y, and mitochondrial scores with particular caution**.
- A concrete SNP whose REF allele in the VCF disagrees with the REF allele
  in the BSgenome because of **a simple REF/ALT swap is warned about, 
  not rejected**; e.g.
  | BSgenome | VCF | Action |
  | --- | --- | --- |
  | C | REF=A, ALT=C | swap-compatible, continue and warn |
  | G | REF=A, ALT=C | unresolved, abort|
  Reference-only (i.e. with a symbolic ALT allele like `<NON_REF>`) hom-ref 
  blocks and no-call interval starts abort in such case. Resolve genome-build 
  and orientation warnings before interpreting a PRS, especially with 
  `--missing_genotype reference`.
- An absent score row assigned `zero` can have neither score allele matching 
  BSgenome and it will still pass.
- Strand harmonization is not performed. A palindromic A/T or C/G SNP cannot 
  be reliably oriented from alleles alone: **a strand flip interpreted as 
  swap-compatible will invert its dosage**.
- The default coverage threshold is disabled. Selecting a zero or
  reference policy can therefore yield a **numeric PRS with substantial
  unobserved input**. Review QC and set an appropriate coverage
  requirement for the intended use.
- Several "counts" **count score rows, not unique SNP positions or VCF records**. 
  This includes n_vcf_records_found, partial-GT counts, and coverage; duplicated 
  score coordinates contribute repeatedly.
- A multiallelic SNP split into several biallelic records —for example by 
  `bcftools norm -m -`— will generally fail at a scored position. Multiallelic 
  sites should remain represented jointly, or implement deliberate reconciliation 
  of split records.

Maintainer: Miguel Martín Álvarez, PhD  
Contact: miguel.m.alvarez3[--at--]gmail[--dot--]com
