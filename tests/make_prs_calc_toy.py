#!/usr/bin/env python3

from decimal import Decimal
from pathlib import Path
import shutil
import subprocess
import sys


# ---------------------------------------------------------------------------
# Output and compression tools
# ---------------------------------------------------------------------------

outdir = Path(sys.argv[1] if len(sys.argv) > 1 else "prs_calc_toy").expanduser()
outdir.mkdir(parents=True, exist_ok=True)
outdir = outdir.resolve()

bgzip = shutil.which("bgzip")
tabix = shutil.which("tabix")
bcftools = shutil.which("bcftools")

if not ((bgzip and tabix) or bcftools):
    raise SystemExit(
        "This generator requires either:\n"
        "  * bgzip and tabix, or\n"
        "  * bcftools\n\n"
        "For example:\n"
        "  mamba install -c bioconda htslib bcftools"
    )


# ---------------------------------------------------------------------------
# Score definition
#
# Tuple fields:
# chromosome, hg19 position, effect allele, other allele, weight,
# hg38 position, true hg19 reference allele
# ---------------------------------------------------------------------------

VARIANTS = [
    ("1",  77967507,  "A", "T", "0.127764345127619",  77501822,  "T"),
    ("3",  189357199, "G", "T", "0.110739353736705",  189639410, "G"),
    ("5",  1285974,   "A", "C", "0.222222086896245",  1285859,   "C"),
    ("5",  1320247,   "G", "A", "0.14128747381231",   1320132,   "G"),
    ("6",  26285867,  "C", "T", "0.0758448628853302", 26285639,  "C"),
    ("6",  29781049,  "T", "G", "0.0951996282391802", 29813272,  "T"),
    ("6",  31434111,  "G", "A", "0.223154591253269",  31466334,  "A"),
    ("8",  27344719,  "G", "A", "0.141137604413745",  27487202,  "G"),
    ("8",  32410110,  "G", "A", "0.124309867253109",  32552592,  "G"),
    ("9",  21830157,  "G", "A", "0.15454663000007",   21830158,  "A"),
    ("9",  22052068,  "G", "A", "0.165877847012963",  22052069,  "A"),
    ("10", 105687632, "C", "A", "0.150557804602624",  103927874, "A"),
    ("11", 118125625, "T", "C", "0.102138630956978",  118254910, "T"),
    ("12", 998819,    "G", "C", "0.145612965618366",  889653,    "G"),
    ("13", 32972626,  "T", "A", "0.472187430515042",  32398489,  "A"),
    ("15", 49376624,  "T", "G", "0.155057925298457",  49084427,  "T"),
    ("15", 78857986,  "G", "C", "0.260378679112338",  78565644,  "C"),
    ("15", 79185180,  "G", "A", "0.0918388891129015", 78892838,  "G"),
    ("19", 41353107,  "C", "T", "0.12279251169745",   40847202,  "T"),
]


# ---------------------------------------------------------------------------
# Synthetic sample design
#
# The order of each list is exactly the order in VARIANTS.
# Genotypes are effect-coded rather than literal VCF genotypes:
#
#   0 = other allele
#   1 = effect allele
#
# They are converted to true VCF REF/ALT indices when VCF records are
# written. Therefore, an effect-coded 1/1 can become VCF GT=0/0 when the
# effect allele is the hg19 reference allele.
#
# None means that the score locus is omitted entirely from that VCF.
# ---------------------------------------------------------------------------

SAMPLE_ORDER = [
    "sample1",
    "sample2",
    "sample3",
    "sample4",
    "sample5",
    "sample6",
    "sample7",
    "sample8",
    "sample9",
    "sample10",
]

PHENOTYPES = {
    "sample1": "case",
    "sample2": "case",
    "sample3": "control",
    "sample4": "case",
    "sample5": "control",
    "sample6": "control",
    "sample7": "case",
    "sample8": "control",
    "sample9": "case",
    "sample10": "control",
}

PROFILES = {
    # Maximum synthetic score: all 19 effect alleles homozygous.
    "sample1": ["1/1"] * 19,

    # High score, with two heterozygous and two entirely absent score sites.
    "sample2": [
        "1/1", "1/1", "1/1", "1/1", "0/1",
        "1/1", "1/1", "1/1", "0/1", "1/1",
        "1/1", None,  "1/1", "1/1", "1/1",
        "1/1", "1/1", None,  "1/1",
    ],

    # All score sites present and homozygous non-effect.
    "sample3": ["0/0"] * 19,

    # Intermediate/high synthetic case score; one site omitted.
    "sample4": [
        "0/1", "0/1", "1/1", "0/1", "0/0",
        "0/1", "1/1", "0/1", None,  "0/1",
        "0/1", "0/1", "0/0", "0/1", "1/1",
        "0/1", "1/1", "0/1", "0/1",
    ],

    # Low score: four heterozygous effect alleles and two omitted sites.
    "sample5": [
        "0/1", "0/0", "0/0", None,  "0/1",
        "0/0", "0/0", "0/0", "0/0", "0/0",
        "0/0", None,  "0/1", "0/0", "0/0",
        "0/0", "0/0", "0/1", "0/0",
    ],

    # Fifteen score sites omitted. One called effect allele is present.
    # Omitted effect-is-reference sites will contribute two effect copies
    # under missing-genotype=reference.
    "sample6": [
        None,  None,  "0/0", None,  "0/1",
        None,  None,  None,  None,  None,
        None,  None,  None,  None,  "0/0",
        None,  "0/0", None,  None,
    ],

    # Intermediate score: every site is present and heterozygous.
    "sample7": ["0/1"] * 19,

    # Low-score control with three heterozygous and two omitted sites.
    "sample8": [
        "0/0", "0/1", "0/0", "0/0", "0/0",
        "0/0", None,  "0/0", "0/1", "0/0",
        "0/0", "0/0", None,  "0/0", "0/0",
        "0/1", "0/0", "0/0", "0/0",
    ],

    # High-score case with five heterozygous and one omitted site.
    "sample9": [
        "1/1", "0/1", "1/1", "1/1", "0/1",
        "1/1", "1/1", None,  "1/1", "0/1",
        "1/1", "1/1", "0/1", "1/1", "1/1",
        "1/1", "0/1", "1/1", "1/1",
    ],

    # All score sites omitted. The expected result is the PRS of a
    # homozygous-hg19-reference sample.
    "sample10": [None] * 19,
}

# Additional synthetic variants not included in the score.
# Fields: chromosome, position, REF, ALT, GT, ID
EXTRAS = {
    "sample1": [
        ("2",  200000,  "G", "A", "0/1", "toy_extra_sample1_1"),
        ("20", 500000,  "C", "G", "1/1", "toy_extra_sample1_2"),
    ],
    "sample2": [
        ("1", 100000,   "G", "C", "0/1", "toy_extra_sample2_1"),
        ("7", 700000,   "A", "G", "1/1", "toy_extra_sample2_2"),
    ],
    "sample3": [
        ("4", 400000,   "C", "T", "0/1", "toy_extra_sample3_1"),
    ],
    "sample4": [
        ("14", 1400000, "T", "C", "0/1", "toy_extra_sample4_1"),
        ("21", 2100000, "G", "A", "1/1", "toy_extra_sample4_2"),
    ],
    "sample5": [
        ("5",  2000000, "T", "C", "0/1", "toy_extra_sample5_1"),
        ("18", 1800000, "A", "G", "0/1", "toy_extra_sample5_2"),
    ],
    "sample6": [
        ("2",  600000,  "A", "C", "1/1", "toy_extra_sample6_1"),
        ("22", 2200000, "C", "T", "0/1", "toy_extra_sample6_2"),
    ],
    "sample7": [
        ("3",  300000,  "A", "G", "0/1", "toy_extra_sample7_1"),
        ("16", 1600000, "C", "T", "1/1", "toy_extra_sample7_2"),
    ],
    "sample8": [
        ("8",  800000,  "G", "A", "0/1", "toy_extra_sample8_1"),
        ("17", 1700000, "T", "C", "0/1", "toy_extra_sample8_2"),
    ],
    "sample9": [
        ("9",  900000,  "C", "G", "1/1", "toy_extra_sample9_1"),
        ("20", 2000000, "A", "T", "0/1", "toy_extra_sample9_2"),
    ],
    "sample10": [
        ("10", 1000000, "G", "T", "0/1", "toy_extra_sample10_1"),
        ("22", 2500000, "C", "A", "1/1", "toy_extra_sample10_2"),
    ],
}

# Dosage in the effect-coded PROFILES representation.
DOSAGE = {
    "0/0": 0,
    "0/1": 1,
    "1/0": 1,
    "1/1": 2,
}

DEFAULT_PLOIDY = 2


# If the effect allele is VCF REF, profile allele indices must be swapped:
#
# profile index 0 = other allele -> VCF index 1
# profile index 1 = effect allele -> VCF index 0
#
EFFECT_REFERENCE_GT = {
    "0/0": "1/1",
    "0/1": "1/0",
    "1/0": "0/1",
    "1/1": "0/0",
}


def profile_gt_to_vcf_gt(gt, effect_is_reference):
    if effect_is_reference:
        return EFFECT_REFERENCE_GT[gt]
    return gt


# Ensure every score model is biallelic and that its declared reference
# allele is one of the two score alleles.
for variant in VARIANTS:
    (
        chrom,
        pos19,
        effect,
        other,
        _weight,
        _pos38,
        reference,
    ) = variant

    if effect == other:
        raise ValueError(
            f"{chrom}:{pos19} has identical effect and other alleles"
        )

    if reference not in (effect, other):
        raise ValueError(
            f"{chrom}:{pos19} hg19 reference allele {reference} "
            f"is not one of the score alleles {effect}/{other}"
        )


for sample in SAMPLE_ORDER:
    if len(PROFILES[sample]) != len(VARIANTS):
        raise ValueError(
            f"{sample} does not have exactly "
            f"{len(VARIANTS)} profile entries"
        )
    invalid = [
        gt for gt in PROFILES[sample]
        if gt is not None and gt not in DOSAGE
    ]
    if invalid:
        raise ValueError(f"Invalid genotypes for {sample}: {invalid}")


# ---------------------------------------------------------------------------
# Write score file
# ---------------------------------------------------------------------------

score_metadata = """###PGS CATALOG SCORING FILE - see https://www.pgscatalog.org/downloads/#dl_ftp_scoring for additional information
#format_version=2.0
##POLYGENIC SCORE (PGS) INFORMATION
#pgs_id=PGS000392
#pgs_name=PRSWEB_PHECODE165.1_GWAS-Catalog-r2019-05-03-X165.1_PT_UKB_20200608
#trait_reported=Lung and bronchus cancer
#trait_mapped=tracheal cancer|bronchus cancer|lung cancer
#trait_efo=MONDO_0001407|MONDO_0001672|MONDO_0008903
#genome_build=GRCh37
#variants_number=19
#weight_type=NR
##SOURCE INFORMATION
#pgp_id=PGP000118
#citation=Fritsche LG et al. Am J Hum Genet (2020). doi:10.1016/j.ajhg.2020.08.025
##HARMONIZATION DETAILS
#HmPOS_build=GRCh38
#HmPOS_date=2022-07-29
#HmPOS_match_chr={"True": null, "False": null}
#HmPOS_match_pos={"True": null, "False": null}
"""

scorefile = outdir / "PGS000392_hg19.txt"

with scorefile.open("w") as handle:
    handle.write(score_metadata)
    handle.write(
        "\t".join([
            "chr_name",
            "chr_position_hg19",
            "effect_allele",
            "other_allele",
            "effect_weight",
            "hm_source",
            "hm_rsID",
            "hm_chr",
            "chr_position_hg38",
            "hm_inferOtherAllele",
        ]) + "\n"
    )

    for (
        chrom,
        pos19,
        effect,
        other,
        weight,
        pos38,
        _reference,
    ) in VARIANTS:
        fields = [
            chrom,
            pos19,
            effect,
            other,
            weight,
            "liftover",
            "",
            chrom,
            pos38,
            "",
        ]
        handle.write("\t".join(str(value) for value in fields) + "\n")


# ---------------------------------------------------------------------------
# Write phenotype table
# ---------------------------------------------------------------------------

phenotype_file = outdir / "phenotypes.tsv"

with phenotype_file.open("w") as handle:
    for sample in SAMPLE_ORDER:
        handle.write(f"{sample}\t{PHENOTYPES[sample]}\n")


# ---------------------------------------------------------------------------
# Write genotype audit and calculate expected raw additive scores
# ---------------------------------------------------------------------------

genotype_design = outdir / "genotype_design.tsv"
summary = {}

with genotype_design.open("w") as handle:
    handle.write(
        "\t".join([
            "sample_id",
            "phenotype",
            "chrom",
            "position_hg19",
            "other_allele",
            "effect_allele",
            "vcf_ref_if_present",
            "vcf_alt_if_present",
            "GT",
            "site_status",
            "expected_effect_dosage",
            "expected_contribution",
        ]) + "\n"
    )

    for sample in SAMPLE_ORDER:
        total = Decimal("0")
        total_dosage = 0
        omitted = 0

        for variant, gt in zip(VARIANTS, PROFILES[sample]):
            (
                chrom,
                pos19,
                effect,
                other,
                weight,
                _pos38,
                reference,
            ) = variant

            effect_is_reference = effect == reference
            vcf_alt = (
                other
                if effect_is_reference
                else effect
            )

            if gt is None:
                # An absent variant-only VCF row is interpreted as
                # homozygous genomic reference.
                dosage = (
                    DEFAULT_PLOIDY
                    if effect_is_reference
                    else 0
                )
                gt_text = "ABSENT"
                status = "absent_assumed_hom_ref"
                omitted += 1
            else:
                dosage = DOSAGE[gt]
                gt_text = profile_gt_to_vcf_gt(
                    gt,
                    effect_is_reference,
                )
                status = "present"

            contribution = Decimal(weight) * dosage
            total += contribution
            total_dosage += dosage

            handle.write(
                "\t".join([
                    sample,
                    PHENOTYPES[sample],
                    chrom,
                    str(pos19),
                    other,
                    effect,
                    reference,
                    vcf_alt,
                    gt_text,
                    status,
                    str(dosage),
                    format(contribution, "f"),
                ]) + "\n"
            )

        summary[sample] = {
            "score": total,
            "dosage": total_dosage,
            "omitted": omitted,
        }


expected_scores = outdir / "expected_scores.tsv"

with expected_scores.open("w") as handle:
    handle.write(
        "sample_id\tphenotype\texpected_raw_prs\t"
        "effect_allele_count\tomitted_score_sites\n"
    )
    for sample in SAMPLE_ORDER:
        values = summary[sample]
        handle.write(
            f"{sample}\t{PHENOTYPES[sample]}\t"
            f"{format(values['score'], 'f')}\t"
            f"{values['dosage']}\t{values['omitted']}\n"
        )


extra_variants = outdir / "extra_variants.tsv"

with extra_variants.open("w") as handle:
    handle.write("sample_id\tchrom\tposition\tREF\tALT\tGT\tID\n")
    for sample in SAMPLE_ORDER:
        for chrom, pos, ref, alt, gt, identifier in EXTRAS[sample]:
            handle.write(
                f"{sample}\t{chrom}\t{pos}\t{ref}\t{alt}\t"
                f"{gt}\t{identifier}\n"
            )


# ---------------------------------------------------------------------------
# VCF generation helpers
# ---------------------------------------------------------------------------

def contig_key(chrom):
    if chrom.isdigit():
        return (0, int(chrom))
    return (1, chrom)


def remove_if_present(path):
    if path.exists():
        path.unlink()


def compress_and_index(vcf_path):
    gz_path = Path(str(vcf_path) + ".gz")

    remove_if_present(gz_path)
    remove_if_present(Path(str(gz_path) + ".tbi"))
    remove_if_present(Path(str(gz_path) + ".csi"))

    if bgzip and tabix:
        subprocess.run(
            [bgzip, "-f", str(vcf_path)],
            check=True,
        )
        subprocess.run(
            [tabix, "-f", "-p", "vcf", str(gz_path)],
            check=True,
        )
    else:
        subprocess.run(
            [bcftools, "view", "-Oz", "-o", str(gz_path), str(vcf_path)],
            check=True,
        )
        subprocess.run(
            [bcftools, "index", "-f", "-t", str(gz_path)],
            check=True,
        )
        vcf_path.unlink()

    return gz_path


# ---------------------------------------------------------------------------
# Write one minimal single-sample VCF per sample
# ---------------------------------------------------------------------------

vcf_dir = outdir / "vcfs"
vcf_dir.mkdir(parents=True, exist_ok=True)

score_contigs = {variant[0] for variant in VARIANTS}
compressed_vcfs = {}

for sample in SAMPLE_ORDER:
    records = []

    for index, variant in enumerate(VARIANTS):
        (
            chrom,
            pos19,
            effect,
            other,
            _weight,
            _pos38,
            reference,
        ) = variant
        gt = PROFILES[sample][index]

        if gt is None:
            # No VCF row is emitted. This tests absent-site handling.
            continue

        effect_is_reference = effect == reference
        vcf_alt = (
            other
            if effect_is_reference
            else effect
        )
        vcf_gt = profile_gt_to_vcf_gt(
            gt,
            effect_is_reference,
        )

        records.append((
            chrom,                   # CHROM
            pos19,                   # POS
            f"{chrom}:{pos19}",      # ID
            reference,               # true hg19 REF
            vcf_alt,                 # non-reference allele
            ".",                     # QUAL
            "PASS",                  # FILTER
            ".",                     # INFO
            "GT",                    # FORMAT
            vcf_gt,                  # true VCF-indexed sample value
        ))

    for chrom, pos, ref, alt, gt, identifier in EXTRAS[sample]:
        records.append((
            chrom,
            pos,
            identifier,
            ref,
            alt,
            ".",
            "PASS",
            ".",
            "GT",
            gt,
        ))

    records.sort(key=lambda row: (contig_key(row[0]), row[1]))

    # Include all score chromosomes in every header, including chromosomes
    # on which that particular sample has no score records.
    header_contigs = score_contigs | {
        extra[0] for extra in EXTRAS[sample]
    }

    vcf_path = vcf_dir / f"{sample}.vcf"

    with vcf_path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=prs_calc_synthetic_test_fixture\n")

        for chrom in sorted(header_contigs, key=contig_key):
            handle.write(f"##contig=<ID={chrom}>\n")

        handle.write(
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
        )
        handle.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,'
            'Description="Genotype">\n'
        )
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            f"{sample}\n"
        )

        for record in records:
            handle.write("\t".join(str(value) for value in record) + "\n")

    compressed_vcfs[sample] = compress_and_index(vcf_path)


# ---------------------------------------------------------------------------
# Explicit two-column, headerless sample sheet
# ---------------------------------------------------------------------------

samplesheet = outdir / "samplesheet.tsv"

with samplesheet.open("w") as handle:
    for sample in SAMPLE_ORDER:
        handle.write(f"{sample}\t{compressed_vcfs[sample]}\n")


print(f"Created fixture under: {outdir}")
print()
print("Expected raw additive scores:")
for sample in SAMPLE_ORDER:
    values = summary[sample]
    print(
        f"  {sample:7s} {PHENOTYPES[sample]:7s} "
        f"score={format(values['score'], 'f')} "
        f"dosage={values['dosage']} "
        f"omitted={values['omitted']}"
    )
