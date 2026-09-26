#!/usr/bin/env python3
import argparse
import csv
import subprocess
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="prs_calc_toy")
ap.add_argument(
    "--cohort",
    action="store_true",
    help="Also generate ten tiny VCFs and phenotypes to smoke-test cohort plots",
)
args = ap.parse_args()

root = Path(args.out).resolve()
root.mkdir(parents=True, exist_ok=True)

START = 20_000_000
END = START + 255
CHR22_LENGTH = 50_818_468
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
NONPAL_ALT = {"A": "C", "C": "A", "G": "A", "T": "C"}

# Use precisely the hg38 package used by calculate_prs.R.
r_expr = (
    "suppressPackageStartupMessages({"
    "library(Biostrings);library(BSgenome.Hsapiens.UCSC.hg38)"
    "});"
    "g <- get('Hsapiens',"
    "as.environment('package:BSgenome.Hsapiens.UCSC.hg38'));"
    "cat(as.character(Biostrings::getSeq("
    "g,names='chr22',start=20000000L,end=20000255L)))"
)
dna = subprocess.check_output(
    ["Rscript", "-e", r_expr], text=True
).strip().upper()

if len(dna) != END - START + 1:
    raise RuntimeError("Did not retrieve the expected 256 hg38 bases")

(root / "hg38_chr22_chunk.fa").write_text(
    ">chr22:20000000-20000255 hg38\n"
    + "\n".join(dna[i:i + 60] for i in range(0, len(dna), 60))
    + "\n"
)

def base_at(position):
    return dna[position - START]

# Space sites out. Ensure one real A/T-reference position and one real
# C/G-reference position for the two palindromic allele pairs.
positions = []
cursor = START + 8
for i in range(1, 13):
    allowed = (
        "AT" if i == 4
        else "CG" if i == 5
        else "ACGT"
    )
    found = next(
        (
            p for p in range(cursor, END - 7)
            if base_at(p) in allowed
        ),
        None,
    )
    if found is None:
        raise RuntimeError(
            "Could not place all sites in this hg38 interval; "
            "choose a different 256-bp interval"
        )
    positions.append(found)
    cursor = found + 12

sites = {}
for i, p in enumerate(positions, start=1):
    ref = base_at(p)
    alt = COMP[ref] if i in (4, 5) else NONPAL_ALT[ref]
    third = next(x for x in "ACGT" if x not in (ref, alt))
    sites[i] = {
        "pos": p, "ref": ref, "alt": alt, "third": third
    }

extra_pos = next(
    p for p in range(START, END + 1)
    if p not in positions and base_at(p) in "ACGT"
)

weights = [
    0.10, -0.20, 0.30, 0.40, -0.50, 0.60,
    0.70, -0.80, 0.90, 0.15, 0.05, 0.25,
]
ref_effect_rows = {2, 5, 6, 8, 11, 12}

score_fields = [
    "rsID", "chr_name", "chr_position_hg38",
    "effect_allele", "other_allele", "effect_weight",
    "hm_inferOtherAllele",
]
score_rows = []
for i in range(1, 13):
    s = sites[i]
    effect = s["ref"] if i in ref_effect_rows else s["alt"]
    other = (
        "" if i == 10
        else s["alt"] if i in ref_effect_rows
        else s["ref"]
    )
    score_rows.append({
        "rsID": f"s{i:02}",
        "chr_name": "22",
        "chr_position_hg38": str(s["pos"]),
        "effect_allele": effect,
        "other_allele": other,
        "effect_weight": f"{weights[i - 1]:.2f}",
        # Deliberately ambiguous inferred other alleles in row 10.
        "hm_inferOtherAllele": "A/C/G" if i == 10 else "",
    })

def write_score(name, rows, fields=score_fields):
    with (root / name).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t",
            lineterminator="\n", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

write_score("score.tsv", score_rows)

pal_swapped = [dict(x) for x in score_rows]
pal_swapped[3]["effect_allele"], pal_swapped[3]["other_allele"] = (
    pal_swapped[3]["other_allele"],
    pal_swapped[3]["effect_allele"],
)
write_score("score_pal_effect_swapped.tsv", pal_swapped)

bad_position = [dict(x) for x in score_rows]
bad_position[0]["chr_position_hg38"] = "not_an_integer"
write_score("score_bad_position.tsv", bad_position)

non_snv = [dict(x) for x in score_rows]
non_snv[0]["effect_allele"] = "AT"
write_score("score_non_snv.tsv", non_snv)

duplicate = [dict(x) for x in score_rows]
duplicate.append({**score_rows[0], "rsID": "s01_duplicate"})
write_score("score_duplicate_position.tsv", duplicate)

write_score("score_no_call_only.tsv", [score_rows[7]])
write_score(
    "score_missing_weight.tsv",
    score_rows,
    fields=[x for x in score_fields if x != "effect_weight"],
)

with (root / "sites.tsv").open("w", newline="") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow([
        "score_id", "hg38_position", "hg38_ref",
        "effect_allele", "other_allele", "first_alt",
    ])
    for i in range(1, 13):
        w.writerow([
            f"s{i:02}", sites[i]["pos"], sites[i]["ref"],
            score_rows[i - 1]["effect_allele"],
            score_rows[i - 1]["other_allele"],
            sites[i]["alt"],
        ])

# VCF record: POS, REF, ALT, INFO, GT.
def record(i, gt, *, nonref=False):
    s = sites[i]
    alts = [s["alt"]]
    if i == 9:
        alts.append(s["third"])  # Real multiallelic GT=1/2.
    if nonref:
        alts.append("<NON_REF>")
    return (s["pos"], s["ref"], ",".join(alts), ".", gt)

def extra_record(*, nonref=False):
    ref = base_at(extra_pos)
    alt = NONPAL_ALT[ref]
    if nonref:
        alt += ",<NON_REF>"
    return (extra_pos, ref, alt, ".", "0/1")

BASE_CALLS = {
    1: "0/1", 2: "1/1", 3: "1/1", 4: "1/1",
    5: "0/1", 8: "0/0", 9: "1/2", 10: "0/1",
    11: "0/0",
}

def plain_records(calls=None):
    if calls is None:
        calls = BASE_CALLS
    return [extra_record()] + [
        record(i, gt) for i, gt in calls.items()
    ]

def write_vcf(
    relpath, records, *, gvcf=False, star=False,
    contigs=("chr22",), samples=None, deletion_alt=False,
):
    dest = root / relpath
    dest.parent.mkdir(parents=True, exist_ok=True)
    raw = dest if dest.suffix == ".vcf" else dest.with_suffix("")
    sample_names = tuple(samples) if samples is not None else (
        dest.name.split(".")[0],
    )

    with raw.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig in contigs:
            out.write(
                f"##contig=<ID={contig},length={CHR22_LENGTH}>\n"
            )
        if gvcf:
            out.write(
                '##INFO=<ID=END,Number=1,Type=Integer,'
                'Description="End of reference block">\n'
            )
            if not star:
                out.write(
                    '##ALT=<ID=NON_REF,'
                    'Description="Unspecified alternative allele">\n'
                )
        if deletion_alt:
            out.write(
                '##ALT=<ID=DEL,Description="Symbolic deletion">\n'
            )
        out.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,'
            'Description="Genotype">\n'
        )
        out.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
            "INFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
        )
        for p, ref, alt, info, gt in sorted(
            records, key=lambda row: row[0]
        ):
            fields = [
                "chr22", str(p), ".", ref, alt, ".", "PASS",
                info, "GT",
            ] + [gt] * len(sample_names)
            out.write("\t".join(fields) + "\n")

    if dest != raw:
        subprocess.run(
            ["bcftools", "view", "-Oz", "-o", str(dest), str(raw)],
            check=True,
        )
        subprocess.run(
            ["bcftools", "index", "-t", "-f", str(dest)],
            check=True,
        )
        raw.unlink()

clean = plain_records()

no_call = dict(BASE_CALLS)
no_call[8] = "./."

swappable = {sites[i]["pos"] for i in (3, 4, 5)}
swap = []
for row in clean:
    p, ref, alt, info, gt = row
    if p in swappable:
        # Deliberately INVALID against hg38. Recode GT 0<->1 so the
        # underlying called allele dosage is unchanged.
        recoded = gt.translate(str.maketrans("01", "10"))
        swap.append((p, alt, ref, info, recoded))
    else:
        swap.append(row)

mismatch = [
    (
        p, ref,
        sites[1]["third"] if p == sites[1]["pos"] else alt,
        info, gt,
    )
    for p, ref, alt, info, gt in clean
]

symbolic = clean + [(
    sites[7]["pos"], sites[7]["ref"], "<DEL>", ".", "0/0"
)]
symbolic_ref_effect = clean + [(
    sites[6]["pos"], sites[6]["ref"], "<DEL>", ".", "0/0"
)]

write_vcf("plain/clean.vcf.gz", clean)
write_vcf("plain/no_call.vcf.gz", plain_records(no_call))
write_vcf("plain/ref_alt_swap.vcf.gz", swap)
write_vcf("plain/allele_mismatch.vcf.gz", mismatch)
write_vcf(
    "plain/tied_records.vcf.gz",
    clean + [record(1, "1/1")],
)
write_vcf(
    "plain/ambiguous_contigs.vcf.gz",
    clean, contigs=("22", "chr22"),
)
write_vcf(
    "plain/symbolic_homref.vcf.gz",
    symbolic, deletion_alt=True,
)
write_vcf(
    "plain/symbolic_ref_effect.vcf.gz",
    symbolic_ref_effect, deletion_alt=True,
)
write_vcf(
    "plain/two_samples.vcf.gz",
    [record(1, "0/1")], samples=("vcf_A", "vcf_B"),
)
# Leave this one uncompressed: do not allow a conversion step to
# repair/change its deliberately malformed header.
write_vcf(
    "plain/no_contigs.vcf",
    [record(1, "0/1")], contigs=(),
)

def gvcf_records(*, star=False, gap=False, bad_block_ref=False):
    nonref = not star
    specials = {
        row[0]: row
        for row in (
            [extra_record(nonref=nonref)]
            + [
                record(i, BASE_CALLS[i], nonref=nonref)
                for i in (1, 2, 3, 4, 5, 9, 10, 11)
            ]
        )
    }

    if gap:
        # A reference-only NO-CALL: current R code loses this record.
        p7 = sites[7]["pos"]
        specials[p7] = (
            p7, sites[7]["ref"], ".", f"END={p7}", "./."
        )
        # A no-call at a record with concrete ALT: it is recognized.
        specials[sites[8]["pos"]] = record(
            8, "./.", nonref=nonref
        )

    uncovered = {sites[12]["pos"]} if gap else set()
    rows = []
    p = START
    while p <= END:
        if p in specials:
            rows.append(specials[p])
            p += 1
        elif p in uncovered:
            p += 1
        else:
            block_start = p
            while (
                p <= END
                and p not in specials
                and p not in uncovered
            ):
                p += 1
            rows.append((
                block_start,
                base_at(block_start),
                "." if star else "<NON_REF>",
                f"END={p - 1}",
                "0/0",
            ))

    if bad_block_ref:
        # Corrupt the REF at the START of a block that contains s06
        # internally. This negative control is not hg38-valid.
        p6 = sites[6]["pos"]
        k = next(
            k for k, row in enumerate(rows)
            if row[2] == "<NON_REF>"
            and row[0] < p6 <= int(row[3][4:])
        )
        old = rows[k]
        wrong = next(x for x in "ACGT" if x != old[1])
        rows[k] = (old[0], wrong, old[2], old[3], old[4])

    return rows

write_vcf(
    "gvcf/full.vcf.gz",
    gvcf_records(), gvcf=True,
)
write_vcf(
    "gvcf/gap_and_no_calls.vcf.gz",
    gvcf_records(gap=True), gvcf=True,
)
write_vcf(
    "gvcf/starling_style.vcf.gz",
    gvcf_records(star=True), gvcf=True, star=True,
)
write_vcf(
    "gvcf/bad_block_ref.vcf.gz",
    gvcf_records(bad_block_ref=True), gvcf=True,
)

def sheet(name, entries):
    (root / name).write_text(
        "".join(f"{sample}\t{vcf}\n" for sample, vcf in entries)
    )

sheets = {
    "plain_reference.tsv": [
        ("clean", "plain/clean.vcf.gz"),
        ("no_call", "plain/no_call.vcf.gz"),
    ],
    "plain_zero.tsv": [
        ("clean", "plain/clean.vcf.gz"),
        ("swap", "plain/ref_alt_swap.vcf.gz"),
    ],
    "plain_single.tsv": [("clean", "plain/clean.vcf.gz")],
    "plain_no_call.tsv": [("no_call", "plain/no_call.vcf.gz")],
    "plain_mismatch.tsv": [
        ("mismatch", "plain/allele_mismatch.vcf.gz")
    ],
    "plain_tied.tsv": [("tied", "plain/tied_records.vcf.gz")],
    "plain_ambiguous_contigs.tsv": [
        ("ambiguous_contigs", "plain/ambiguous_contigs.vcf.gz")
    ],
    "plain_symbolic.tsv": [
        ("symbolic", "plain/symbolic_homref.vcf.gz")
    ],
    "plain_symbolic_ref_effect.tsv": [
        ("symbolic_ref", "plain/symbolic_ref_effect.vcf.gz")
    ],
    "plain_two_samples.tsv": [
        ("two_samples", "plain/two_samples.vcf.gz")
    ],
    "plain_no_contigs.tsv": [
        ("no_contigs", "plain/no_contigs.vcf")
    ],
    "plain_numeric_ids.tsv": [
        ("001", "plain/clean.vcf.gz"),
        ("1", "plain/no_call.vcf.gz"),
    ],
    "gvcf_full.tsv": [("full", "gvcf/full.vcf.gz")],
    "gvcf_gap.tsv": [
        ("gap", "gvcf/gap_and_no_calls.vcf.gz")
    ],
    "gvcf_star.tsv": [
        ("star", "gvcf/starling_style.vcf.gz")
    ],
    "gvcf_bad_block.tsv": [
        ("bad_block", "gvcf/bad_block_ref.vcf.gz")
    ],
}
for name, entries in sheets.items():
    sheet(name, entries)

if args.cohort:
    entries = []
    phenotypes = []
    gt = ("0/0", "0/1", "1/1")
    for j in range(10):
        name = f"cohort_{j + 1:02}"
        calls = dict(BASE_CALLS)
        calls[1] = gt[j % 3]
        calls[3] = gt[(j // 3) % 3]
        calls[4] = gt[(j + 1) % 3]
        calls[5] = gt[(j // 2) % 3]
        rel = f"plain/{name}.vcf.gz"
        write_vcf(rel, plain_records(calls), samples=(name,))
        entries.append((name, rel))
        phenotypes.append((
            name, "control" if j < 5 else "case"
        ))
    sheet("cohort.tsv", entries)
    with (root / "phenotypes.tsv").open(
        "w", newline=""
    ) as handle:
        w = csv.writer(
            handle, delimiter="\t", lineterminator="\n"
        )
        w.writerow(("sample_id", "phenotype"))
        w.writerows(phenotypes)

print(f"Wrote fixtures to {root}")
print("Inspect sites.tsv for the actual hg38 coordinates and bases.")