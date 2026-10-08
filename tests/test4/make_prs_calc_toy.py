#!/usr/bin/env python3

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path


# Adapt these names if your score reader uses different column names.
# The values written below follow this exact column order.
SCORE_COLUMNS = [
    "chr_name",
    "chr_position",
    "effect_allele",
    "other_allele",
    "effect_weight",
    "effect_genotype",
]


parser = argparse.ArgumentParser(
    description="Generate single-sample PRS test fixtures and pipeline commands."
)
parser.add_argument("--out", default="toy_prs")
parser.add_argument(
    "--command",
    required=True,
    help=(
        "An argv-style command template. Available placeholders: "
        "{vcf}, {score}, {score_dir}, {reference}, {out}, {work}, "
        "{non_additive}, {genotype_calls}, {no_calls}."
    ),
)
args = parser.parse_args()

for tool in ("bgzip", "tabix"):
    if shutil.which(tool) is None:
        parser.error(f"Required executable not found: {tool}")

for placeholder in (
    "{non_additive}",
    "{genotype_calls}",
    "{no_calls}",
):
    if placeholder not in args.command:
        parser.error(f"--command must contain {placeholder}")

command_tokens = shlex.split(args.command)
root = Path(args.out).resolve()
root.mkdir(parents=True, exist_ok=True)
(root / "logs").mkdir(exist_ok=True)

# A real reference makes the bases at internal reference-block positions
# unambiguous. All fixtures use REF=A on contig 1.
reference = root / "toy.fa"
reference.write_text(
    ">1\n" + "\n".join("A" * min(60, 1000 - i)
                      for i in range(0, 1000, 60)) + "\n"
)
Path(str(reference) + ".fai").write_text(
    "1\t1000\t3\t60\t61\n"
)


def V(pos, alt, gt, gp, info="."):
    """Create a VCF row. None means the FORMAT field is absent."""
    keys, values = [], []

    if gt is not None:
        keys.append("GT")
        values.append(gt)

    if gp is not None:
        keys.append("GP")
        values.append(gp)

    return [
        "1", str(pos), f"v{pos}", "A", alt, "60", "PASS", info,
        ":".join(keys), ":".join(values),
    ]


def S(pos, effect="A", genotypes="AA", other=None):
    """Create a score row, always with weight 0.24."""
    if other is None:
        other = "T" if effect == "A" else "A"

    return [
        "1", str(pos), effect, other, "0.24", genotypes,
    ]


def E(hard_add, soft_add, hard_nonadd, soft_nonadd):
    """Independent expected raw totals. None means expected error."""
    return {
        ("false", "hard"): hard_add,
        ("false", "soft"): soft_add,
        ("true", "hard"): hard_nonadd,
        ("true", "soft"): soft_nonadd,
    }


fixtures = []


def add(name, rows, scores, expected, *,
        gvcf=False, gt_header=True, gp_header=True,
        zero=None, pass_kind="PASS", note=""):
    fixtures.append({
        "name": name,
        "rows": rows,
        "scores": scores,
        "expected": expected,
        "gvcf": gvcf,
        "gt_header": gt_header,
        "gp_header": gp_header,
        "zero": zero,
        "pass_kind": pass_kind,
        "note": note,
    })


# Five biallelic SNPs exercise:
# - reference and alternative effect alleles;
# - recessive, dominant and heterozygote-only genotype sets;
# - phased GT;
# - TA versus AT canonicalisation;
# - no inference of unlisted genotypes.
core = [
    (100, "0/0", "0.83,0.16,0.01"),
    (110, "0|1", "0.10,0.70,0.20"),
    (120, "1/1", "0.05,0.15,0.80"),
    (130, "1|0", "0.20,0.60,0.20"),
    (140, "0/1", "0.10,0.70,0.20"),
]
core_scores = [
    S(100, "A", "AA"),
    S(110, "T", "TA/TT"),
    S(120, "T", "TT"),
    S(130, "T", "TA"),
    S(140, "A", "AA"),
]

add(
    "complete_gp",
    [V(pos, "T", gt, gp) for pos, gt, gp in core],
    core_scores,
    E("1.6800", "1.5768", "0.9600", "0.7752"),
    note="Complete concrete SNPs; all GP vectors usable.",
)

add(
    "no_gp_header",
    [V(pos, "T", gt, None) for pos, gt, gp in core],
    core_scores,
    E("1.6800", "1.6800", "0.9600", "0.9600"),
    gp_header=False,
    note="GP undeclared: soft must fall back to GT.",
)

add(
    "mixed_gp",
    [
        V(pos, "T", gt,
          "." if pos == 100 else None if pos == 110 else gp)
        for pos, gt, gp in core
    ],
    core_scores,
    E("1.6800", "1.5960", "0.9600", "0.8400"),
    note="GP=. at 100; GP absent from FORMAT at 110; otherwise usable.",
)

nocall_scores = [
    S(100, "A", "AA/AT"),
    S(110, "A", "AA/AT"),
]

add(
    "nocalls_valid_gp",
    [
        V(100, "T", "0/.", "0.83,0.16,0.01"),
        V(110, "T", "./.", "0.83,0.16,0.01"),
    ],
    nocall_scores,
    E(None, "0.8736", None, "0.4752"),
    zero=E("0.2400", "0.8736", "0.0000", "0.4752"),
    note=(
        "GP-first policy: soft uses GP despite partial/full missing GT. "
        "Hard non-additive incomplete GT gets multiplier zero with no_calls=zero."
    ),
)

add(
    "nocalls_missing_gp",
    [
        V(100, "T", "0/.", "."),
        V(110, "T", "./.", "."),
    ],
    nocall_scores,
    E(None, None, None, None),
    zero=E("0.2400", "0.2400", "0.0000", "0.0000"),
    note="No usable GP: soft follows the same no-call rules as hard.",
)

# For A,C,G, diploid VCF GP ordering is:
# AA, AC, CC, AG, CG, GG.
# Expected A dosage = 2*.10 + .20 + .25 = .65.
# P(AA or AG) = .10 + .25 = .35.
# other_allele is intentionally unspecified for this multiallelic model.
add(
    "multiallelic",
    [V(100, "C,G", "0/2", "0.10,0.20,0.10,0.25,0.15,0.20")],
    [S(100, "A", "AA/GA", other="")],
    E("0.2400", "0.1560", "0.2400", "0.0840"),
    note="Tests multiallelic Number=G ordering and explicit genotype membership.",
)

add(
    "haploid",
    [V(100, "T", "0", "0.90,0.10")],
    [S(100)],
    E("0.2400", "0.2160", None, None),
    note="Additive supports haploid dosage; this non-additive model is diploid-only.",
)

add(
    "gp_only",
    [V(pos, "T", None, gp) for pos, gt, gp in core],
    core_scores,
    E(None, "1.5768", None, "0.7752"),
    gt_header=False,
    pass_kind="OPTIONAL",
    note="Soft totals apply only if GP-only input and preflight are supported.",
)

# A gVCF SNP record is not a reference block.
# GP ordering for A,T,<NON_REF>:
# AA, AT, TT, A/<NON_REF>, T/<NON_REF>, <NON_REF>/<NON_REF>.
add(
    "gvcf_snp_zero_symbolic_gp",
    [V(100, "T,<NON_REF>", "0/0", "0.83,0.16,0.01,0,0,0")],
    [S(100)],
    E("0.4800", "0.4368", "0.2400", "0.1992"),
    gvcf=True,
    note="Concrete SNP in gVCF: soft uses GP; symbolic states have zero mass.",
)

add(
    "gvcf_reference_block",
    [
        V(
            100, "<NON_REF>", "0/0",
            "0.83,0.16,0.01", "END=120"
        ),
    ],
    [S(100), S(110)],
    E("0.9600", "0.9600", "0.4800", "0.4800"),
    gvcf=True,
    pass_kind="POLICY",
    note=(
        "Totals assume fully called reference blocks are accepted as reference "
        "evidence. GP is not reused at 100 or internal position 110."
    ),
)

add(
    "gvcf_positive_symbolic_gp",
    [V(100, "T,<NON_REF>", "0/0", "0.82,0.16,0.01,0.01,0,0")],
    [S(100)],
    E("0.4800", None, "0.2400", None),
    gvcf=True,
    note="Soft must reject positive probability on unsupported symbolic states.",
)

add(
    "gvcf_symbolic_gt",
    [V(100, "T,<NON_REF>", "0/2", "0.83,0.16,0.01,0,0,0")],
    [S(100)],
    E(None, None, None, None),
    gvcf=True,
    note="GT selects <NON_REF>; reject even if GP gives symbolic states zero mass.",
)

for name, gp in [
    ("bad_gp_length", "0.90,0.10"),
    ("bad_gp_negative", "1.10,-0.10,0.00"),
]:
    add(
        name,
        [V(100, "T", "0/0", gp)],
        [S(100)],
        E("0.4800", None, "0.2400", None),
        note="Hard ignores GP; soft rejects malformed GP rather than falling back.",
    )


def write_tsv(path, header, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


# Write, compress and index each fixture.
manifest_rows = []

for case in fixtures:
    directory = root / "data" / case["name"]
    scores_directory = directory / "scores"
    scores_directory.mkdir(parents=True, exist_ok=True)

    suffix = ".g.vcf" if case["gvcf"] else ".vcf"
    plain_vcf = directory / (case["name"] + suffix)
    compressed_vcf = Path(str(plain_vcf) + ".gz")
    score_path = scores_directory / (case["name"] + ".tsv")

    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=1000>",
        '##ALT=<ID=NON_REF,Description="Any possible non-reference allele">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Record end position">',
    ]

    if case["gt_header"]:
        header.append(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        )

    if case["gp_header"]:
        header.append(
            '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">'
        )

    header.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTOY"
    )

    plain_vcf.write_text(
        "\n".join(header)
        + "\n"
        + "\n".join("\t".join(row) for row in case["rows"])
        + "\n"
    )

    with compressed_vcf.open("wb") as handle:
        subprocess.run(
            ["bgzip", "-c", str(plain_vcf)],
            stdout=handle,
            check=True,
        )

    subprocess.run(
        ["tabix", "-f", "-p", "vcf", str(compressed_vcf)],
        check=True,
    )

    write_tsv(score_path, SCORE_COLUMNS, case["scores"])

    case["vcf"] = compressed_vcf
    case["score"] = score_path

    manifest_rows.append([
        case["name"], str(compressed_vcf), str(score_path), case["note"]
    ])

write_tsv(
    root / "fixtures.tsv",
    ["scenario", "vcf", "score", "note"],
    manifest_rows,
)


# Generate commands and expected outcomes.
exit_codes = root / "exit_codes.tsv"
quote = shlex.quote

run_lines = [
    "#!/usr/bin/env bash",
    # Deliberately no set -e: negative tests must not stop the test suite.
    "set -uo pipefail",
    f"mkdir -p {quote(str(root / 'logs'))}",
    f"printf 'test_id\\texit_code\\n' > {quote(str(exit_codes))}",
]

expected_rows = []
command_rows = []

for case in fixtures:
    policies = [("error", case["expected"])]

    if case["zero"] is not None:
        policies.append(("zero", case["zero"]))

    for no_calls, expectations in policies:
        for non_additive in ("false", "true"):
            for genotype_calls in ("hard", "soft"):
                test_id = (
                    f"{case['name']}__nc_{no_calls}"
                    f"__na_{non_additive}__gc_{genotype_calls}"
                )

                total = expectations[(non_additive, genotype_calls)]
                outcome = "ERROR" if total is None else case["pass_kind"]

                values = {
                    "vcf": str(case["vcf"]),
                    "score": str(case["score"]),
                    "score_dir": str(case["score"].parent),
                    "reference": str(reference),
                    "out": str(root / "results" / test_id),
                    "work": str(root / "work" / test_id),
                    "non_additive": non_additive,
                    "genotype_calls": genotype_calls,
                    "no_calls": no_calls,
                }

                try:
                    argv = [
                        token.format_map(values)
                        for token in command_tokens
                    ]
                except KeyError as exc:
                    parser.error(f"Unknown command-template placeholder: {exc}")

                command = shlex.join(argv)
                log = root / "logs" / f"{test_id}.log"

                expected_rows.append([
                    test_id,
                    case["name"],
                    no_calls,
                    non_additive,
                    genotype_calls,
                    outcome,
                    "" if total is None else total,
                    case["note"],
                ])

                command_rows.append([test_id, command])

                run_lines.extend([
                    f"printf '\\n=== %s ===\\n' {quote(test_id)}",
                    f"if {command} > {quote(str(log))} 2>&1; then",
                    "  rc=0",
                    "else",
                    "  rc=$?",
                    "fi",
                    (
                        f"printf '%s\\t%s\\n' {quote(test_id)} \"$rc\""
                        f" >> {quote(str(exit_codes))}"
                    ),
                ])

write_tsv(
    root / "expected.tsv",
    [
        "test_id", "scenario", "no_calls", "non_additive",
        "genotype_calls", "expected_outcome", "expected_total", "note",
    ],
    expected_rows,
)

write_tsv(
    root / "commands.tsv",
    ["test_id", "command"],
    command_rows,
)

runner = root / "run_all.sh"
runner.write_text("\n".join(run_lines) + "\n")
runner.chmod(0o755)

print(f"Generated {len(fixtures)} fixtures and {len(command_rows)} test commands.")
print(f"Fixtures:     {root / 'fixtures.tsv'}")
print(f"Expectations: {root / 'expected.tsv'}")
print(f"Commands:     {root / 'commands.tsv'}")
print(f"Runner:       {runner}")
print("No pipeline commands have been executed.")
