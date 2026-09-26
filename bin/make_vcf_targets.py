#!/usr/bin/env python3

import argparse
import csv
import re
import sys
from collections import defaultdict


def canonical_contig(value: str) -> str:
    value = value.strip()

    if value.lower().startswith("chr"):
        value = value[3:]

    value = value.upper()

    aliases = {
        "23": "X",
        "24": "Y",
        "25": "MT",
        "M": "MT",
        "MTDNA": "MT",
    }

    return aliases.get(value, value)


def read_score_positions(scorefile):
    positions = defaultdict(set)

    with open(scorefile, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"chrom", "position"}
        missing = required.difference(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                "Normalised score file is missing column(s): "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            chrom = canonical_contig(row["chrom"])
            position = int(row["position"])
            positions[chrom].add(position)

    if not positions:
        raise RuntimeError("The normalised score file has no positions")

    return positions


def read_header_contigs(header_file):
    contigs = []

    pattern = re.compile(r"^##contig=<ID=([^,>]+)")

    with open(
        header_file,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            match = pattern.match(line)

            if match:
                contig = match.group(1).strip().strip('"')

                if contig not in contigs:
                    contigs.append(contig)

    if not contigs:
        raise RuntimeError(
            "The VCF header does not contain ##contig declarations. "
            "Contig declarations are needed to map score chromosomes to "
            "the chromosome naming scheme used by the VCF."
        )

    return contigs


def select_vcf_contigs(score_contig, vcf_contigs):
    """Return every declared VCF spelling of this score chromosome."""
    return [
        contig
        for contig in vcf_contigs
        if canonical_contig(contig) == score_contig
    ]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Map normalised score chromosomes to VCF contig names and "
            "write a bcftools target file."
        )
    )

    parser.add_argument("--scorefile", required=True)
    parser.add_argument("--vcf-header", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)

    args = parser.parse_args()

    score_positions = read_score_positions(args.scorefile)
    vcf_contigs = read_header_contigs(args.vcf_header)

    primary_chromosomes = {
        *(str(number) for number in range(1, 23)),
        "X",
        "Y",
        "MT",
    }

    primary_contigs = [
        contig
        for contig in vcf_contigs
        if canonical_contig(contig) in primary_chromosomes
    ]

    prefixed = [
        contig.lower().startswith("chr")
        for contig in primary_contigs
    ]

    if any(prefixed) and any(not value for value in prefixed):
        print(
            "WARNING: VCF header mixes primary chromosome names "
            "with and without 'chr'. Both spellings will be queried "
            "and conflicting records at equivalent positions will "
            "be rejected.",
            file=sys.stderr,
        )

    contig_order = {
        contig: index
        for index, contig in enumerate(vcf_contigs)
    }

    targets = []
    mapping_report = []

    for score_contig in sorted(score_positions):
        matched_contigs = select_vcf_contigs(
            score_contig,
            vcf_contigs,
        )

        n_targets = len(score_positions[score_contig])

        mapping_report.append(
            {
                "score_chrom": score_contig,
                "vcf_chrom": ",".join(matched_contigs),
                "n_targets": n_targets,
                "mapped": "TRUE" if matched_contigs else "FALSE",
            }
        )

        if not matched_contigs:
            print(
                f"WARNING: score chromosome '{score_contig}' is not "
                f"declared in the VCF header",
                file=sys.stderr,
            )
            continue

        for vcf_contig in matched_contigs:
            for position in score_positions[score_contig]:
                targets.append(
                    (vcf_contig, position, position)
                )

    if not targets:
        raise RuntimeError(
            "None of the score chromosomes could be mapped to VCF contigs"
        )

    targets.sort(
        key=lambda row: (
            contig_order.get(row[0], len(contig_order)),
            row[1],
        )
    )

    with open(args.output, "w", encoding="utf-8") as handle:
        for chrom, start, end in targets:
            handle.write(f"{chrom}\t{start}\t{end}\n")

    with open(
        args.report,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "score_chrom",
                "vcf_chrom",
                "n_targets",
                "mapped",
            ],
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(mapping_report)

    print(
        f"Wrote {len(targets)} VCF target interval(s) for "
        f"{sum(len(values) for values in score_positions.values())} "
        f"unique score position(s)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
