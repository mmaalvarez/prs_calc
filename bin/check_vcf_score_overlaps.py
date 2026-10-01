#!/usr/bin/env python3

"""Reject undeclared INFO/END and non-SNV records spanning score SNPs."""

import argparse
import csv
import re
import sys
from bisect import bisect_left
from collections import defaultdict


def canonical_contig(value):
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]

    value = value.upper()
    return {
        "23": "X",
        "24": "Y",
        "25": "MT",
        "M": "MT",
        "MTDNA": "MT",
    }.get(value, value)


def read_score_positions(filename):
    positions = defaultdict(set)

    with open(filename, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"chrom", "position"}.issubset(reader.fieldnames or []):
            raise ValueError(
                "Normalised score file needs chrom and position columns"
            )

        for row in reader:
            chrom = canonical_contig(row["chrom"])
            position = int(row["position"])
            if position < 1:
                raise ValueError("Score positions must be positive")
            positions[chrom].add(position)

    if not positions:
        raise ValueError("Normalised score file contains no positions")

    return {
        chrom: sorted(values)
        for chrom, values in positions.items()
    }


def record_end(info, declared, chrom, position):
    end_fields = [
        field
        for field in info.split(";")
        if field.partition("=")[0] == "END"
    ]

    if not end_fields:
        return None

    if not declared:
        raise ValueError(
            f"VCF record {chrom}:{position} uses INFO/END, but the "
            "VCF header has no ##INFO=<ID=END,...> declaration"
        )

    if len(end_fields) != 1 or not end_fields[0].startswith("END="):
        raise ValueError(
            f"Invalid INFO/END field at {chrom}:{position}"
        )

    value = end_fields[0][4:]
    if value == ".":
        return None

    if not re.fullmatch(r"[0-9]+", value):
        raise ValueError(
            f"Invalid INFO/END={value!r} at {chrom}:{position}"
        )

    end = int(value)
    if end < position or end > 2_147_483_647:
        raise ValueError(
            f"INFO/END={end} is out of range at {chrom}:{position}"
        )

    return end


def is_reference_only_alt(alt):
    return alt in {".", "<NON_REF>", "<*>"}


def is_snv_compatible(ref, alt, position, end):
    if len(ref) != 1 or ref not in "ACGT":
        return False

    # An extended, non-reference-only record is not a single-site SNP.
    if end is not None and end > position and not is_reference_only_alt(alt):
        return False

    if alt == ".":
        return True

    # <*> and <NON_REF> are gVCF placeholders. <DEL>, <INS>, and
    # bare * (a spanning-deletion allele) are not SNP alleles.
    return all(
        allele in {"<NON_REF>", "<*>"}
        or (len(allele) == 1 and allele in "ACGT")
        for allele in alt.split(",")
    )


def split_header_fields(body):
    """Split VCF metadata fields without splitting quoted descriptions."""
    fields = []
    field = []
    in_quotes = False
    escaped = False

    for char in body:
        if escaped:
            field.append(char)
            escaped = False
        elif char == "\\" and in_quotes:
            field.append(char)
            escaped = True
        elif char == '"':
            field.append(char)
            in_quotes = not in_quotes
        elif char == "," and not in_quotes:
            fields.append("".join(field).strip())
            field = []
        else:
            field.append(char)

    if in_quotes:
        raise ValueError("Unterminated quoted VCF header field")

    fields.append("".join(field).strip())
    return fields


def header_declares_end(filename):
    declarations = []

    with open(filename, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith("##INFO=<"):
                continue

            if not line.endswith(">"):
                raise ValueError(f"Malformed INFO declaration: {line}")

            body = line[len("##INFO=<"):-1]
            attributes = {}

            for field in split_header_fields(body):
                key, separator, value = field.partition("=")
                if not separator:
                    raise ValueError(f"Malformed INFO declaration: {line}")
                key = key.strip()
                if key in attributes:
                    raise ValueError(
                        f"Duplicate {key} in INFO declaration: {line}"
                    )
                attributes[key] = value.strip()

            if attributes.get("ID") == "END":
                declarations.append(attributes)

    if len(declarations) > 1:
        raise ValueError("VCF header declares INFO/END more than once")

    if not declarations:
        return False

    declaration = declarations[0]
    if (
        declaration.get("Number") != "1"
        or declaration.get("Type") != "Integer"
    ):
        raise ValueError(
            "INFO/END must be declared with Number=1,Type=Integer"
        )

    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scorefile", required=True)
    parser.add_argument("--vcf-header", required=True)
    parser.add_argument("--format-report", required=True)
    parser.add_argument("--expected-records", required=True)
    args = parser.parse_args()

    positions = read_score_positions(args.scorefile)
    end_declared = header_declares_end(args.vcf_header)
    has_gvcf_evidence = False

    with open(
        args.expected_records, "w", encoding="utf-8", newline=""
    ) as expected_handle:
        writer = csv.writer(expected_handle, delimiter="\t")
        writer.writerow(["chrom", "position", "ref", "alt", "end"])

        # stdin is the *entire* VCF record stream.
        for line_number, line in enumerate(sys.stdin, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise ValueError(
                    f"Malformed VCF record on streamed line {line_number}"
                )

            vcf_chrom = fields[0]
            position = int(fields[1])
            if position < 1:
                raise ValueError(
                    f"Invalid VCF position at {vcf_chrom}:{position}"
                )

            ref = fields[3].upper()
            alt = fields[4].upper()
            end = record_end(
                fields[7], end_declared, vcf_chrom, position
            )

            # An ALT header declaration alone is not evidence that any
            # record actually uses a gVCF marker.
            if (
                any(
                    allele in {"<NON_REF>", "<*>"}
                    for allele in alt.split(",")
                )
                or (
                    alt == "."
                    and end is not None
                    and len(ref) == 1
                    and ref in "ACGTN"
                )
            ):
                has_gvcf_evidence = True

            score_chrom = canonical_contig(vcf_chrom)
            chromosome_positions = positions.get(score_chrom)
            if chromosome_positions is None:
                continue

            compatible = is_snv_compatible(
                ref, alt, position, end
            )
            span_end = max(
                position + max(len(ref), 1) - 1,
                end if end is not None else position,
            )

            candidate = bisect_left(
                chromosome_positions, position
            )
            intersects_score = (
                candidate < len(chromosome_positions)
                and chromosome_positions[candidate] <= span_end
            )
            if not intersects_score:
                continue

            if not compatible:
                score_position = chromosome_positions[candidate]
                raise ValueError(
                    f"Score SNP {score_chrom}:{score_position} overlaps "
                    f"non-SNV VCF record {vcf_chrom}:{position}-{span_end} "
                    f"(REF={ref}, ALT={alt}, "
                    f"INFO/END={end if end is not None else '.'}). "
                    "Cannot score it as a SNP or treat it as absent."
                )

            # Write once per intersecting VCF record, not once per
            # score position covered by a long hom-ref block.
            writer.writerow(
                [
                    vcf_chrom,
                    position,
                    ref,
                    alt,
                    end if end is not None else ".",
                ]
            )

    with open(
        args.format_report, "w", encoding="utf-8"
    ) as handle:
        handle.write(
            ("gvcf" if has_gvcf_evidence else "ambiguous")
            + "\t"
            + ("true" if end_declared else "false")
            + "\n"
        )
        

if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)