#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import Counter


COLUMNS = {"chrom", "position", "ref", "alt", "end"}


def read_records(filename):
    records = Counter()

    with open(filename, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not COLUMNS.issubset(reader.fieldnames or []):
            raise ValueError(
                f"{filename} needs columns: "
                + ", ".join(sorted(COLUMNS))
            )

        for row in reader:
            key = (
                row["chrom"],
                int(row["position"]),
                row["ref"].upper(),
                row["alt"].upper(),
                row["end"] if row["end"] not in ("", None) else ".",
            )
            records[key] += 1

    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", required=True)
    parser.add_argument("--queried", required=True)
    args = parser.parse_args()

    expected = read_records(args.expected)
    queried = read_records(args.queried)

    for record, count in expected.items():
        if queried[record] < count:
            chrom, position, ref, alt, end = record
            raise ValueError(
                "Targeted bcftools query omitted a VCF record "
                f"intersecting a score position: {chrom}:{position} "
                f"REF={ref}, ALT={alt}, END={end}; "
                f"expected {count}, queried {queried[record]}. "
                "Refusing to treat observed VCF evidence as missing."
            )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
        