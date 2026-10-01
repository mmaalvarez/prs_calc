#!/usr/bin/env python3

import argparse
import csv
import shutil
import sys


ID_HEADERS = {"sample", "sample_id", "id"}
PHENOTYPE_HEADERS = {
    "phenotype", "status", "label", "case_control", "case-control"
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--sample-ids", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.sample_ids, encoding="utf-8") as handle:
        sample_ids = {
            line.strip()
            for line in handle
            if line.strip()
        }

    rows = []
    with open(
        args.phenotypes, encoding="utf-8", newline=""
    ) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or (len(row) == 1 and not row[0].strip()):
                continue
            if row[0].lstrip().startswith("#"):
                continue
            if len(row) != 2:
                raise ValueError(
                    "Phenotype file must have exactly two TSV columns"
                )
            rows.append([value.strip() for value in row])

    if not rows:
        raise ValueError("Phenotype file contains no records")

    if (
        rows[0][0].lower() in ID_HEADERS
        and rows[0][1].lower() in PHENOTYPE_HEADERS
    ):
        rows.pop(0)

    if not rows:
        raise ValueError("Phenotype file contains only a header")

    phenotypes = {}
    for row_number, (sample_id, label) in enumerate(rows, start=1):
        label = label.lower()

        if not sample_id or sample_id == "NA":
            raise ValueError(
                f"Missing or reserved phenotype sample ID at row {row_number}"
            )
        if label not in {"case", "control"}:
            raise ValueError(
                f"Invalid phenotype {label!r} at row {row_number}"
            )
        if sample_id in phenotypes:
            raise ValueError(
                f"Duplicate phenotype sample ID: {sample_id}"
            )

        phenotypes[sample_id] = label

    missing = sample_ids - phenotypes.keys()
    if missing:
        raise ValueError(
            "Missing phenotypes for pipeline sample ID(s): "
            + ", ".join(sorted(missing)[:20])
        )

    scored_labels = {
        phenotypes[sample_id]
        for sample_id in sample_ids
    }
    if scored_labels != {"case", "control"}:
        raise ValueError(
            "Cohort plotting requires both a case and a control "
            "among the scored samples"
        )

    extra = phenotypes.keys() - sample_ids
    if extra:
        print(
            "WARNING: Ignoring phenotype ID(s) without a scored sample: "
            + ", ".join(sorted(extra)[:20]),
            file=sys.stderr,
        )

    shutil.copyfile(args.phenotypes, args.output)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
