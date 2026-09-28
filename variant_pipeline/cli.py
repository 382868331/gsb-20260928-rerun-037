from __future__ import annotations

import argparse

from .alignments import parse_sam
from .caller import call_candidates
from .config import load_manifest
from .io import read_fasta
from .report import build_report, write_csv, write_json


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Review candidate variants and supporting/opposing evidence."
    )
    parser.add_argument("manifest", help="JSON manifest containing sample, reference, and inputs")
    parser.add_argument("--json-out", required=True, help="Path for detailed JSON report")
    parser.add_argument("--csv-out", help="Optional path for tabular CSV report")
    args = parser.parse_args()

    config = load_manifest(args.manifest)
    references = read_fasta(config.reference_fasta)
    if len(references) != 1:
        raise ValueError("This foundation workflow expects exactly one FASTA reference sequence")
    reference = next(iter(references.values()))
    alignments = parse_sam(config.reads_sam, config.parameters)
    candidates = call_candidates(alignments, reference, config.parameters)
    report = build_report(config, candidates)
    write_json(report, args.json_out)
    if args.csv_out:
        write_csv(candidates, args.csv_out)


if __name__ == "__main__":
    main()
