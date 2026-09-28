from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from .caller import Candidate
from .models import SampleConfig


def build_report(config: SampleConfig, candidates: list[Candidate]) -> dict[str, object]:
    return {
        "research_only": True,
        "disclaimer": (
            "Research workflow output only; not validated for diagnosis, treatment, "
            "or other medical decisions."
        ),
        "sample": {
            "sample_id": config.sample_id,
            "reference_name": config.reference_name,
            "reference_version": config.reference_version,
        },
        "parameters": asdict(config.parameters),
        "candidates": [asdict(candidate) for candidate in candidates],
    }


def write_json(report: dict[str, object], path: str | Path) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def write_csv(candidates: list[Candidate], path: str | Path) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "variant_id",
        "chromosome",
        "position",
        "reference_allele",
        "alternate_allele",
        "variant_type",
        "status",
        "score",
        "support_templates",
        "reference_templates",
        "other_allele_templates",
        "conflicting_templates",
        "low_quality_templates",
        "template_depth",
        "allele_frequency",
        "strand_bias_fisher_p",
        "coverage_note",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for candidate in candidates:
            coverage = candidate.evidence["coverage"]
            supports = candidate.evidence["supports"]
            opposes = candidate.evidence["opposes"]
            writer.writerow(
                {
                    "variant_id": candidate.variant_id,
                    "chromosome": candidate.chromosome,
                    "position": candidate.position,
                    "reference_allele": candidate.reference_allele,
                    "alternate_allele": candidate.alternate_allele,
                    "variant_type": candidate.variant_type,
                    "status": candidate.status,
                    "score": candidate.score,
                    "support_templates": supports["templates"],
                    "reference_templates": opposes["reference_templates"],
                    "other_allele_templates": opposes["other_allele_templates"],
                    "conflicting_templates": opposes["conflicting_templates"],
                    "low_quality_templates": candidate.evidence["low_quality_templates"],
                    "template_depth": coverage["template_depth"],
                    "allele_frequency": coverage["allele_frequency"],
                    "strand_bias_fisher_p": coverage["strand_bias_fisher_p"],
                    "coverage_note": candidate.evidence["coverage_note"] or "",
                }
            )
