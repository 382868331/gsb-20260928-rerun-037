from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import comb

from .alignments import Alignment, observations_and_coverage
from .models import Observation, Parameters


@dataclass
class Candidate:
    variant_id: str
    chromosome: str
    position: int
    reference_allele: str
    alternate_allele: str
    variant_type: str
    score: float
    status: str
    evidence: dict[str, object]


def call_candidates(
    alignments: list[Alignment],
    reference: str,
    parameters: Parameters,
) -> list[Candidate]:
    observations, base_quality, orientation, intervals = observations_and_coverage(
        alignments, reference
    )
    coverage = _template_coverage(intervals, len(reference))
    groups: dict[tuple[str, int, str, str, str], list[Observation]] = defaultdict(list)
    for observation in observations:
        key = (
            observation.variant_id,
            observation.chromosome,
            observation.position,
            observation.ref,
            observation.alt,
        )
        groups[key].append(observation)

    candidates: list[Candidate] = []
    for (variant_id, chromosome, position, ref, alt), grouped in sorted(
        groups.items(), key=lambda item: (item[0][2], item[0][4])
    ):
        variant_type = grouped[0].variant_type
        allele_span = max(len(ref), len(alt))
        coverage_span = 1 if variant_type in {"insertion", "deletion"} else allele_span
        anchor0 = position - 1
        depth_templates = set()
        covered_positions = 0
        covered_bases = 0
        for offset in range(coverage_span):
            ref0 = anchor0 + offset
            if 0 <= ref0 < len(coverage):
                depth_templates.update(coverage[ref0])
                if coverage[ref0]:
                    covered_positions += 1
                    covered_bases += len(coverage[ref0])
        depth = len(depth_templates)
        mean_depth = covered_bases / coverage_span if coverage_span else 0.0

        by_template: dict[str, list[Observation]] = defaultdict(list)
        for observation in grouped:
            by_template[observation.template_id].append(observation)

        alt_templates: set[str] = set()
        low_quality_alt_templates: set[str] = set()
        forward_support = reverse_support = 0
        for template_id, template_observations in by_template.items():
            best = max(template_observations, key=lambda item: item.base_quality)
            if (
                best.base_quality >= parameters.min_base_quality
                and best.mapping_quality >= parameters.min_mapping_quality
            ):
                alt_templates.add(template_id)
                if best.is_reverse:
                    reverse_support += 1
                else:
                    forward_support += 1
            else:
                low_quality_alt_templates.add(template_id)

        reference_templates: set[str] = set()
        low_quality_reference_templates: set[str] = set()
        forward_reference = reverse_reference = 0
        anchor_quality_by_template = base_quality.get(anchor0, {})
        for template_id in depth_templates:
            if template_id in by_template:
                continue
            if anchor_quality_by_template.get(template_id, -1) < parameters.min_base_quality:
                low_quality_reference_templates.add(template_id)
            else:
                reference_templates.add(template_id)
                anchor_orientation = orientation.get(anchor0, {}).get(template_id)
                if anchor_orientation is not None and not anchor_orientation[1]:
                    forward_reference += 1
                else:
                    reverse_reference += 1

        other_allele_templates: set[str] = set()
        conflicting_templates: set[str] = set()
        for other in observations:
            if (
                other.chromosome == chromosome
                and other.position == position
                and (other.ref, other.alt) != (ref, alt)
            ):
                if other.template_id in alt_templates:
                    conflicting_templates.add(other.template_id)
                else:
                    other_allele_templates.add(other.template_id)

        strand_bias_p = fisher_exact_two_sided(
            forward_support,
            reverse_support,
            forward_reference,
            reverse_reference,
        )
        support_count = len(alt_templates)
        ref_count = len(reference_templates)
        other_count = len(other_allele_templates)
        conflict_count = len(conflicting_templates)
        low_quality_count = len(low_quality_alt_templates | low_quality_reference_templates)
        denominator = support_count + ref_count + other_count
        allele_frequency = support_count / denominator if denominator else 0.0
        quality_observations = [item.base_quality for item in grouped]
        mean_alt_quality = (
            sum(quality_observations) / len(quality_observations)
            if quality_observations
            else 0.0
        )

        vaf_component = min(35.0, allele_frequency * 70.0)
        support_component = min(25.0, support_count * 5.0)
        quality_component = min(15.0, mean_alt_quality / 40.0 * 15.0)
        depth_component = min(15.0, depth / max(parameters.min_depth, 1) * 15.0)
        strand_component = max(0.0, 10.0 - (1.0 - strand_bias_p) * 10.0)
        penalty = conflict_count * 10.0 + other_count * 3.0
        score = round(
            max(
                0.0,
                vaf_component
                + support_component
                + quality_component
                + depth_component
                + strand_component
                - penalty,
            ),
            2,
        )

        coverage_note = None
        if covered_positions < coverage_span:
            coverage_note = "Allele interval is not completely covered by aligned reads"
        elif mean_depth < parameters.min_depth:
            coverage_note = (
                f"Mean template depth {mean_depth:.1f} is below minimum "
                f"{parameters.min_depth}"
            )

        hard_pass = (
            support_count >= parameters.min_alt_templates
            and depth >= parameters.min_depth
            and coverage_note is None
            and strand_bias_p >= parameters.min_strand_bias_p
        )
        status = "candidate" if score >= parameters.min_pass_score and hard_pass else "review"

        candidates.append(
            Candidate(
                variant_id=variant_id,
                chromosome=chromosome,
                position=position,
                reference_allele=ref,
                alternate_allele=alt,
                variant_type=variant_type,
                score=score,
                status=status,
                evidence={
                    "supports": {
                        "templates": support_count,
                        "forward_templates": forward_support,
                        "reverse_templates": reverse_support,
                        "mean_base_quality": round(mean_alt_quality, 2),
                    },
                    "opposes": {
                        "reference_templates": ref_count,
                        "other_allele_templates": other_count,
                        "conflicting_templates": conflict_count,
                    },
                    "low_quality_templates": low_quality_count,
                    "coverage": {
                        "template_depth": depth,
                        "mean_allele_depth": round(mean_depth, 2),
                        "covered_allele_positions": covered_positions,
                        "required_allele_positions": coverage_span,
                        "allele_frequency": round(allele_frequency, 4),
                        "strand_bias_fisher_p": round(strand_bias_p, 4),
                    },
                    "coverage_note": coverage_note,
                },
            )
        )
    return candidates


def _template_coverage(
    intervals: list[tuple[str, int, int]], reference_length: int
) -> list[set[str]]:
    coverage = [set() for _ in range(reference_length)]
    seen: set[tuple[str, int, int]] = set()
    for template_id, start, end in intervals:
        interval = (template_id, start, end)
        if interval in seen:
            continue
        seen.add(interval)
        for ref0 in range(max(0, start), min(reference_length, end)):
            coverage[ref0].add(template_id)
    return coverage


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    row1 = a + b
    col1 = a + c
    total = a + b + c + d
    if total == 0:
        return 1.0

    def probability(x: int) -> float:
        return (
            comb(row1, x)
            * comb(total - row1, col1 - x)
            / comb(total, col1)
        )

    observed = probability(a)
    low = max(0, col1 - (total - row1))
    high = min(row1, col1)
    return min(
        1.0,
        sum(probability(x) for x in range(low, high + 1) if probability(x) <= observed + 1e-12),
    )
