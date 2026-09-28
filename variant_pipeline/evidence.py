"""Evidence statistics: strand bias via Fisher's exact test and confidence.

A candidate is never decided by raw occurrence count. The confidence score
combines alt fraction, base quality, strand balance, low-quality evidence,
mate overlap, coverage and conflict evidence with transparent weights.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for a 2x2 table [[a,b],[c,d]]."""
    n = a + b + c + d
    r1 = a + b
    r2 = c + d
    col1 = a + c
    if n == 0:
        return 1.0
    lo = max(0, col1 - r2)
    hi = min(r1, col1)

    def prob(x: int) -> float:
        return math.exp(
            _log_comb(r1, x)
            + _log_comb(r2, col1 - x)
            - _log_comb(n, col1)
        )

    observed = prob(a)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1) if prob(x) <= observed + 1e-12))


SCORE_WEIGHTS = {
    "alt_fraction": 0.40,
    "base_quality": 0.20,
    "strand_balance": 0.15,
    "coverage": 0.10,
    "quality_consistency": 0.10,
    "conflict": 0.05,
}


@dataclass
class CandidateEvidence:
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    kind: str
    # supporting evidence
    alt_fragments: int
    alt_forward: int
    alt_reverse: int
    alt_mean_base_quality: float
    low_quality_alt_fragments: int
    # opposing evidence
    ref_fragments: int
    other_alt_fragments: int
    conflicting_overlap_fragments: int
    strand_bias_p: float
    # depth
    depth: int                 # fragments spanning the variant anchor
    passing_depth: int         # fragments eligible at the anchor (MAPQ/base)
    alt_fraction: float
    informative_depth: int     # alt + ref (decidable fragments)
    # score components 0..1 and overall
    components: dict[str, float] = field(default_factory=dict)
    confidence: float = 0.0
    flags: list[str] = field(default_factory=list)
    status: str = "LOW_COVERAGE"  # PASS | SUSPECT | REJECTED | LOW_COVERAGE

    def supporting_summary(self) -> dict:
        return {
            "alt_fragments": self.alt_fragments,
            "forward": self.alt_forward,
            "reverse": self.alt_reverse,
            "mean_base_quality": round(self.alt_mean_base_quality, 2),
            "low_quality_alt_fragments": self.low_quality_alt_fragments,
            "alt_fraction_of_informative": round(self.alt_fraction, 4),
        }

    def opposing_summary(self) -> dict:
        return {
            "reference_fragments": self.ref_fragments,
            "other_allele_fragments": self.other_alt_fragments,
            "overlap_conflict_fragments": self.conflicting_overlap_fragments,
            "strand_bias_fisher_p": (
                round(self.strand_bias_p, 6) if self.strand_bias_p is not None else None
            ),
        }


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def score_candidate(
    ev: CandidateEvidence,
    min_alt_count: int,
    min_alt_fraction: float,
    min_coverage: int,
    strand_bias_p_threshold: float,
    low_quality_fraction_threshold: float,
) -> CandidateEvidence:
    components: dict[str, float] = {}
    flags: list[str] = []

    # Alt fraction relative to fragments that actually show REF or this ALT.
    components["alt_fraction"] = _clamp01(ev.alt_fraction)

    # Base quality: Q40 -> ~1.0, Q20 -> ~0.68.
    q = ev.alt_mean_base_quality
    components["base_quality"] = _clamp01(1.0 - 10.0 ** (-q / 10.0)) if q > 0 else 0.0

    # Strand balance: both directions support -> 1.0; one direction -> 0.0.
    if ev.alt_forward + ev.alt_reverse > 0:
        balance = min(ev.alt_forward, ev.alt_reverse) / max(
            ev.alt_forward, ev.alt_reverse
        )
    else:
        balance = 0.0
    components["strand_balance"] = balance
    if ev.alt_fragments > 0 and (
        ev.alt_forward == 0 or ev.alt_reverse == 0
    ):
        flags.append("ONE_STRAND_ONLY")

    # Coverage saturates at 4x the calling minimum.
    components["coverage"] = _clamp01(ev.depth / max(1, 4 * min_coverage))

    low_q_fraction = (
        ev.low_quality_alt_fragments / ev.alt_fragments if ev.alt_fragments else 0.0
    )
    components["quality_consistency"] = 1.0 - low_q_fraction
    if low_q_fraction > low_quality_fraction_threshold:
        flags.append("LOW_QUALITY_EVIDENCE")

    # Conflict/opposition penalty.
    opposition = ev.other_alt_fragments + ev.conflicting_overlap_fragments
    if ev.alt_fragments + opposition > 0:
        components["conflict"] = 1.0 - opposition / (ev.alt_fragments + opposition)
    else:
        components["conflict"] = 0.0

    if ev.strand_bias_p < strand_bias_p_threshold and ev.alt_fragments > 0:
        flags.append("STRAND_BIAS")
    if ev.conflicting_overlap_fragments > 0:
        flags.append("MATE_OVERLAP_CONFLICT")
    if ev.passing_depth < min_coverage or ev.depth < min_coverage:
        flags.append("INSUFFICIENT_COVERAGE")

    confidence = sum(
        SCORE_WEIGHTS[name] * components[name] for name in SCORE_WEIGHTS
    )
    ev.components = {k: round(v, 4) for k, v in components.items()}
    ev.confidence = round(confidence, 4)
    ev.flags = flags

    if ev.depth < min_coverage or ev.passing_depth < min_coverage:
        ev.status = "LOW_COVERAGE"
    elif ev.alt_fragments < min_alt_count or ev.alt_fraction < min_alt_fraction:
        ev.status = "REJECTED"
    elif {
        "STRAND_BIAS",
        "LOW_QUALITY_EVIDENCE",
        "MATE_OVERLAP_CONFLICT",
    } & set(flags):
        ev.status = "SUSPECT"
    else:
        ev.status = "PASS"
    return ev
