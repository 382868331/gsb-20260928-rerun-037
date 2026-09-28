from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class EvidencePosition(str, Enum):
    SUPPORTS = "supports"
    OPPOSES_REFERENCE = "opposes_reference"
    OPPOSES_OTHER_ALLELE = "opposes_other_allele"
    LOW_QUALITY = "low_quality"
    CONFLICTING_TEMPLATE = "conflicting_template"


@dataclass(frozen=True)
class Parameters:
    min_base_quality: int = 20
    min_mapping_quality: int = 20
    min_depth: int = 8
    min_alt_templates: int = 3
    min_pass_score: float = 70.0
    min_strand_bias_p: float = 0.10
    skip_duplicates: bool = True


@dataclass(frozen=True)
class SampleConfig:
    sample_id: str
    reference_name: str
    reference_version: str
    reference_fasta: str
    reads_sam: str
    parameters: Parameters = field(default_factory=Parameters)


@dataclass(frozen=True)
class Observation:
    template_id: str
    chromosome: str
    position: int
    ref: str
    alt: str
    variant_type: str
    is_reverse: bool
    base_quality: int
    mapping_quality: int
    alleles_in_template: tuple[str, ...] = ()

    @property
    def variant_id(self) -> str:
        return f"{self.chromosome}:{self.position}:{self.ref}>{self.alt}"
