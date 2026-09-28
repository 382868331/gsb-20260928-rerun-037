"""Typed, file-backed run configuration.

Sample identity, reference version and analysis parameters live in a JSON
config file rather than in ad-hoc scripts, so downstream annotation stages
can consume the same run metadata and variant identity scheme.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path


@dataclass(frozen=True)
class SampleConfig:
    sample_id: str
    read_group: str | None = None
    description: str | None = None


@dataclass(frozen=True)
class ReferenceConfig:
    fasta_path: str
    name: str
    assembly: str  # reference version, must never be implicit
    description: str | None = None


@dataclass(frozen=True)
class Parameters:
    min_base_quality: int = 20          # Phred, inclusive
    min_mapping_quality: int = 1        # MAPQ, inclusive
    min_coverage_for_calling: int = 5   # passing fragment depth
    min_alt_count: int = 2              # passing alt observations
    min_alt_fraction: float = 0.2       # alt / informative observations
    strand_bias_p_threshold: float = 0.01
    low_quality_fraction_threshold: float = 0.2
    low_coverage_depth: int = 5         # separate coverage-report cutoff


@dataclass(frozen=True)
class RunConfig:
    sample: SampleConfig
    reference: ReferenceConfig
    parameters: Parameters = field(default_factory=Parameters)
    reads_path: str | None = None

    @staticmethod
    def from_dict(data: dict) -> "RunConfig":
        if "sample" not in data or "reference" not in data:
            raise ValueError("config must contain 'sample' and 'reference' sections")
        sample = SampleConfig(**_filter_kwargs(SampleConfig, data["sample"]))
        reference = ReferenceConfig(
            **_filter_kwargs(ReferenceConfig, data["reference"])
        )
        params = Parameters(
            **_filter_kwargs(Parameters, data.get("parameters", {}))
        )
        return RunConfig(
            sample=sample,
            reference=reference,
            parameters=params,
            reads_path=data.get("reads_path"),
        )

    @staticmethod
    def load(path: str | Path) -> "RunConfig":
        path = Path(path)
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        cfg = RunConfig.from_dict(data)
        return _resolve_paths(cfg, path.parent)

    def to_dict(self) -> dict:
        return asdict(self)


def _filter_kwargs(cls, data: dict) -> dict:
    allowed = {f.name for f in fields(cls)}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(f"unknown config keys for {cls.__name__}: {sorted(unknown)}")
    return {k: data[k] for k in data if k in allowed}


def _resolve_paths(cfg: RunConfig, base_dir: Path) -> RunConfig:
    ref_path = Path(cfg.reference.fasta_path)
    if not ref_path.is_absolute():
        ref_path = base_dir / ref_path
    reads_path = cfg.reads_path
    if reads_path is not None:
        rp = Path(reads_path)
        if not rp.is_absolute():
            reads_path = str((base_dir / rp).resolve())
    reference = ReferenceConfig(
        fasta_path=str(ref_path.resolve()),
        name=cfg.reference.name,
        assembly=cfg.reference.assembly,
        description=cfg.reference.description,
    )
    return RunConfig(
        sample=cfg.sample,
        reference=reference,
        parameters=cfg.parameters,
        reads_path=reads_path,
    )
