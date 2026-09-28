"""End-to-end orchestration: config -> references -> fragments -> calls."""

from __future__ import annotations

from pathlib import Path

from .caller import CallResult, call_variants
from .config import RunConfig
from .fasta import read_fasta
from .samio import build_fragments, read_sam


def run(config: RunConfig, reads_path: str | None = None) -> CallResult:
    reads = reads_path or config.reads_path
    if not reads:
        raise ValueError("no reads path: set reads_path in config or pass --reads")
    references = read_fasta(config.reference.fasta_path)
    if config.reference.name not in references:
        raise ValueError(
            f"reference name {config.reference.name!r} not found in FASTA"
        )
    _header_lengths, records = read_sam(Path(reads))
    fragments = build_fragments(
        records,
        references,
        min_base_quality=config.parameters.min_base_quality,
        min_mapping_quality=config.parameters.min_mapping_quality,
    )
    return call_variants(fragments, references, config.parameters)
