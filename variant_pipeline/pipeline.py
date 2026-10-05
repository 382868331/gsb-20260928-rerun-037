"""End-to-end orchestration: config -> references -> fragments -> calls."""

from __future__ import annotations

from pathlib import Path

from .association import load_association_relations
from .combination import (
    constraints_from_association_doc,
    load_combination_constraints,
)
from .store import sha256_file

from .caller import CallResult, call_variants
from .config import RunConfig
from .fasta import read_fasta
from .samio import build_fragments, read_sam


def load_fragments(config: RunConfig, reads_path: str | None = None):
    reads = reads_path or config.reads_path
    if not reads:
        raise ValueError("no reads path: set reads_path in config or pass --reads")
    references = read_fasta(config.reference.fasta_path)
    if config.reference.name not in references:
        raise ValueError(
            f"reference name {config.reference.name!r} not found in FASTA"
        )
    _header_lengths, records = read_sam(Path(reads))
    source_id = sha256_file(reads)[:16]
    fragments = build_fragments(
        records,
        references,
        min_base_quality=config.parameters.min_base_quality,
        min_mapping_quality=config.parameters.min_mapping_quality,
        source_id=source_id,
        source_name=Path(reads).name,
        batch_id=source_id,
    )
    return fragments, references


def run(config: RunConfig, reads_path: str | None = None,
        associations_path: str | None = None) -> CallResult:
    fragments, references = load_fragments(config, reads_path)
    relations = load_association_relations(
        associations_path or config.associations_path)
    constraints = _load_constraints(config, associations_path)
    return call_variants(fragments, references, config.parameters, relations,
                         constraints)


def _load_constraints(config: RunConfig, associations_path: str | None):
    assoc_doc = associations_path or config.associations_path
    constraints = constraints_from_association_doc(assoc_doc)
    constraints += load_combination_constraints(
        config.combination_constraints_path)
    # De-duplicate by constraint id while staying order independent.
    unique = {c.constraint_id: c for c in constraints}
    return [unique[cid] for cid in sorted(unique)]


def run_many(config: RunConfig, read_paths: list[str],
             associations_path: str | None = None) -> CallResult:
    references = read_fasta(config.reference.fasta_path)
    all_fragments = []
    seen_fragment_ids: set[str] = set()
    for index, reads in enumerate(read_paths, start=1):
        _header_lengths, records = read_sam(Path(reads))
        source_id = sha256_file(reads)[:16]
        fragments = build_fragments(
            records,
            references,
            min_base_quality=config.parameters.min_base_quality,
            min_mapping_quality=config.parameters.min_mapping_quality,
            source_id=source_id,
            source_name=Path(reads).name,
            batch_id=f"batch-{index:04d}-{source_id}",
        )
        for fragment in fragments:
            if fragment.evidence_id not in seen_fragment_ids:
                seen_fragment_ids.add(fragment.evidence_id)
                all_fragments.append(fragment)
    relations = load_association_relations(
        associations_path or config.associations_path)
    constraints = _load_constraints(config, associations_path)
    return call_variants(all_fragments, references, config.parameters, relations,
                         constraints)
