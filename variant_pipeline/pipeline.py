"""End-to-end orchestration: config -> references -> fragments -> calls."""

from __future__ import annotations

import hashlib
from pathlib import Path

from .caller import CallResult, call_variants
from .config import RunConfig
from .fasta import read_fasta
from .samio import build_fragments, read_sam
from .store import StoreError


def _batch_specs(config: RunConfig, reads_path: str | None) -> list[dict]:
    if reads_path:
        return [{"path": reads_path, "label": Path(reads_path).name}]
    if config.read_batches:
        return [dict(batch) for batch in config.read_batches]
    if config.reads_path:
        return [{"path": config.reads_path,
                 "label": Path(config.reads_path).name}]
    raise ValueError("no reads path: set reads_path in config or pass --reads")


def _batch_content_id(path: str) -> str:
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]
    return f"batch-{digest}"


def run(config: RunConfig, reads_path: str | None = None) -> CallResult:
    references = read_fasta(config.reference.fasta_path)
    if config.reference.name not in references:
        raise ValueError(
            f"reference name {config.reference.name!r} not found in FASTA"
        )
    fragments = []
    imported_batches = []
    skipped_batches = []
    seen_batch_ids = set()
    seen_content_ids = set()
    for spec in _batch_specs(config, reads_path):
        batch_path = spec["path"]
        label = spec.get("label") or Path(batch_path).name
        batch = label
        content_id = _batch_content_id(batch_path)
        if batch in seen_batch_ids and content_id in seen_content_ids:
            skipped_batches.append({
                "path": batch_path,
                "batch_id": batch,
                "content_batch_id": content_id,
            })
            continue
        if batch in seen_batch_ids:
            raise StoreError(
                f"batch label {batch!r} changed content; use a new stable label")
        _header_lengths, records = read_sam(Path(batch_path))
        if content_id in seen_content_ids:
            fragments.extend(build_fragments(
                records,
                references,
                min_base_quality=config.parameters.min_base_quality,
                min_mapping_quality=config.parameters.min_mapping_quality,
                batch_id=batch,
                batch_label=label,
            ))
            imported_batches.append({
                "path": batch_path,
                "label": label,
                "batch_id": batch,
                "content_batch_id": content_id,
                "record_count": len(records),
                "distinct_copy_of_content_batch": True,
            })
            seen_batch_ids.add(batch)
            continue
        seen_batch_ids.add(batch)
        seen_content_ids.add(content_id)
        fragments.extend(build_fragments(
            records,
            references,
            min_base_quality=config.parameters.min_base_quality,
            min_mapping_quality=config.parameters.min_mapping_quality,
            batch_id=batch,
            batch_label=label,
        ))
        imported_batches.append({
            "path": batch_path,
            "label": label,
            "batch_id": batch,
            "content_batch_id": content_id,
            "record_count": len(records),
        })
    result = call_variants(fragments, references, config.parameters)
    result.evidence_ledger["source_batches"] = imported_batches
    result.evidence_ledger["duplicate_imports_skipped"] = skipped_batches
    return result
