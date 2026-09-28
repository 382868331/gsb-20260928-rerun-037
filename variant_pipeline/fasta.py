"""Minimal FASTA reader for manual reference sequences."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Reference:
    name: str
    sequence: str


def read_fasta(path: str | Path) -> dict[str, Reference]:
    """Return an ordered mapping of record name to Reference."""
    records: dict[str, Reference] = {}
    name: str | None = None
    chunks: list[str] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = Reference(name, "".join(chunks).upper())
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence line before any header")
                chunks.append(line)
    if name is not None:
        records[name] = Reference(name, "".join(chunks).upper())
    if not records:
        raise ValueError(f"no FASTA records found in {path}")
    return records
