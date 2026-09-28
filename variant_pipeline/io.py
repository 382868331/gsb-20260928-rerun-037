from __future__ import annotations

from pathlib import Path


def read_fasta(path: str | Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current = ""
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                sequences[current] = []
            elif current:
                sequences[current].append(line.upper())
            else:
                raise ValueError("FASTA sequence data appears before a sequence name")
    return {name: "".join(parts) for name, parts in sequences.items()}
