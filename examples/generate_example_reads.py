"""Generate examples/reads.sam with a known substitution and short indels.

Run from the repository root:
    .venv/bin/python examples/generate_example_reads.py

Research fixture only; reads are synthetic and deterministic.
"""

from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tests.helpers import (
    RNAME,
    REF,
    deletion_read,
    insertion_read,
    substitution_read,
)

OUT = Path(__file__).with_name("reads.sam")


def main() -> None:
    lines = [f"@HD\tVN:1.6\n", f"@SQ\tSN:{RNAME}\tLN:{len(REF)}\n"]
    # Known substitution g.15G>A: bidirectional support with opposition.
    lines += [substitution_read(f"snpF{i}", 6 + i, 20, 15, "A") for i in range(3)]
    lines += [substitution_read(f"snpR{i}", 8 + i, 20, 15, "A", reverse=True) for i in range(2)]
    lines += [substitution_read(f"ref{i}", 7 + i, 20, 15) for i in range(3)]
    # Known 2 bp insertion after anchor 20.
    lines += [insertion_read("insF1", 10, 20, 20, "AA")]
    lines += [insertion_read("insR1", 11, 20, 20, "AA", reverse=True)]
    lines += [insertion_read("insF2", 12, 20, 20, "AA")]
    # Known 2 bp deletion after anchor 20.
    lines += [deletion_read("delF1", 24, 18, 30, 2)]
    lines += [deletion_read("delR1", 25, 18, 30, 2, reverse=True)]
    lines += [deletion_read("delF2", 26, 18, 30, 2)]
    OUT.write_text("".join(lines), encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
