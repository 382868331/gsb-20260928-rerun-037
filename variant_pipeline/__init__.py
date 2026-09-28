"""Research-only candidate variant calling pipeline foundation.

This package summarizes candidate variants from aligned reads against a
manual reference sequence and records supporting/opposing evidence for each
candidate. It is intended for research workflow use only and does not emit
clinical or medical conclusions.
"""

__version__ = "0.1.0"

RESEARCH_DISCLAIMER = (
    "Research use only. Output is an unvalidated candidate summary for "
"workflow development; it is not a medical diagnosis, clinical report, or "
"substitute for professional genetic interpretation."
)
