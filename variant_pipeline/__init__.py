"""Research-only candidate variant calling pipeline foundation.

This package summarizes candidate variants from aligned reads against a
manual reference sequence and records supporting/opposing evidence for each
candidate. It is intended for research workflow use only and does not emit
clinical or medical conclusions.
"""

__version__ = "0.4.0"

# Analysis rule version. Bumped when evidence semantics change; outputs
# produced under previous rule versions remain valid historical batches and
# are never rewritten in place.
ANALYSIS_VERSION = "orientation-forward-v3"

RESEARCH_DISCLAIMER = (
    "Research use only. Output is an unvalidated candidate summary for "
"workflow development; it is not a medical diagnosis, clinical report, or "
"substitute for professional genetic interpretation."
)
