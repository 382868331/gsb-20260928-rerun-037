"""Research-only candidate variant calling pipeline foundation.

This package summarizes candidate variants from aligned reads against a
manual reference sequence and records supporting/opposing evidence for each
candidate. It is intended for research workflow use only and does not emit
clinical or medical conclusions.
"""

__version__ = "0.2.0"

# Version of the strand/orientation and counting rules used during analysis.
# Bumped for the fix that keeps SEQ in forward-reference orientation for
# reverse-strand SAM records (no second reverse complement) and treats
# fragments as the single unit of evidence. Old output batches are never
# rewritten; new runs carry this version in their payload.
RULES_VERSION = "orientation-rules-v2"

RESEARCH_DISCLAIMER = (
    "Research use only. Output is an unvalidated candidate summary for "
"workflow development; it is not a medical diagnosis, clinical report, or "
"substitute for professional genetic interpretation."
)
