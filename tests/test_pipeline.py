"""Validate known substitutions and short indels against synthetic data."""

from variant_pipeline.caller import call_variants
from variant_pipeline.config import Parameters
from variant_pipeline.fasta import read_fasta
from variant_pipeline.normalize import normalize
from variant_pipeline.samio import build_fragments, read_sam
from tests.helpers import (
    REF,
    RNAME,
    deletion_read,
    insertion_read,
    substitution_read,
    write_case,
)


PARAMS = Parameters()


def _call(tmp_path, lines, params=PARAMS):
    fasta, sam = write_case(tmp_path, lines)
    refs = read_fasta(fasta)
    _, records = read_sam(sam)
    fragments = build_fragments(records, refs, params.min_base_quality, params.min_mapping_quality)
    return call_variants(fragments, refs, params), refs


def test_known_substitution_detected_with_bidirectional_support(tmp_path):
    pos = 15
    alt = "A"  # REF pos 15 = G (pattern ACGT, idx14 -> G)
    assert REF[pos - 1] == "G"
    lines = []
    # 4 alt (2 forward, 2 reverse) + 3 ref spanning reads
    lines += [substitution_read(f"a{i}", 6 + i, 20, pos, alt, reverse=False) for i in range(2)]
    lines += [substitution_read(f"r{i}", 8 + i, 20, pos, alt, reverse=True) for i in range(2)]
    lines += [substitution_read(f"ref{i}", 7 + i, 20, pos) for i in range(3)]
    result, _ = _call(tmp_path, lines)
    assert len(result.candidates) == 1
    ev = result.candidates[0]
    assert ev.variant_id == f"{RNAME}:g.15G>A"
    assert ev.kind == "SNV"
    assert ev.status == "PASS"
    assert ev.alt_fragments == 4
    assert ev.ref_fragments == 3
    assert ev.alt_forward == 2 and ev.alt_reverse == 2
    assert ev.strand_bias_p > 0.05
    assert ev.depth == 7


def test_strand_bias_is_flagged_not_counted_away(tmp_path):
    pos = 15
    assert REF[pos - 1] == "G"
    lines = []
    lines += [substitution_read(f"f{i}", 6 + i, 20, pos, "A", reverse=False) for i in range(5)]
    lines += [substitution_read(f"r{i}", 7 + i, 20, pos, reverse=True) for i in range(5)]
    result, _ = _call(tmp_path, lines)
    ev = result.candidates[0]
    assert ev.alt_fragments == 5
    assert ev.alt_reverse == 0
    assert "STRAND_BIAS" in ev.flags
    assert "ONE_STRAND_ONLY" in ev.flags
    assert ev.status == "SUSPECT"
    assert ev.strand_bias_p < 0.01


def test_low_quality_alt_lowers_confidence_and_status(tmp_path):
    pos = 15
    assert REF[pos - 1] == "G"
    lines = [
        substitution_read("hq1", 6, 20, pos, "A"),
        substitution_read("hq2", 7, 20, pos, "A", reverse=True),
        substitution_read("lq1", 8, 20, pos, "A", qual="+"),
        substitution_read("lq2", 9, 20, pos, "A", qual="+", reverse=True),
        substitution_read("r1", 10, 20, pos),
        substitution_read("r2", 11, 20, pos),
        substitution_read("r3", 12, 20, pos),
    ]
    result, _ = _call(tmp_path, lines)
    ev = result.candidates[0]
    assert ev.alt_fragments == 4
    assert ev.low_quality_alt_fragments == 2
    assert "LOW_QUALITY_EVIDENCE" in ev.flags
    assert ev.status == "SUSPECT"
    assert ev.components["quality_consistency"] == 0.5


def test_insertion_candidate(tmp_path):
    anchor = 20
    lines = [
        insertion_read("i1", 10, 20, anchor, "AA"),
        insertion_read("i2", 11, 20, anchor, "AA", reverse=True),
        insertion_read("i3", 12, 20, anchor, "AA"),
        substitution_read("r1", 10, 22, anchor),
        substitution_read("r2", 11, 22, anchor, reverse=True),
        substitution_read("r3", 12, 22, anchor),
    ]
    result, _ = _call(tmp_path, lines)
    ev = result.candidates[0]
    assert ev.kind == "INS"
    assert ev.variant_id == f"{RNAME}:g.20_21insAA"
    assert ev.alt_fragments == 3
    assert ev.status == "PASS"


def test_deletion_candidate(tmp_path):
    anchor = 20
    lines = [
        deletion_read("d1", 10, 20, anchor, 2),
        deletion_read("d2", 11, 20, anchor, 2, reverse=True),
        deletion_read("d3", 12, 20, anchor, 2),
        substitution_read("r1", 10, 22, anchor),
        substitution_read("r2", 11, 22, anchor),
        substitution_read("r3", 12, 22, anchor, reverse=True),
    ]
    result, refs = _call(tmp_path, lines)
    ev = result.candidates[0]
    assert ev.kind == "DEL"
    assert ev.variant_id == f"{RNAME}:g.20_22delAC"
    assert ev.alt_fragments == 3
    assert ev.status == "PASS"


def test_left_normalization_collapses_repeat_spellings(tmp_path):
    fasta = tmp_path / "r.fasta"
    fasta.write_text(f">c\nAAAAAC\n")
    from variant_pipeline.fasta import Reference

    ref = Reference("c", "AAAAAC")
    # Insertion of A anchored on the 5th A vs anchored on the 4th A must
    # normalize to the same identity.
    v1 = normalize("c", 5, "A", "AAA", ref)   # ..AA[A-ins]A...
    v2 = normalize("c", 4, "A", "AAA", ref)
    assert v1.canonical_key == v2.canonical_key
    assert v1.vid == "c:g.1_2insAA"


def test_overlapping_mates_counted_once(tmp_path):
    pos = 15
    assert REF[pos - 1] == "G"
    lines = [
        # same fragment, both mates alt, identical call -> one vote
        substitution_read(
            "pair1", 6, 20, pos, "A", paired_flags=0x1 | 0x40
        ),
        substitution_read(
            "pair1", 10, 20, pos, "A", paired_flags=0x1 | 0x80 | 0x10
        ),
        substitution_read("s2", 7, 20, pos, "A"),
        substitution_read("r1", 8, 20, pos),
        substitution_read("r2", 9, 20, pos),
        substitution_read("r3", 10, 20, pos),
    ]
    result, _ = _call(tmp_path, lines)
    ev = result.candidates[0]
    assert ev.alt_fragments == 2  # pair1 merged + s2
    assert ev.depth == 5


def test_conflicting_overlapping_mates_abstain(tmp_path):
    pos = 15
    assert REF[pos - 1] == "G"
    lines = [
        substitution_read("p1", 6, 20, pos, "A", paired_flags=0x1 | 0x40),
        substitution_read("p1", 10, 20, pos, "T", paired_flags=0x1 | 0x80 | 0x10),
        substitution_read("a1", 7, 20, pos, "A"),
        substitution_read("a2", 8, 20, pos, "A", reverse=True),
        substitution_read("r1", 9, 20, pos),
        substitution_read("r2", 10, 20, pos),
        substitution_read("r3", 11, 20, pos),
    ]
    result, _ = _call(tmp_path, lines)
    ev = result.candidates[0]
    assert ev.alt_fragments == 2
    assert ev.conflicting_overlap_fragments == 1
    assert "MATE_OVERLAP_CONFLICT" in ev.flags


def test_low_coverage_reported_separately(tmp_path):
    pos = 30
    lines = [
        substitution_read("a1", pos - 5, 12, pos, "T"),
        substitution_read("a2", pos - 5, 12, pos, "T", reverse=True),
        substitution_read("r1", pos - 5, 12, pos),
    ]
    params = Parameters(low_coverage_depth=5, min_coverage_for_calling=5)
    result, _ = _call(tmp_path, lines, params)
    ev = result.candidates[0]
    assert ev.depth == 3
    assert ev.status == "LOW_COVERAGE"
    assert "INSUFFICIENT_COVERAGE" in ev.flags
    # flanks outside the small reads must be listed as low coverage
    assert any(region[0] == RNAME and region[1] == 1 for region in result.low_coverage)
    assert any(region[2] == len(REF) for region in result.low_coverage)


def test_reverse_read_sequence_is_oriented(tmp_path):
    from variant_pipeline.samio import reverse_complement

    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AACG") == "CGTT"
