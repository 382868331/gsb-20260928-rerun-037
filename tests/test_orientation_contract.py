"""Hand-written standard SAM records verify the orientation contract.

These records are authored independently of examples/generate_example_reads.py
so the checker can never inherit the generator's conventions. Every SEQ is
written in forward reference order; FLAG 0x10 only marks strand. Candidates,
fragment coverage and support counts must be identical regardless of the
read direction flag.
"""

from variant_pipeline.caller import call_variants
from variant_pipeline.config import Parameters
from variant_pipeline.fasta import read_fasta
from variant_pipeline.samio import build_fragments, read_sam

Q = "I"  # Phred 40

# Non-palindromic reference so reverse-complementing a matching SEQ is
# visibly wrong: positions 12-17 read AACGTT.
SIMPLE_REF = "ACGTACGTACGAACGTTCATGC"
SIMPLE_RNAME = "simple-1"

# 31 bp reference with non-tandem CA motif around positions 17-19.
INDEL_REF = "ACGTACGTACGTACGTGCAATGCATGCATGC"
INDEL_RNAME = "indel-1"

CALL_PARAMS = Parameters(
    min_base_quality=20,
    min_mapping_quality=1,
    min_coverage_for_calling=1,
    min_alt_count=1,
    min_alt_fraction=0.0,
    low_coverage_depth=1,
)


def _record(qname, flag, rname, pos, cigar, seq, qual=None):
    if qual is None:
        qual = Q * len(seq)
    return (
        f"{qname}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t"
        f"{seq}\t{qual}"
    )


def _run(tmp_path, ref_name, ref_seq, records, params=CALL_PARAMS):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">{ref_name}\n{ref_seq}\n", encoding="utf-8")
    body = [
        "@HD\tVN:1.6",
        f"@SQ\tSN:{ref_name}\tLN:{len(ref_seq)}",
        *records,
    ]
    sam = tmp_path / "reads.sam"
    sam.write_text("\n".join(body) + "\n", encoding="utf-8")
    refs = read_fasta(fasta)
    _, sam_records = read_sam(sam)
    fragments = build_fragments(
        sam_records, refs, params.min_base_quality, params.min_mapping_quality
    )
    return call_variants(fragments, refs, params), fragments


def test_reverse_matching_aacg_read_has_zero_variants(tmp_path):
    # Headline regression: AACG matching reference positions 12-15 (plus TT)
    # must not produce any candidate merely because FLAG 0x10 is set.
    assert SIMPLE_REF[11:17] == "AACGTT"
    records = [_record("revMatch", 16, SIMPLE_RNAME, 12, "6M", "AACGTT")]
    result, fragments = _run(tmp_path, SIMPLE_RNAME, SIMPLE_REF, records)
    assert result.candidates == []
    frag = fragments[0]
    assert frag.strand == "R"
    assert frag.alt == {}
    assert frag.ref_positions == {12, 13, 14, 15, 16, 17}
    assert frag.covered == [(12, 17)]
    depths = {c.pos: c.fragment_depth for c in result.coverage}
    assert all(depths[p] == 1 for p in range(12, 18))
    assert all(depths[p] == 0 for p in [*range(1, 12), *range(18, 23)])


def test_forward_and_reverse_matching_reads_give_identical_coverage(tmp_path):
    records = [
        _record("fwd", 0, SIMPLE_RNAME, 12, "6M", "AACGTT"),
        _record("rev", 16, SIMPLE_RNAME, 12, "6M", "AACGTT"),
    ]
    result, fragments = _run(tmp_path, SIMPLE_RNAME, SIMPLE_REF, records)
    assert result.candidates == []
    assert {f.strand for f in fragments} == {"F", "R"}
    assert all(not f.alt for f in fragments)
    depths = {c.pos: c.fragment_depth for c in result.coverage}
    assert all(depths[p] == 2 for p in range(12, 18))


def test_reverse_mismatch_is_one_snv_at_correct_position(tmp_path):
    assert SIMPLE_REF[11] == "A"  # position 12
    records = [
        _record("revAlt", 16, SIMPLE_RNAME, 10, "6M", "CGTACG"),
        _record("revRef", 16, SIMPLE_RNAME, 10, "6M", "CGAACG"),
        _record("fwdRef", 0, SIMPLE_RNAME, 10, "6M", "CGAACG"),
    ]
    result, _ = _run(tmp_path, SIMPLE_RNAME, SIMPLE_REF, records)
    assert len(result.candidates) == 1
    ev = result.candidates[0]
    assert ev.variant_id == "simple-1:g.12A>T"
    assert ev.alt_fragments == 1
    assert ev.alt_reverse == 1 and ev.alt_forward == 0
    assert ev.ref_fragments == 2
    assert ev.depth == 3


def test_reverse_insertion_is_placed_at_reference_anchor(tmp_path):
    assert INDEL_REF[11] == "T"  # anchor position 12
    seq = INDEL_REF[4:12] + "TT" + INDEL_REF[12:21]
    assert seq == "ACGTACGTTTACGTGCAAT"
    records = [_record("revIns", 16, INDEL_RNAME, 5, "8M2I9M", seq)]
    result, fragments = _run(tmp_path, INDEL_RNAME, INDEL_REF, records)
    assert len(result.candidates) == 1
    ev = result.candidates[0]
    assert ev.kind == "INS"
    assert ev.variant_id == "indel-1:g.12_13insTT"
    assert ev.pos == 12
    assert ev.alt_fragments == 1 and ev.alt_reverse == 1
    assert fragments[0].covered == [(5, 21)]


def test_reverse_deletion_is_placed_at_reference_anchor(tmp_path):
    assert INDEL_REF[16] == "G"  # anchor position 17
    assert INDEL_REF[17:19] == "CA"  # deleted positions 18-19
    seq = INDEL_REF[4:17] + INDEL_REF[19:31]
    assert len(seq) == 25
    records = [_record("revDel", 16, INDEL_RNAME, 5, "13M2D12M", seq)]
    result, fragments = _run(tmp_path, INDEL_RNAME, INDEL_REF, records)
    assert len(result.candidates) == 1
    ev = result.candidates[0]
    assert ev.kind == "DEL"
    assert ev.variant_id == "indel-1:g.17_19delCA"
    assert ev.alt_fragments == 1 and ev.alt_reverse == 1
    assert fragments[0].covered == [(5, 31)]
    depths = {c.pos: c.fragment_depth for c in result.coverage}
    assert all(depths[p] == 1 for p in range(5, 32))
    assert all(depths[p] == 0 for p in [*range(1, 5), *range(32, 32)])


def test_soft_clips_on_either_end_are_ignored_in_both_directions(tmp_path):
    assert SIMPLE_REF[15] == "T"  # position 16
    records = [
        # leading soft clip on a reverse read; the Ns consume no reference
        _record("leadS", 16, SIMPLE_RNAME, 12, "4S6M", "NNNNAACGCT"),
        # trailing soft clip on a forward read
        _record("trailS", 0, SIMPLE_RNAME, 12, "6M4S", "AACGTTNNNN"),
        _record("ref1", 0, SIMPLE_RNAME, 12, "6M", "AACGTT"),
    ]
    result, fragments = _run(tmp_path, SIMPLE_RNAME, SIMPLE_REF, records)
    assert len(result.candidates) == 1
    ev = result.candidates[0]
    assert ev.variant_id == "simple-1:g.16T>C"
    assert ev.alt_fragments == 1 and ev.alt_reverse == 1
    assert ev.ref_fragments == 2
    for frag in fragments:
        assert frag.covered == [(12, 17)]
    depths = {c.pos: c.fragment_depth for c in result.coverage}
    assert all(depths[p] == 3 for p in range(12, 18))
    assert all(depths[p] == 0 for p in range(8, 12))


def test_overlapping_reverse_mate_votes_once(tmp_path):
    assert SIMPLE_REF[5:15] == "CGTACGAACG"
    assert SIMPLE_REF[9:19] == "CGAACGTTCA"
    records = [
        _record("pair1", 0x1 | 0x40, SIMPLE_RNAME, 6, "10M", "CGTACGTACG"),
        _record("pair1", 0x1 | 0x80 | 0x10, SIMPLE_RNAME, 10, "10M", "CGTACGTTCA"),
        _record("s2", 0, SIMPLE_RNAME, 6, "10M", "CGTACGTACG"),
        _record("r1", 0, SIMPLE_RNAME, 6, "10M", "CGTACGAACG"),
        _record("r2", 0, SIMPLE_RNAME, 6, "10M", "CGTACGAACG"),
    ]
    result, _ = _run(tmp_path, SIMPLE_RNAME, SIMPLE_REF, records)
    assert len(result.candidates) == 1
    ev = result.candidates[0]
    assert ev.variant_id == "simple-1:g.12A>T"
    assert ev.alt_fragments == 2  # merged pair counts once + one single
    assert ev.depth == 4  # pair1 is one fragment, never one vote per mate


def test_conflicting_overlapping_mates_keep_abstention_reason(tmp_path):
    records = [
        _record("pair1", 0x1 | 0x40, SIMPLE_RNAME, 6, "10M", "CGTACGTACG"),
        _record("pair1", 0x1 | 0x80 | 0x10, SIMPLE_RNAME, 10, "10M", "CGAACGTTCA"),
        _record("a1", 0, SIMPLE_RNAME, 6, "10M", "CGTACGTACG"),
        _record("r1", 0, SIMPLE_RNAME, 10, "10M", "CGAACGTTCA"),
    ]
    result, fragments = _run(tmp_path, SIMPLE_RNAME, SIMPLE_REF, records)
    pair = next(f for f in fragments if f.qname == "pair1")
    assert pair.conflict_positions == {12}
    assert pair.alt == {}
    ev = result.candidates[0]
    assert ev.variant_id == "simple-1:g.12A>T"
    assert ev.alt_fragments == 1
    assert ev.ref_fragments == 1
    assert ev.conflicting_overlap_fragments == 1
    assert ev.depth == 3
    assert "MATE_OVERLAP_CONFLICT" in ev.flags
