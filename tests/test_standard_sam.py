"""Validate the parser against independently hand-written standard SAM.

These records are written as literal SAM lines (not via the fixture
generator) so they pin down the SAM storage convention independently:
SEQ/QUAL are forward-reference oriented for every alignment, FLAG 0x10 only
marks a reverse-strand alignment, and CIGAR walks increasing reference
coordinates.
"""

from variant_pipeline.caller import call_variants
from variant_pipeline.config import Parameters
from variant_pipeline.fasta import read_fasta
from variant_pipeline.samio import build_fragments, read_sam

Q = "I" * 60  # Phred 40


def _run(tmp_path, fasta_body, sam_lines, params=Parameters()):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">chr\n{fasta_body}\n", encoding="utf-8")
    sam = tmp_path / "reads.sam"
    header = f"@HD\tVN:1.6\n@SQ\tSN:chr\tLN:{len(fasta_body)}\n"
    sam.write_text(header + "".join(sam_lines), encoding="utf-8")
    refs = read_fasta(fasta)
    _, records = read_sam(sam)
    fragments = build_fragments(
        records, refs, params.min_base_quality, params.min_mapping_quality
    )
    return call_variants(fragments, refs, params), fragments


def _line(qname, flag, pos, cigar, seq, qual=None):
    qual = qual or Q[: len(seq)]
    return f"{qname}\t{flag}\tchr\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n"


def test_reverse_exact_match_yields_zero_variants(tmp_path):
    # AACG is present in the reference (positions 5-8). A reverse-strand
    # read that is an exact match stores the SAME forward-oriented bases in
    # SEQ; re-reverse-complementing it must not create four false SNVs.
    ref = "ACGTAACGTT"  # len 10, AACG at 5-8
    assert ref[4:8] == "AACG"
    lines = [
        _line("rfwd", 0, 5, "4M", "AACG"),
        _line("rrev", 16, 5, "4M", "AACG"),
    ]
    result, fragments = _run(tmp_path, ref, lines)
    assert result.candidates == []
    depth = {p.pos: p.fragment_depth for p in result.coverage}
    for pos in range(5, 9):
        assert depth[pos] == 2
    assert depth[1] == 0 and depth[10] == 0
    assert {f.strand for f in fragments} == {"F", "R"}


def test_standard_sam_reverse_snv_observed_on_correct_strand(tmp_path):
    # 18 bp reference; SNV at pos 9: ref C -> alt T.
    ref = "AAAAAAAACTTTTTTTTT"
    assert ref[8] == "C"
    forward_alt = ref[:8] + "T" + ref[9:]
    reverse_alt = forward_alt  # standard SAM: same forward-oriented SEQ
    lines = [
        _line("f1", 0, 1, "18M", forward_alt),
        _line("r1", 16, 1, "18M", reverse_alt),
        _line("ref1", 0, 1, "18M", ref),
        _line("ref2", 16, 1, "18M", ref),
        _line("ref3", 0, 1, "18M", ref),
    ]
    params = Parameters(min_coverage_for_calling=5, low_coverage_depth=5)
    result, _ = _run(tmp_path, ref, lines, params)
    assert len(result.candidates) == 1
    ev = result.candidates[0]
    assert ev.variant_id == "chr:g.9C>T"
    assert ev.alt_fragments == 2
    assert ev.alt_forward == 1 and ev.alt_reverse == 1
    assert ev.ref_fragments == 3
    assert ev.depth == 5
    assert ev.status == "PASS"


def test_reverse_soft_clip_and_deletion_use_reference_coordinates(tmp_path):
    # Reverse read covering pos 6..19: 2 bp soft clip at the stored-sequence
    # start (reference-left end) and a 2 bp deletion of reference pos 11-12.
    ref = "AAAAAGGGGGTTCCCCCCCC"  # positions 1..20
    assert ref[5:10] == "GGGGG" and ref[10:12] == "TT" and ref[12:20] == "CCCCCCCC"
    seq = "NN" + "GGGGG" + "CCCCCCC"
    lines = [_line("rdel", 16, 6, "2S5M2D7M", seq, "#" + Q[1:14])]
    result, fragments = _run(tmp_path, ref, lines)
    frag = fragments[0]
    assert frag.covered == [(6, 19)]
    assert len(result.candidates) == 1
    ev = result.candidates[0]
    assert ev.kind == "DEL"
    assert ev.variant_id in ("chr:g.10_12delGTT", "chr:g.10_12delTT")
    depth = {p.pos: p.fragment_depth for p in result.coverage}
    assert depth[5] == 0
    assert depth[8] == 1
    # the fragment spans the gap, but the gap bases carry no base observation
    assert depth[11] == 1 and depth[12] == 1
    assert not ({11, 12} & fragments[0].ref_positions)
    assert depth[15] == 1
    assert depth[20] == 0


def test_paired_overlap_reverse_reads_counted_once(tmp_path):
    ref = "AAAAAAAACTTTTTTTTT"  # SNV ref C -> T at pos 9
    alt = ref[:8] + "T" + ref[9:]
    lines = [
        # one fragment, both mates show the same ALT -> one vote
        _line("p1", 0x1 | 0x40, 1, "18M", alt),
        _line("p1", 0x1 | 0x80 | 0x10, 1, "18M", alt),
        _line("s1", 0, 1, "18M", alt),
        _line("r1", 0, 1, "18M", ref),
        _line("r2", 16, 1, "18M", ref),
        _line("r3", 0, 1, "18M", ref),
    ]
    params = Parameters(min_coverage_for_calling=5, low_coverage_depth=5)
    result, fragments = _run(tmp_path, ref, lines, params)
    assert len(fragments) == 5
    ev = result.candidates[0]
    assert ev.variant_id == "chr:g.9C>T"
    assert ev.alt_fragments == 2
    assert ev.alt_forward == 2 and ev.alt_reverse == 0
    assert ev.depth == 5


def test_paired_overlap_conflict_abstains_on_reverse_mate(tmp_path):
    ref = "AAAAAAAACTTTTTTTTT"  # pos 9: read1 ALT T vs read2 REF C
    alt = ref[:8] + "T" + ref[9:]
    lines = [
        _line("p1", 0x1 | 0x40, 1, "18M", alt),
        _line("p1", 0x1 | 0x80 | 0x10, 1, "18M", ref),
        _line("a1", 0, 1, "18M", alt),
        _line("a2", 16, 1, "18M", alt),
        _line("r1", 0, 1, "18M", ref),
        _line("r2", 0, 1, "18M", ref),
        _line("r3", 0, 1, "18M", ref),
    ]
    params = Parameters(min_coverage_for_calling=5, low_coverage_depth=5)
    result, _ = _run(tmp_path, ref, lines, params)
    ev = result.candidates[0]
    assert ev.alt_fragments == 2
    assert ev.ref_fragments == 3
    assert ev.conflicting_overlap_fragments == 1
    assert ev.depth == 6
    assert "MATE_OVERLAP_CONFLICT" in ev.flags
