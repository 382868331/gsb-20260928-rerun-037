"""Helpers to synthesize SAM reads against a simple ACGT-repeat reference.

All reads are authored directly on the forward reference axis: SEQ follows
the reference 5'->3' order and CIGAR is expressed in reference coordinates,
regardless of FLAG 0x10. The strand flag only labels read direction; the
record stores SEQ/QUAL as written, exactly like standard records produced
under the pipeline's documented orientation convention.
"""

REF = "ACGT" * 13  # 52 bp
RNAME = "synth-1"
Q40 = "I"  # Phred 40
Q10 = "+"  # Phred 10


def ref_slice(start: int, length: int) -> str:
    return REF[start - 1:start - 1 + length]


def substitution_read(
    qname: str,
    start: int,
    length: int,
    pos: int,
    base: str | None = None,
    reverse: bool = False,
    mapq: int = 60,
    qual: str = Q40,
    paired_flags: int | None = None,
):
    seq = list(ref_slice(start, length))
    if base is not None:
        seq[pos - start] = base
    seq = "".join(seq)
    flag = paired_flags if paired_flags is not None else (16 if reverse else 0)
    return _record(qname, flag, start, f"{length}M", seq, qual * length, mapq)


def insertion_read(
    qname: str,
    start: int,
    length: int,
    anchor_pos: int,
    inserted: str,
    reverse: bool = False,
):
    before = REF[start - 1:anchor_pos]
    after = REF[anchor_pos:start - 1 + length]
    seq = before + inserted + after
    x = anchor_pos - start + 1
    cigar = f"{x}M{len(inserted)}I{length - x}M"
    flag = 16 if reverse else 0
    return _record(qname, flag, start, cigar, seq, Q40 * len(seq))


def deletion_read(
    qname: str,
    start: int,
    read_match: int,
    anchor_pos: int,
    deleted_len: int,
    reverse: bool = False,
):
    before = REF[start - 1:anchor_pos]
    after = REF[anchor_pos + deleted_len:anchor_pos + deleted_len + read_match - len(before)]
    seq = before + after
    x = anchor_pos - start + 1
    y = read_match - x
    cigar = f"{x}M{deleted_len}D{y}M"
    flag = 16 if reverse else 0
    return _record(qname, flag, start, cigar, seq, Q40 * len(seq))


def _record(qname, flag, start, cigar, seq, qual, mapq=60):
    return (
        f"{qname}\t{flag}\t{RNAME}\t{start}\t{mapq}\t{cigar}\t*\t0\t0\t"
        f"{seq}\t{qual}\n"
    )


def write_case(tmp_path, lines, fasta_name="ref.fasta", sam_name="reads.sam"):
    fasta = tmp_path / fasta_name
    fasta.write_text(f">{RNAME}\n{REF}\n", encoding="utf-8")
    sam = tmp_path / sam_name
    header = f"@HD\tVN:1.6\n@SQ\tSN:{RNAME}\tLN:{len(REF)}\n"
    sam.write_text(header + "".join(lines), encoding="utf-8")
    return fasta, sam
