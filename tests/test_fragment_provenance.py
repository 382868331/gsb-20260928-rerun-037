"""Shared ambiguity, source provenance, and tentative interpretation tests."""

import json

from variant_pipeline.cli import main
from variant_pipeline.diff import compare_versions
from variant_pipeline.selection import ambiguity_groups, load_payload, project
from variant_pipeline.store import RunStore
from tests.helpers import Q40, REF, RNAME, write_case
from variant_pipeline.caller import call_variants
from variant_pipeline.config import Parameters
from variant_pipeline.fasta import Reference
from variant_pipeline.normalize import Variant
from variant_pipeline.samio import Fragment, Observation


def _record(qname, start, seq, cigar):
    return (
        f"{qname}\t0\t{RNAME}\t{start}\t60\t{cigar}\t*\t0\t0\t"
        f"{seq}\t{Q40 * len(seq)}\n"
    )


def _sub(qname, start, pos, alt=None):
    seq = list(REF[start - 1:start - 1 + 20])
    if alt is not None:
        seq[pos - start] = alt
    return _record(qname, start, "".join(seq), "20M")


def _ins(qname, start, anchor):
    before = REF[start - 1:anchor]
    after = REF[anchor:start - 1 + 18]
    seq = before + "AA" + after
    x = anchor - start + 1
    return _record(qname, start, seq, f"{x}M2I{18 - x}M")


def _ambiguous_anchor(qname, start, anchor):
    before = REF[start - 1:anchor]
    after = REF[anchor + 1:start - 1 + 19]
    seq = before + "T" + "A" + after
    x = anchor - start + 1
    return _record(qname, start, seq, f"{x}M1I{19 - x}M")


def _variant(chrom=RNAME, pos=14, ref="C", alt="CT", kind="INS"):
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, kind=kind)


def _ambiguous_fragments():
    ins = _variant(alt="CT", kind="INS")
    snv = _variant(alt="T", kind="SNV")
    shared = Fragment(qname="shared", strand="F", uid="source-1\0shared",
                      source_batch_id="source-1", source_label="source-1")
    shared.covered = [(10, 29)]
    shared.alt = {
        ins.canonical_key: Observation(ins, (40,), "F", False,
                                       shared.uid, "source-1", "source-1"),
        snv.canonical_key: Observation(snv, (40,), "F", False,
                                       shared.uid, "source-1", "source-1"),
    }
    persistent = Fragment(qname="persist", strand="F",
                          uid="source-1\0persist", source_batch_id="source-1")
    persistent.covered = [(10, 29)]
    persistent.alt = {
        key: Observation(obs.variant, obs.base_quals, obs.strand, False,
                         persistent.uid, "source-1", "source-1")
        for key, obs in shared.alt.items()
    }
    only_ins = Fragment(qname="ins-only", strand="F",
                        uid="source-1\0ins-only", source_batch_id="source-1")
    only_ins.covered = [(10, 29)]
    only_ins.alt = {ins.canonical_key: Observation(
        ins, (40,), "F", False, only_ins.uid, "source-1", "source-1")}
    only_snv = Fragment(qname="snv-only", strand="F",
                        uid="source-1\0snv-only", source_batch_id="source-1")
    only_snv.covered = [(10, 29)]
    only_snv.alt = {snv.canonical_key: Observation(
        snv, (40,), "F", False, only_snv.uid, "source-1", "source-1")}
    refs = []
    for name in ("r1", "r2", "r3"):
        fragment = Fragment(qname=name, strand="F",
                            uid=f"source-1\0{name}", source_batch_id="source-1")
        fragment.covered = [(10, 29)]
        fragment.ref_positions = {14}
        refs.append(fragment)
    return [shared, persistent, only_ins, only_snv, *refs]


def _config(case_dir, fasta, batches, labels=None):
    config = {
        "sample": {"sample_id": "S1"},
        "reference": {
            "fasta_path": str(fasta),
            "name": RNAME,
            "assembly": "synth-v1",
        },
        "read_batches": [
            {"path": str(path), "label": label}
            for path, label in zip(batches, labels or [f"batch-{i + 1}"
                                                       for i in range(len(batches))])
        ],
    }
    config_path = case_dir / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    return config_path


def _write_batch(case_dir, name, lines):
    path = case_dir / name
    header = f"@HD\tVN:1.6\n@SQ\tSN:{RNAME}\tLN:{len(REF)}\n"
    path.write_text(header + "".join(lines), encoding="utf-8")
    return path


def test_shared_fragment_is_not_counted_as_new_independent_support(tmp_path):
    fragments = _ambiguous_fragments()
    result = call_variants(fragments, {RNAME: Reference(RNAME, REF)}, Parameters())
    payload = {"candidates": [
        {"variant_id": ev.variant_id, "kind": ev.kind,
         "supporting_evidence": ev.supporting_summary()}
        for ev in result.candidates
    ], "evidence_ledger": result.evidence_ledger}
    sub = next(c for c in payload["candidates"] if c["kind"] == "SNV")
    ins = next(c for c in payload["candidates"] if c["kind"] == "INS")
    assert sub["supporting_evidence"]["alt_fragments"] == 3
    assert sub["supporting_evidence"]["independent_alt_fragments"] == 1
    assert sub["supporting_evidence"]["shared_ambiguous_alt_fragments"] == 2
    assert ins["supporting_evidence"]["independent_alt_fragments"] == 1
    assert ins["supporting_evidence"]["shared_ambiguous_alt_fragments"] == 2
    groups = ambiguity_groups(payload)
    assert len(groups) == 2
    assert any(group["fragment_uid"].endswith("shared") for group in groups)


def test_same_qname_from_different_source_is_distinct_evidence(tmp_path):
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    fasta, _ = write_case(case_dir, [])
    b1 = _write_batch(case_dir, "b1.sam", [
        _sub("same-name", 10, 15, "A"),
        _sub("r1", 10, 15, "A"),
        _sub("r2", 10, 15, "A"),
        _sub("r3", 10, 15, "A"),
        _sub("ref1", 10, 15),
        _sub("ref2", 10, 15),
        _sub("ref3", 10, 15),
    ])
    b2 = _write_batch(case_dir, "b2.sam", [
        _ins("same-name", 10, 14),
        _ins("i1", 10, 14),
        _ins("i2", 10, 14),
        _ins("i3", 10, 14),
    ])
    out = tmp_path / "out"
    assert main(["call", "--config", str(_config(
        case_dir, fasta, [b1, b2], ["source-a", "source-b"])),
                 "--outdir", str(out)]) == 0
    payload = load_payload(RunStore(out), "v0001")
    same_names = [
        record for record in payload["evidence_ledger"]["fragment_records"]
        if record["qname"] == "same-name"
    ]
    assert len(same_names) == 2
    assert {record["source_label"] for record in same_names} == {
        "source-a", "source-b"
    }
    assert all(not record["ambiguous_between"] for record in same_names)


def test_duplicate_batch_import_is_skipped_and_does_not_inflate_support(tmp_path):
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    fasta, _ = write_case(case_dir, [])
    batch = _write_batch(case_dir, "b1.sam", [
        _sub("a", 10, 15, "A"),
        _sub("b", 10, 15, "A"),
        _sub("r1", 10, 15),
        _sub("r2", 10, 15),
        _sub("r3", 10, 15),
    ])
    out = tmp_path / "out"
    assert main(["call", "--config", str(_config(
        case_dir, fasta, [batch, batch], ["dup", "dup"])),
                 "--outdir", str(out)]) == 0
    payload = load_payload(RunStore(out), "v0001")
    assert payload["candidates"][0]["supporting_evidence"]["alt_fragments"] == 2
    assert len(payload["evidence_ledger"]["fragment_records"]) == 5
    assert len(payload["evidence_ledger"]["duplicate_imports_skipped"]) == 1


def test_adopt_whatif_withdraw_leaves_run_unchanged(tmp_path):
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    fasta, _ = write_case(case_dir, [])
    direct = call_variants(_ambiguous_fragments(),
                           {RNAME: Reference(RNAME, REF)}, Parameters())
    direct_payload = {
        "candidates": [
            {"variant_id": ev.variant_id, "kind": ev.kind,
             "supporting_evidence": ev.supporting_summary()}
            for ev in direct.candidates
        ],
        "evidence_ledger": direct.evidence_ledger,
    }
    direct_sub = next(c["variant_id"] for c in direct_payload["candidates"]
                      if c["kind"] == "SNV")
    direct_ins = next(c["variant_id"] for c in direct_payload["candidates"]
                      if c["kind"] == "INS")
    direct_view = project(direct_payload, {direct_sub})
    direct_by_id = {c["variant_id"]: c for c in direct_view["candidates"]}
    assert direct_by_id[direct_sub]["projected_alt_fragments"] == 3
    assert direct_by_id[direct_ins]["projected_alt_fragments"] == 1

    sam = _write_batch(case_dir, "reads.sam", [
        _sub("sub-only", 10, 15, "A"),
        _sub("ref1", 10, 15),
    ])
    out = tmp_path / "out"
    main(["call", "--config", str(_config(case_dir, fasta, [sam])),
          "--outdir", str(out)])
    before = sorted(p.read_bytes() for p in (out / "runs" / "v0001").iterdir())
    payload = load_payload(RunStore(out), "v0001")
    sub = next(c["variant_id"] for c in payload["candidates"]
               if c["kind"] == "SNV")
    assert main(["what-if", "--outdir", str(out), "--version", "v0001",
                 "--variant", sub,
                 "--json-out", str(tmp_path / "proj.json")]) == 0
    projection_path = tmp_path / "proj.json"
    assert main(["what-if", "--outdir", str(out), "--version", "v0001",
                 "--json-out", str(projection_path)]) == 0
    projection = json.loads(projection_path.read_text(encoding="utf-8"))
    by_id = {c["variant_id"]: c for c in projection["candidates"]}
    assert by_id[sub]["projected_alt_fragments"] == 1
    after = sorted(p.read_bytes() for p in (out / "runs" / "v0001").iterdir())
    assert before == after
    assert project(payload, set())["selected_variant_ids"] == []


def test_compare_reports_resolved_and_persistent_ambiguity(tmp_path):
    from variant_pipeline.diff import _ambiguity_reconciliation

    result = call_variants(_ambiguous_fragments(),
                           {RNAME: Reference(RNAME, REF)}, Parameters())
    old = {"evidence_ledger": result.evidence_ledger}
    resolved_fragment = next(
        record for record in result.evidence_ledger["fragment_records"]
        if record["qname"] == "shared"
    )
    new_ledger = {
        "fragment_records": [
            {**resolved_fragment, "ambiguity_groups": [],
             "ambiguous_between": False},
            *[
                record for record in result.evidence_ledger["fragment_records"]
                if record["qname"] != "shared"
            ],
        ]
    }
    reconciliation = _ambiguity_reconciliation(
        old, {"evidence_ledger": new_ledger})
    assert any("shared" in x for x in reconciliation["resolved_ambiguous_fragments"])
    assert any("persist" in x
               for x in reconciliation["persistent_ambiguous_fragments"])
