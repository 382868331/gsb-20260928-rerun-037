"""Linked fragments from one original molecule as a single evidence unit."""

import json

import pytest

from variant_pipeline.caller import call_variants
from variant_pipeline.config import (
    Parameters,
    ReferenceConfig,
    RunConfig,
    SampleConfig,
)
from variant_pipeline.fasta import Reference
from variant_pipeline.interpretation import (
    link_view,
    list_link_selections,
    load_link_selection,
    save_link_selection,
    withdraw_link_selection,
)
from variant_pipeline.linkage import LinkSpec, recount_payload
from variant_pipeline.normalize import normalize
from variant_pipeline.report import build_payload
from variant_pipeline.samio import Fragment, Observation
from variant_pipeline.store import RunStore


PARAMS = Parameters(min_coverage_for_calling=1, min_alt_count=1)
REF = Reference("c", "ACGTACGTACGT")


def _fragment(qname, source, allele_specs, strand="F", span=None):
    frag = Fragment(qname=qname, strand=strand, source_id=source,
                    source_name=source, batch_id=source)
    frag.chroms = {"c"}
    for pos, alt in allele_specs:
        ref_base = REF.sequence[pos - 1]
        variant = normalize("c", pos, ref_base, alt, REF)
        frag.alt[variant.canonical_key] = Observation(
            variant=variant, base_quals=(40,), strand=strand,
            low_quality=False)
    if span is None:
        positions = [p for p, _ in allele_specs] or [1]
        span = (min(positions), max(positions))
    frag.covered.append(span)
    return frag


def _link(link_id, molecule, members, source=None, orientation="any",
          max_distance=None, same_chrom=False):
    return LinkSpec(
        link_id=link_id,
        molecule_id=molecule,
        members=tuple(members),
        source=source,
        expected_orientation=orientation,
        max_distance=max_distance,
        same_chrom=same_chrom,
    )


def _unit(result, link_id):
    return next(u for u in result.linkage_units if u["link_id"] == link_id)


def test_linked_fragments_are_one_molecule_not_two_supports():
    fragments = [
        _fragment("m1a", "batch-a", [(4, "A")], span=(1, 4)),
        _fragment("m1b", "batch-a", [(4, "A")], span=(4, 9)),
        _fragment("solo", "batch-a", [(4, "A")], span=(2, 5)),
    ]
    links = [_link("L1", "mol-1", ["m1a", "m1b"], source="batch-a")]
    result = call_variants(fragments, {"c": REF}, PARAMS, links=links)
    ev = result.candidates[0]
    # Fragment-level record is unchanged...
    assert ev.alt_fragments == 3
    assert ev.independent_alt_fragments == 3
    # ...but the two linked fragments are one molecule, not two supports.
    assert ev.alt_molecules == 2
    assert ev.independent_alt_molecules == 2
    assert _unit(result, "L1")["state"] == "complete"


def test_combination_of_linked_matches_excludes_interpretations():
    shared = _fragment("m2a", "batch-a", [(4, "A"), (4, "G")], span=(1, 4))
    partner = _fragment("m2b", "batch-a", [(4, "A")], span=(4, 9))
    links = [_link("L2", "mol-2", ["m2a", "m2b"])]
    result = call_variants([shared, partner], {"c": REF}, PARAMS, links=links)
    by_id = {c.variant_id: c for c in result.candidates}
    a = by_id["c:g.4T>A"]
    g = by_id["c:g.4T>G"]
    # Fragment-level ambiguity is still reported as before...
    assert a.shared_ambiguous_alt_fragments == 1
    assert g.shared_ambiguous_alt_fragments == 1
    # ...but the combination keeps only the interpretation both ends share.
    assert a.linkage_resolved_molecules == 1
    assert a.independent_alt_molecules == 1
    assert g.independent_alt_molecules == 0
    assert g.shared_ambiguous_molecules == 0


def test_disjoint_linked_matches_are_a_molecule_conflict():
    first = _fragment("m3a", "batch-a", [(4, "A")], span=(1, 4))
    second = _fragment("m3b", "batch-a", [(4, "G")], span=(4, 9))
    links = [_link("L3", "mol-3", ["m3a", "m3b"])]
    result = call_variants([first, second], {"c": REF}, PARAMS, links=links)
    assert result.molecule_conflicts
    conflict = result.molecule_conflicts[0]
    assert conflict["reason"] == "members_support_disjoint_candidates"
    assert set(conflict["candidate_variant_ids"]) == {"c:g.4T>A", "c:g.4T>G"}


def test_pending_link_keeps_single_end_as_one_explainable_unit():
    fragments = [_fragment("m4a", "batch-a", [(4, "A")], span=(1, 4))]
    links = [_link("L4", "mol-4", ["m4a", "m4b"])]
    result = call_variants(fragments, {"c": REF}, PARAMS, links=links)
    unit = _unit(result, "L4")
    assert unit["state"] == "pending"
    assert unit["missing_members"] == ["m4b"]
    assert unit["reasons"][0]["kind"] == "fragment_matching_insufficient"
    # The partial unit is still exactly one molecule of evidence.
    assert result.candidates[0].alt_molecules == 1

    # The other end arrives later (new analysis): now complete, and the
    # molecule still counts once -- late linkage never inflates support.
    fragments.append(_fragment("m4b", "batch-a", [(4, "A")], span=(4, 9)))
    result2 = call_variants(fragments, {"c": REF}, PARAMS, links=links)
    assert _unit(result2, "L4")["state"] == "complete"
    assert result2.candidates[0].alt_molecules == 1
    assert result2.candidates[0].alt_fragments == 2


def test_same_qname_from_different_source_is_not_auto_paired():
    frag_a = _fragment("dup", "batch-a", [(4, "A")], span=(1, 4))
    frag_b = _fragment("dup", "batch-b", [(4, "A")], span=(4, 9))
    links = [_link("L5", "mol-5", ["dup", "other"])]
    result = call_variants([frag_a, frag_b], {"c": REF}, PARAMS, links=links)
    unit = _unit(result, "L5")
    assert unit["state"] == "ambiguous_identity"
    assert unit["reasons"][0]["kind"] == "identity_ambiguous"
    # Neither end is merged or dropped: both stay fragment-level evidence.
    ev = result.candidates[0]
    assert ev.alt_fragments == 2
    assert ev.independent_alt_fragments == 2
    assert ev.alt_molecules == 2

    # Declaring the source pairs exactly the intended fragment.
    links = [_link("L5", "mol-5", ["dup", "other"], source="batch-a")]
    result = call_variants([frag_a, frag_b], {"c": REF}, PARAMS, links=links)
    unit = _unit(result, "L5")
    assert unit["state"] == "pending"
    assert [m["fragment_id"] for m in unit["members"]] == ["batch-a::dup"]


def test_contradictory_distance_is_reported_and_both_ends_kept():
    frag_a = _fragment("m6a", "batch-a", [(4, "A")], span=(1, 4))
    frag_b = _fragment("m6b", "batch-a", [], span=(10, 12))
    links = [_link("L6", "mol-6", ["m6a", "m6b"], max_distance=2)]
    result = call_variants([frag_a, frag_b], {"c": REF}, PARAMS, links=links)
    unit = _unit(result, "L6")
    assert unit["state"] == "inconsistent"
    kinds = {r["kind"] for r in unit["reasons"]}
    assert kinds == {"linkage_evidence_inconsistent"}
    # The inconvenient end is not silently dropped: both members stay
    # listed on the unit, and the alt end keeps its own evidence.
    assert {m["fragment_id"] for m in unit["members"]} == {
        "batch-a::m6a", "batch-a::m6b"}
    ev = result.candidates[0]
    assert ev.alt_fragments == 1
    assert ev.independent_alt_fragments == 1
    assert ev.alt_molecules == 1


def test_contradictory_orientation_is_reported():
    frag_a = _fragment("m7a", "batch-a", [(4, "A")], strand="F", span=(1, 4))
    frag_b = _fragment("m7b", "batch-a", [(4, "A")], strand="F", span=(4, 9))
    links = [_link("L7", "mol-7", ["m7a", "m7b"], orientation="FR")]
    result = call_variants([frag_a, frag_b], {"c": REF}, PARAMS, links=links)
    unit = _unit(result, "L7")
    assert unit["state"] == "inconsistent"
    assert "orientation" in unit["reasons"][0]["detail"]
    # Both ends remain fragment-level evidence despite the contradiction.
    ev = result.candidates[0]
    assert ev.alt_fragments == 2
    assert ev.independent_alt_fragments == 2
    assert ev.alt_molecules == 2


def test_duplicate_molecule_import_merges_into_one_unit():
    frag_a = _fragment("m8a", "batch-a", [(4, "A")], span=(1, 4))
    frag_b = _fragment("m8b", "batch-a", [(4, "A")], span=(4, 9))
    links = [
        _link("L8a", "mol-8", ["m8a"]),
        _link("L8b", "mol-8", ["m8b"]),
    ]
    result = call_variants([frag_a, frag_b], {"c": REF}, PARAMS, links=links)
    units = [u for u in result.linkage_units if u["molecule_id"] == "mol-8"]
    assert len(units) == 1
    assert units[0]["reasons"][0]["kind"] == "duplicate_import"
    assert result.candidates[0].alt_molecules == 1


def _publish_linked_payload(store, tmp_path, links, name="case"):
    from variant_pipeline.report import write_json

    fragments = [
        _fragment("m2a", "batch-a", [(4, "A"), (4, "G")], span=(1, 4)),
        _fragment("m2b", "batch-a", [(4, "A")], span=(4, 9)),
    ]
    result = call_variants(fragments, {"c": REF}, PARAMS, links=links)
    case_dir = tmp_path / name
    case_dir.mkdir(exist_ok=True)
    fasta = case_dir / "ref.fasta"
    sam = case_dir / "reads.sam"
    fasta.write_text(">c" + chr(10) + REF.sequence + chr(10), encoding="utf-8")
    sam.write_text("@HD" + chr(10), encoding="utf-8")
    cfg = RunConfig(
        sample=SampleConfig(sample_id="s"),
        reference=ReferenceConfig(
            fasta_path=str(fasta), name="c", assembly="ref-1"),
        parameters=PARAMS,
    )
    payload = build_payload(cfg, result)
    staging = case_dir / "staging"
    staging.mkdir()
    write_json(payload, staging / "candidates.json")
    for fname in ("candidates.tsv", "coverage.tsv", "report.txt"):
        (staging / fname).write_text("", encoding="utf-8")
    fp = store.fingerprint(config_dict=cfg.to_dict(), reads_path=sam,
                           fasta_path=fasta)
    return store.publish(
        {fname: staging / fname for fname in
         ("candidates.json", "candidates.tsv", "coverage.tsv", "report.txt")},
        config_dict=cfg.to_dict(), reads_path=sam, fasta_path=fasta,
        fingerprint=fp,
    ), payload


def test_link_selection_keeps_excluded_candidates_and_rationale(tmp_path):
    store = RunStore(tmp_path / "store")
    links = [_link("L2", "mol-2", ["m2a", "m2b"])]
    info, _payload = _publish_linked_payload(store, tmp_path, links)
    record = save_link_selection(
        store, info.version, "L2", "c:g.4T>A", "关联组合仅共同支持 A")
    assert record["status"] == "provisional"
    assert record["excluded_variant_ids"] == ["c:g.4T>G"]
    assert record["rationale"] == "关联组合仅共同支持 A"
    with pytest.raises(ValueError):
        save_link_selection(store, info.version, "L2", "c:g.4T>G", "覆盖")
    with pytest.raises(ValueError):
        save_link_selection(
            store, info.version, "L2", "c:g.1A>C", "不在该关联单元中")

    view = link_view(store.verify(info.version), "L2", adoption=record)
    states = {c["variant_id"]: c.get("selection_state")
              for c in view["candidates"]}
    assert states["c:g.4T>A"] == "provisionally_selected"
    assert states["c:g.4T>G"] == "excluded_by_provisional_choice"
    assert view["does_not_modify_base_version"] is True

    withdrawn = withdraw_link_selection(store, info.version, "L2")
    assert withdrawn["status"] == "withdrawn"
    assert load_link_selection(store, info.version, "L2") is None
    # Withdrawal restores the original ambiguity view; nothing is invented.
    restored = link_view(store.verify(info.version), "L2", adoption=None)
    assert all(c.get("selection_state") is None
               for c in restored["candidates"])
    # History (including the withdrawn record) is preserved.
    history = list_link_selections(store, info.version)
    assert len(history) == 1 and history[0]["status"] == "withdrawn"


def test_new_version_does_not_overwrite_old_link_selection(tmp_path):
    store = RunStore(tmp_path / "store")
    links = [_link("L2", "mol-2", ["m2a", "m2b"])]
    info1, _ = _publish_linked_payload(store, tmp_path, links, name="c1")
    save_link_selection(store, info1.version, "L2", "c:g.4T>A", "第一版选择")
    info2, _ = _publish_linked_payload(store, tmp_path, links, name="c2")
    assert info2.version != info1.version
    old = load_link_selection(store, info1.version, "L2")
    assert old["status"] == "provisional"
    assert old["selected_variant_id"] == "c:g.4T>A"
    assert load_link_selection(store, info2.version, "L2") is None


def test_export_linkage_audit_traces_molecule_to_fragments(tmp_path):
    store = RunStore(tmp_path / "store")
    links = [_link("L2", "mol-2", ["m2a", "m2b"])]
    info, _ = _publish_linked_payload(store, tmp_path, links)
    save_link_selection(store, info.version, "L2", "c:g.4T>A", "采用A")
    withdraw_link_selection(store, info.version, "L2")
    dest = tmp_path / "export"
    store.export(info.version, dest)
    audit = json.loads((dest / "linkage_audit.json").read_text())
    unit = next(u for u in audit["linkage_units"] if u["link_id"] == "L2")
    member_ids = {m["fragment_id"] for m in unit["members"]}
    assert member_ids == {"batch-a::m2a", "batch-a::m2b"}
    assert audit["link_adoptions"][0]["status"] == "withdrawn"
    mol = audit["candidate_molecule_evidence"]["c:g.4T>A"]
    assert mol["linkage_resolved_molecules"] == 1
    assert mol["molecule_ids"] == ["unit-mol-2"]


def test_recount_under_alternative_linkage_rule(tmp_path):
    store = RunStore(tmp_path / "store")
    links = [_link("L2", "mol-2", ["m2a", "m2b"])]
    _info, payload = _publish_linked_payload(store, tmp_path, links)
    # Alternative rule: no linkage at all -> both fragments count alone and
    # the shared ambiguity returns instead of the resolved molecule.
    alternative = recount_payload(payload, [])
    counts = alternative["molecule_counts"]
    assert counts["c:g.4T>A"]["independent_alt_molecules"] == 1
    assert counts["c:g.4T>A"]["shared_ambiguous_molecules"] == 1
    assert counts["c:g.4T>A"]["linkage_resolved_molecules"] == 0
    # Stored rule: the combination resolved the molecule to A.
    stored = payload["candidates"]
    a = next(c for c in stored if c["variant_id"] == "c:g.4T>A")
    assert a["molecule_evidence"]["linkage_resolved_molecules"] == 1


def test_cli_linkage_workflow(tmp_path):
    import json as _json

    from variant_pipeline.cli import main
    from tests.helpers import RNAME, substitution_read, write_case

    lines = []
    lines += [substitution_read(f"a{i}", 6 + i, 20, 15, "A")
              for i in range(2)]
    lines += [substitution_read(f"r{i}", 8 + i, 20, 15, "A", reverse=True)
              for i in range(2)]
    lines += [substitution_read(f"ref{i}", 7 + i, 20, 15)
              for i in range(3)]
    fasta, sam = write_case(tmp_path, lines)
    config = {
        "sample": {"sample_id": "S1"},
        "reference": {
            "fasta_path": str(fasta), "name": RNAME, "assembly": "synth-v1"},
        "reads_path": str(sam),
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(_json.dumps(config), encoding="utf-8")
    links_path = tmp_path / "links.json"
    links_path.write_text(_json.dumps({"links": [
        {"link_id": "L1", "molecule_id": "mol-1",
         "members": ["a0", "a1"],
         "expected": {"max_distance": 100}},
        {"link_id": "L2", "molecule_id": "mol-2",
         "members": ["r0", "late1"]},
    ]}), encoding="utf-8")
    empty_links = tmp_path / "links-none.json"
    empty_links.write_text(_json.dumps({"links": []}), encoding="utf-8")
    outdir = tmp_path / "out"

    assert main(["call", "--config", str(config_path),
                 "--outdir", str(outdir), "--links", str(links_path)]) == 0
    payload = _json.loads(
        (outdir / "runs" / "v0001" / "candidates.json").read_text())
    cand = payload["candidates"][0]
    assert cand["supporting_evidence"]["alt_fragments"] == 4
    # a0+a1 are one molecule; r0 is a pending unit; r1 is a singleton.
    assert cand["molecule_evidence"]["alt_molecules"] == 3
    states = {u["link_id"]: u["state"] for u in payload["linkage_units"]
              if u["link_id"]}
    assert states == {"L1": "complete", "L2": "pending"}

    vid = "synth-1:g.15G>A"
    assert main(["link-view", "--outdir", str(outdir),
                 "--version", "v0001", "--link", "L1"]) == 0
    assert main(["link-select", "--outdir", str(outdir),
                 "--version", "v0001", "--link", "L1",
                 "--variant", vid, "--rationale", "两端共同支持"]) == 0
    assert main(["link-compare", "--outdir", str(outdir),
                 "--version", "v0001", "--links", str(empty_links)]) == 0
    assert main(["link-withdraw", "--outdir", str(outdir),
                 "--version", "v0001", "--link", "L1"]) == 0
    dest = tmp_path / "exp"
    assert main(["export", "--outdir", str(outdir),
                 "--version", "v0001", "--dest", str(dest)]) == 0
    audit = _json.loads((dest / "linkage_audit.json").read_text())
    assert audit["link_adoptions"][0]["status"] == "withdrawn"
    assert audit["link_adoptions"][0]["rationale"] == "两端共同支持"
    assert main(["verify", "--outdir", str(outdir)]) == 0
