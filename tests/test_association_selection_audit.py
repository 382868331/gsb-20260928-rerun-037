import json

from variant_pipeline.caller import call_variants
from variant_pipeline.config import Parameters
from variant_pipeline.fasta import Reference
from variant_pipeline.interpretation import (
    load_selection,
    save_selection,
    selected_view,
    withdraw_selection,
)
from variant_pipeline.normalize import normalize
from variant_pipeline.report import build_payload, write_json
from variant_pipeline.samio import Fragment, Observation
from variant_pipeline.association import AssociationRelation
from variant_pipeline.store import RunStore


PARAMS = Parameters(min_coverage_for_calling=1, min_alt_count=1)
REF = Reference("c", "ACGTACGT")


def _fragment(qname, source, specs, strand="F"):
    frag = Fragment(qname=qname, strand=strand, source_id=source,
                    source_name=source, batch_id=source)
    for pos, alt in specs:
        variant = normalize("c", pos, REF.sequence[pos - 1], alt, REF)
        frag.alt[variant.canonical_key] = Observation(
            variant=variant, base_quals=(40,), strand=strand, low_quality=False)
        frag.covered.append((pos, pos))
    return frag


def test_selection_history_and_export_trace_to_molecule(tmp_path):
    from variant_pipeline.config import ReferenceConfig, RunConfig, SampleConfig

    store = RunStore(tmp_path / "store")
    fragments = [
        _fragment("a", "s1", [(4, "A"), (8, "A")], "F"),
        _fragment("b", "s1", [(4, "G"), (8, "A")], "R"),
    ]
    relation = AssociationRelation(
        "unit-1", "mol-1", "s1::a", "s1::b", "s1", "s1",
        orientation="opposite", min_distance=4, max_distance=4,
        relation_id="rule-1")
    result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
    cfg = RunConfig(
        sample=SampleConfig(sample_id="s"),
        reference=ReferenceConfig(
            fasta_path=str(tmp_path / "ref.fasta"), name="c", assembly="r1"),
        parameters=PARAMS,
    )
    (tmp_path / "ref.fasta").write_text(">c\nACGTACGT\n")
    reads = tmp_path / "reads.sam"
    reads.write_text("@HD\n")
    staging = tmp_path / "stage"
    staging.mkdir()
    write_json(build_payload(cfg, result), staging / "candidates.json")
    for name in ("candidates.tsv", "coverage.tsv", "report.txt"):
        (staging / name).write_text("")
    fp = store.fingerprint(config_dict=cfg.to_dict(), reads_path=reads,
                           fasta_path=tmp_path / "ref.fasta")
    info = store.publish(
        {n: staging / n for n in (
            "candidates.json", "candidates.tsv", "coverage.tsv", "report.txt")},
        config_dict=cfg.to_dict(), reads_path=reads,
        fasta_path=tmp_path / "ref.fasta", fingerprint=fp)

    # First choice picks the linked, compatible c:g.8T>A allele.
    save_selection(store, info.version, "c:g.4T>A", "unit-1", "rule-1", "first")
    first_decision = load_selection(store, info.version)
    assert first_decision["excluded_variant_ids"] == ["c:g.8T>A"]
    assert first_decision["conflicting_variant_ids"] == [
        "c:g.4T>A", "c:g.4T>G"]
    view = selected_view(store.verify(info.version), "c:g.4T>A", first_decision)
    view_by_id = {c["variant_id"]: c for c in view["candidates"]}
    # The other same-site allele is an intra-molecule conflict, never silently
    # dropped or marked as excluded by import order; distance-excluded 8T>A
    # is marked under the provisional relation choice.
    assert view_by_id["c:g.4T>G"]["selection_state"] == "intra_molecule_conflict"
    assert view_by_id["c:g.8T>A"]["selection_state"] == "excluded_by_provisional_choice"
    assert view_by_id["c:g.4T>A"]["supporting_evidence"]["alt_fragments"] == 0
    # Swap the provisional choice to the other end's allele: conflict basis,
    # excluded combinations and support counts stay identical and consistent.
    save_selection(store, info.version, "c:g.4T>G", "unit-1", "rule-2", "second")
    decision = load_selection(store, info.version)
    assert decision["selected_variant_id"] == "c:g.4T>G"
    assert decision["supersedes"] is not None
    assert len(decision["history"]) == 1
    assert decision["history"][0]["selected_variant_id"] == "c:g.4T>A"
    assert decision["excluded_variant_ids"] == ["c:g.8T>A"]
    assert decision["conflicting_variant_ids"] == [
        "c:g.4T>A", "c:g.4T>G"]
    conflict = decision["conflicting_combinations"][0]
    assert conflict["a_variant_id"] == "c:g.4T>A"
    assert conflict["b_variant_id"] == "c:g.4T>G"
    view2 = selected_view(store.verify(info.version), "c:g.4T>G", decision)
    view2_by_id = {c["variant_id"]: c for c in view2["candidates"]}
    assert view2_by_id["c:g.4T>A"]["selection_state"] == "intra_molecule_conflict"
    assert view2_by_id["c:g.4T>G"]["supporting_evidence"]["alt_fragments"] == 0
    withdraw_selection(store, info.version)
    assert load_selection(store, info.version) is None

    dest = tmp_path / "export"
    store.export(info.version, dest)
    audit = json.loads((dest / "shared_evidence_audit.json").read_text())
    assert audit["evidence_units"][0]["molecule_id"] == "mol-1"
    assert audit["evidence_units"][0]["fragment_ids"] == ["s1::a", "s1::b"]
    assert audit["provisional_selection"]["status"] == "withdrawn"
