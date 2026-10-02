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
from variant_pipeline.samio import Fragment, Observation
from variant_pipeline.store import RunStore


PARAMS = Parameters(min_coverage_for_calling=1, min_alt_count=1)
REF = Reference("c", "ACGTACGT")


def _fragment(qname, source, allele_specs):
    fragment = Fragment(qname=qname, strand="F", source_id=source,
                        source_name=source, batch_id=source)
    for pos, alt in allele_specs:
        ref_base = REF.sequence[pos - 1]
        variant = normalize("c", pos, ref_base, alt, REF)
        fragment.alt[variant.canonical_key] = Observation(
            variant=variant, base_quals=(40,), strand="F", low_quality=False)
        fragment.covered.append((pos, pos))
    return fragment


def test_same_fragment_shared_by_candidates_is_not_independent_twice():
    fragments = [
        _fragment("shared", "batch-a", [(4, "A"), (4, "G")]),
        _fragment("only-a", "batch-a", [(4, "A")]),
        _fragment("only-b", "batch-a", [(4, "G")]),
    ]
    result = call_variants(fragments, {"c": REF}, PARAMS)
    by_id = {item.variant_id: item for item in result.candidates}
    a = by_id["c:g.4T>A"]
    b = by_id["c:g.4T>G"]
    assert a.alt_fragments == 2
    assert a.independent_alt_fragments == 1
    assert a.shared_ambiguous_alt_fragments == 1
    assert b.independent_alt_fragments == 1
    assert b.shared_ambiguous_alt_fragments == 1
    assert result.ambiguity_groups[0]["shared_fragments"] == [
        {
            "fragment_id": "batch-a::shared",
            "supports": ["c:g.4T>A", "c:g.4T>G"],
        }
    ]


def test_same_qname_from_different_sources_remains_distinct_evidence():
    fragments = [
        _fragment("same-name", "batch-a", [(4, "A")]),
        _fragment("same-name", "batch-b", [(4, "A")]),
    ]
    result = call_variants(fragments, {"c": REF}, PARAMS)
    candidate = result.candidates[0]
    assert candidate.alt_fragments == 2
    assert candidate.independent_alt_fragments == 2
    assert candidate.alt_fragment_ids == [
        "batch-a::same-name", "batch-b::same-name"]
    assert result.fragment_sources[0]["source_id"] == "batch-a"
    assert result.fragment_sources[1]["source_id"] == "batch-b"


def _publish_payload(store, tmp_path):
    from variant_pipeline.report import build_payload, write_json
    from variant_pipeline.config import RunConfig, ReferenceConfig, SampleConfig

    fragments = [
        _fragment("shared", "batch-a", [(4, "A"), (4, "G")]),
        _fragment("only-a", "batch-a", [(4, "A")]),
        _fragment("only-b", "batch-a", [(4, "G")]),
    ]
    result = call_variants(fragments, {"c": REF}, PARAMS)
    fasta = tmp_path / "ref.fasta"
    sam = tmp_path / "reads.sam"
    fasta.write_text(">c\nACGTACGT\n", encoding="utf-8")
    sam.write_text("@HD\n", encoding="utf-8")
    cfg = RunConfig(
        sample=SampleConfig(sample_id="s"),
        reference=ReferenceConfig(
            fasta_path=str(fasta), name="c", assembly="ref-1"),
        parameters=PARAMS,
    )
    payload = build_payload(cfg, result)
    staging = tmp_path / "staging"
    staging.mkdir()
    write_json(payload, staging / "candidates.json")
    for name in ("candidates.tsv", "coverage.tsv", "report.txt"):
        (staging / name).write_text("", encoding="utf-8")
    fp = store.fingerprint(config_dict=cfg.to_dict(), reads_path=sam,
                           fasta_path=fasta)
    return store.publish(
        {name: staging / name for name in
         ("candidates.json", "candidates.tsv", "coverage.tsv", "report.txt")},
        config_dict=cfg.to_dict(), reads_path=sam, fasta_path=fasta,
        fingerprint=fp,
    )


def test_selection_withdraw_and_export_audit(tmp_path):
    store = RunStore(tmp_path / "store")
    info = _publish_payload(store, tmp_path)
    before = (info.path / "candidates.json").read_bytes()
    save_selection(store, info.version, "c:g.4T>A")
    view = selected_view(store.verify(info.version), "c:g.4T>A")
    selected, other = [
        c for c in view["candidates"]
        if c["variant_id"] in ("c:g.4T>A", "c:g.4T>G")
    ]
    assert selected["supporting_evidence"][
        "provisional_independent_alt_fragments"] == 2
    assert other["supporting_evidence"][
        "provisional_independent_alt_fragments"] == 1
    assert (info.path / "candidates.json").read_bytes() == before
    assert load_selection(store, info.version)["status"] == "provisional"
    withdraw_selection(store, info.version)
    assert load_selection(store, info.version) is None
    assert (info.path / "candidates.json").read_bytes() == before

    dest = tmp_path / "export"
    store.export(info.version, dest)
    audit = json.loads((dest / "shared_evidence_audit.json").read_text())
    record = audit["ambiguity_groups"][0]["shared_fragments"][0]
    assert record["fragment_id"] == "batch-a::shared"
    assert record["source"]["source_name"] == "batch-a"


def test_cli_selection_is_side_record_and_export_includes_audit(tmp_path):
    from variant_pipeline.cli import main

    store = RunStore(tmp_path / "store")
    info = _publish_payload(store, tmp_path)
    original = (info.path / "candidates.json").read_bytes()
    assert main([
        "select", "--outdir", str(store.root),
        "--version", info.version, "--variant", "c:g.4T>A",
    ]) == 0
    assert main([
        "selected-view", "--outdir", str(store.root),
        "--version", info.version, "--json-out",
        str(tmp_path / "view.json"),
    ]) == 0
    assert main([
        "withdraw", "--outdir", str(store.root),
        "--version", info.version,
    ]) == 0
    assert (info.path / "candidates.json").read_bytes() == original
    assert main([
        "export", "--outdir", str(store.root),
        "--version", info.version, "--dest", str(tmp_path / "out"),
    ]) == 0
    assert (tmp_path / "out" / "shared_evidence_audit.json").is_file()
