import json

from variant_pipeline.caller import call_variants
from variant_pipeline.config import Parameters
from variant_pipeline.fasta import Reference
from variant_pipeline.selection import adopt, load_payload, load_selections
from variant_pipeline.store import RunStore, sha256_file
from tests.helpers import REF, RNAME
from tests.test_fragment_provenance import _ambiguous_fragments


def test_selection_files_are_outside_immutable_runs(tmp_path):
    store = RunStore(tmp_path / "out")
    fragments = _ambiguous_fragments()
    result = call_variants(fragments, {RNAME: Reference(RNAME, REF)},
                           Parameters())
    payload = {
        "candidates": [
            {"variant_id": ev.variant_id, "kind": ev.kind,
             "supporting_evidence": ev.supporting_summary()}
            for ev in result.candidates
        ],
        "evidence_ledger": result.evidence_ledger,
    }
    chosen = next(c["variant_id"] for c in payload["candidates"]
                  if c["kind"] == "SNV")
    other = next(c["variant_id"] for c in payload["candidates"]
                 if c["kind"] == "INS")

    store.root.mkdir()
    run_dir = store.root / "runs" / "v0001"
    run_dir.mkdir(parents=True)
    (run_dir / "candidates.json").write_text(
        json.dumps(payload), encoding="utf-8")
    candidate_hash = sha256_file(run_dir / "candidates.json")
    manifest = {"version": "v0001",
                "outputs": {"candidates.json": candidate_hash}}
    (run_dir / "manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8")

    adopt(store, "v0001", chosen)
    before = load_selections(store)
    assert load_payload(store, "v0001")
    assert load_selections(store) == before
    try:
        adopt(store, "v0001", other)
    except Exception as exc:
        assert "manual choice already exists" in str(exc)
    assert load_selections(store) == before

    exported = tmp_path / "exported"
    store.export("v0001", exported)
    decision = json.loads((exported / "selected_interpretation.json")
                          .read_text(encoding="utf-8"))
    assert decision["version"] == "v0001"
    assert chosen in json.dumps(decision, ensure_ascii=False)
    assert (exported / "selection_snapshot.json").is_file()
