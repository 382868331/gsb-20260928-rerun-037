"""Versioned batches, cross-version diff and frozen historical export.

End-to-end tests through the CLI so the storage contract is exercised the
same way users invoke it. Each run reads/writes only its own batch and the
shared index; historical batch files are never re-analyzed or rewritten.
"""

from __future__ import annotations

import json

from variant_pipeline import ANALYSIS_VERSION
from variant_pipeline.cli import main as cli

from tests.helpers import (
    REF,
    RNAME,
    substitution_read,
)


def _setup_inputs(base):
    fasta = base / "reference.fasta"
    fasta.write_text(f">{RNAME}\n{REF}\n", encoding="utf-8")
    config = base / "config.json"
    config.write_text(
        json.dumps(
            {
                "sample": {"sample_id": "S1"},
                "reference": {
                    "fasta_path": "reference.fasta",
                    "name": RNAME,
                    "assembly": "asm-v1",
                },
                "reads_path": "reads.sam",
            }
        ),
        encoding="utf-8",
    )
    return fasta, config


def _write_reads(base, pos, alt="A", n_alt=4, n_ref=3):
    lines = [
        f"@HD\tVN:1.6\n",
        f"@SQ\tSN:{RNAME}\tLN:{len(REF)}\n",
    ]
    for i in range(n_alt // 2):
        lines.append(substitution_read(f"f{i}", 6 + i, 20, pos, alt))
        lines.append(substitution_read(f"r{i}", 8 + i, 20, pos, alt, reverse=True))
    for i in range(n_ref):
        lines.append(substitution_read(f"ref{i}", 7 + i, 20, pos))
    (base / "reads.sam").write_text("".join(lines), encoding="utf-8")


def _call(outdir, config_path):
    rc = cli(["call", "--config", str(config_path), "--outdir", str(outdir)])
    assert rc == 0


def _batch_ids(outdir):
    index = json.loads((outdir / "index.json").read_text(encoding="utf-8"))
    return [b["batch_id"] for b in index["batches"]]


def _payload(outdir, batch_id):
    path = outdir / "batches" / batch_id / "candidates.json"
    return json.loads(path.read_text(encoding="utf-8"))


def test_same_output_root_keeps_distinct_immutable_batches(tmp_path):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    first = _batch_ids(out)[0]
    first_json = (out / "batches" / first / "candidates.json").read_bytes()
    first_report = (out / "batches" / first / "report.txt").read_bytes()

    # New inputs analyzed into the same output root: old batch must survive.
    _write_reads(tmp_path, pos=19)
    _call(out, tmp_path / "config.json")
    ids = _batch_ids(out)
    assert len(ids) == 2 and first == ids[0]
    assert (out / "batches" / first / "candidates.json").read_bytes() == first_json
    assert (out / "batches" / first / "report.txt").read_bytes() == first_report

    # The old batch still describes the original position 15 variant.
    old = _payload(out, first)
    assert old["candidates"][0]["variant_id"] == f"{RNAME}:g.15G>A"
    new = _payload(out, ids[1])
    assert new["candidates"][0]["variant_id"] == f"{RNAME}:g.19G>A"
    assert new["batch"]["input_hashes"] != old["batch"]["input_hashes"]


def test_duplicate_request_creates_new_batch_without_touching_old(tmp_path):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    first = _batch_ids(out)[0]
    first_files = {
        p.relative_to(out / "batches" / first): p.read_bytes()
        for p in (out / "batches" / first).rglob("*")
        if p.is_file()
    }
    _call(out, tmp_path / "config.json")
    ids = _batch_ids(out)
    assert len(ids) == 2
    second = ids[1]
    assert first != second
    index = json.loads((out / "index.json").read_text(encoding="utf-8"))
    assert index["batches"][1]["duplicate_of"] == [first]
    for rel, data in first_files.items():
        assert (out / "batches" / first / rel).read_bytes() == data
    # Duplicate outputs are identical analysis content, stored separately.
    assert _payload(out, first)["candidates"] == _payload(out, second)["candidates"]


def test_input_snapshot_is_preserved_per_batch(tmp_path):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    batch_id = _batch_ids(out)[0]
    batch = out / "batches" / batch_id
    assert (batch / "input" / "reads.sam").read_bytes() == (
        tmp_path / "reads.sam"
    ).read_bytes()
    assert (batch / "input" / "reference.fasta").read_bytes() == (
        tmp_path / "reference.fasta"
    ).read_bytes()
    manifest = json.loads((batch / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["analysis_version"] == ANALYSIS_VERSION
    assert manifest["content_id"]
    assert manifest["files"]["input/reads.sam"]
    # Soft-clipped/shared-fragment evidence lives only in this batch's SAM.
    assert "input/reads.sam" in manifest["files"]


def test_diff_attributes_appeared_disappeared_to_input_change(tmp_path):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    first = _batch_ids(out)[0]
    _write_reads(tmp_path, pos=19)
    _call(out, tmp_path / "config.json")
    second = _batch_ids(out)[1]

    rc = cli(["diff", "--outdir", str(out), "--from", first,
              "--to", second, "--json"])
    assert rc == 0
    # Text rendering must also run without re-analysis.
    assert cli(["diff", "--outdir", str(out), "--to", second]) == 0

    from variant_pipeline.diff import diff_batches

    result = diff_batches(out, first, second)
    assert result["inputs_changed"] is True
    assert result["interpretation_rules_changed"] is False
    assert [e["variant_id"] for e in result["removed"]] == [
        f"{RNAME}:g.15G>A"
    ]
    assert [e["variant_id"] for e in result["added"]] == [f"{RNAME}:g.19G>A"]
    for entry in result["removed"] + result["added"]:
        assert entry["attribution"] == ["input_change"]


def test_diff_tracks_support_count_changes_from_inputs(tmp_path):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15, n_alt=4, n_ref=3)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    first = _batch_ids(out)[0]
    _write_reads(tmp_path, pos=15, n_alt=4, n_ref=9)
    _call(out, tmp_path / "config.json")
    second = _batch_ids(out)[1]

    from variant_pipeline.diff import diff_batches

    result = diff_batches(out, first, second)
    changed = {e["variant_id"]: e for e in result["changed"]}
    vid = f"{RNAME}:g.15G>A"
    entry = changed[vid]
    fields = {c["field"]: c for c in entry["field_changes"]}
    assert fields["reference_fragments"]["from"] == 3
    assert fields["reference_fragments"]["to"] == 9
    assert entry["from"]["alt_fragments"] == 4
    assert entry["to"]["alt_fragments"] == 4
    # Side-by-side snapshots, never spliced support records.
    assert entry["from"]["variant_id"] == entry["to"]["variant_id"] == vid


def test_export_and_verify_use_stored_bytes_only(tmp_path, capsys):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    batch_id = _batch_ids(out)[0]
    original = (out / "batches" / batch_id / "candidates.json").read_bytes()

    # A later batch with different inputs.
    _write_reads(tmp_path, pos=19)
    _call(out, tmp_path / "config.json")

    export_dir = tmp_path / "exported"
    rc = cli(["export", "--outdir", str(out), "--batch", batch_id,
              "--dest", str(export_dir)])
    assert rc == 0
    assert (export_dir / "candidates.json").read_bytes() == original
    assert (export_dir / "input" / "reads.sam").read_bytes() == (
        out / "batches" / batch_id / "input" / "reads.sam"
    ).read_bytes()

    assert cli(["verify", "--outdir", str(out), "--batch", batch_id]) == 0
    # Tampering is detected and never repaired by re-analyzing.
    target = out / "batches" / batch_id / "candidates.json"
    target.write_text(original.decode("utf-8").replace("S1", "S1-TAMPERED"))
    assert cli(["verify", "--outdir", str(out), "--batch", batch_id]) == 1


def test_latest_prefix_resolution_and_listing(tmp_path):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    _write_reads(tmp_path, pos=19)
    _call(out, tmp_path / "config.json")
    ids = _batch_ids(out)
    assert cli(["list", "--outdir", str(out)]) == 0
    from variant_pipeline.store import resolve_batch_id

    assert resolve_batch_id(out, "latest") == ids[1]
    assert resolve_batch_id(out, ids[0][:20]) == ids[0]


def test_publish_failure_leaves_existing_batches_intact(tmp_path, monkeypatch):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    first = _batch_ids(out)[0]
    first_bytes = {
        p.relative_to(out / "batches" / first): p.read_bytes()
        for p in (out / "batches" / first).rglob("*")
        if p.is_file()
    }
    index_bytes = (out / "index.json").read_bytes()

    # New inputs, but staging fails after writing partial files.
    _write_reads(tmp_path, pos=19)
    from variant_pipeline import store

    real_replace = store.os.replace

    def fail_replace(src, dst):
        if str(src).startswith(str(out / "batches")):
            raise OSError("simulated publish failure")
        return real_replace(src, dst)

    monkeypatch.setattr(store.os, "replace", fail_replace)
    try:
        cli(["call", "--config", str(tmp_path / "config.json"),
             "--outdir", str(out)])
        raised = False
    except OSError:
        raised = True
    assert raised

    # Index untouched, existing batch untouched, no staged junk left behind.
    assert (out / "index.json").read_bytes() == index_bytes
    for rel, data in first_bytes.items():
        assert (out / "batches" / first / rel).read_bytes() == data
    leftovers = [
        p.name for p in (out / "batches").iterdir()
        if p.name.startswith(".stage_")
    ]
    assert leftovers == []


def test_diff_attributes_differences_to_rule_change_when_inputs_identical(tmp_path):
    _setup_inputs(tmp_path)
    _write_reads(tmp_path, pos=15)
    out = tmp_path / "out"
    _call(out, tmp_path / "config.json")
    first = _batch_ids(out)[0]

    # Simulate re-interpretation under new orientation rules: same input
    # snapshot bytes, but the stored historical result carries a different
    # analysis_version and different support counts/status.
    bdir = out / "batches" / first
    second_id = f"{first}_reinterpreted"
    second_dir = out / "batches" / second_id
    import shutil

    shutil.copytree(bdir, second_dir)
    manifest_path = second_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["batch_id"] = second_id
    manifest["analysis_version"] = "orientation-forward-v3"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    json_path = second_dir / "candidates.json"
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    payload["analysis_version"] = "orientation-forward-v3"
    payload["candidates"][0]["status"] = "REJECTED"
    payload["candidates"][0]["supporting_evidence"]["alt_fragments"] = 0
    json_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    from variant_pipeline.diff import diff_batches

    result = diff_batches(out, first, second_id)
    assert result["inputs_changed"] is False
    assert result["interpretation_rules_changed"] is True
    entry = result["changed"][0]
    assert entry["attribution"] == ["interpretation_rule_change"]
    changed_fields = {c["field"]: c for c in entry["field_changes"]}
    assert changed_fields["alt_fragments"]["from"] == 4
    assert changed_fields["alt_fragments"]["to"] == 0
    assert changed_fields["status"]["to"] == "REJECTED"
