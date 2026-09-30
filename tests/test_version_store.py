"""Versioned run store: immutability, dedup, diff, export, verify."""

import json
from pathlib import Path

import pytest

from variant_pipeline.cli import main
from variant_pipeline.diff import compare_versions
from variant_pipeline.store import (
    OUTPUT_FILENAMES,
    RunStore,
    VerificationError,
)
from tests.helpers import REF, RNAME, substitution_read, write_case


def _make_case(tmp_path, alt="A", alt_count=4, name="case"):
    case_dir = tmp_path / name
    case_dir.mkdir()
    pos = 15
    lines = []
    lines += [substitution_read(f"a{i}", 6 + i, 20, pos, alt, reverse=False)
              for i in range(alt_count // 2)]
    lines += [substitution_read(f"r{i}", 8 + i, 20, pos, alt, reverse=True)
              for i in range(alt_count - alt_count // 2)]
    lines += [substitution_read(f"ref{i}", 7 + i, 20, pos) for i in range(3)]
    fasta, sam = write_case(case_dir, lines)
    config = {
        "sample": {"sample_id": "S1"},
        "reference": {
            "fasta_path": str(fasta),
            "name": RNAME,
            "assembly": "synth-v1",
        },
        "reads_path": str(sam),
    }
    config_path = case_dir / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    return config_path


def _call(config_path, outdir, *extra):
    argv = ["call", "--config", str(config_path), "--outdir", str(outdir)]
    argv += list(extra)
    assert main(argv) == 0


def _snapshot(version_dir: Path) -> dict:
    return {
        name: (version_dir / name).read_bytes()
        for name in OUTPUT_FILENAMES
    }


def test_same_outdir_creates_distinct_versions_and_keeps_old(tmp_path):
    outdir = tmp_path / "out"
    _call(_make_case(tmp_path, alt="A", name="c1"), outdir)
    before = _snapshot(outdir / "runs" / "v0001")
    _call(_make_case(tmp_path, alt="T", name="c2"), outdir)
    assert (outdir / "runs" / "v0002").is_dir()
    # Old batch was not overwritten or recomputed in place.
    assert _snapshot(outdir / "runs" / "v0001") == before
    # No top-level candidates.json clobbering previous results.
    assert not (outdir / "candidates.json").exists()


def test_duplicate_request_returns_existing_version(tmp_path, capsys):
    outdir = tmp_path / "out"
    config = _make_case(tmp_path, alt="A", name="c1")
    _call(config, outdir)
    _call(config, outdir)
    versions = [i.version for i in RunStore(outdir).list_versions()]
    assert versions == ["v0001"]
    assert "duplicate" in capsys.readouterr().out.lower()


def test_force_duplicate_publishes_new_version(tmp_path):
    outdir = tmp_path / "out"
    config = _make_case(tmp_path, alt="A", name="c1")
    _call(config, outdir)
    _call(config, outdir, "--force")
    versions = [i.version for i in RunStore(outdir).list_versions()]
    assert versions == ["v0001", "v0002"]


def test_failed_publish_preserves_existing_versions(tmp_path, monkeypatch):
    outdir = tmp_path / "out"
    _call(_make_case(tmp_path, alt="A", name="c1"), outdir)
    before = _snapshot(outdir / "runs" / "v0001")

    import variant_pipeline.store as store_mod

    def boom(src, dst):
        raise OSError("simulated publish failure")

    monkeypatch.setattr(store_mod.os, "rename", boom)
    with pytest.raises(OSError):
        _call(_make_case(tmp_path, alt="T", name="c2"), outdir)
    store = RunStore(outdir)
    assert [i.version for i in store.list_versions()] == ["v0001"]
    assert _snapshot(outdir / "runs" / "v0001") == before
    # No staging leftovers after the failed publish.
    leftovers = [p for p in (outdir / "runs").iterdir()
                 if p.name.startswith(".staging-")]
    assert leftovers == []
    store.verify("v0001")


def test_compare_attributes_input_change(tmp_path):
    outdir = tmp_path / "out"
    _call(_make_case(tmp_path, alt="A", alt_count=4, name="c1"), outdir)
    _call(_make_case(tmp_path, alt="A", alt_count=6, name="c2"), outdir)
    diff = compare_versions(RunStore(outdir), "v0001", "v0002")
    assert diff["attribution"]["cause"] == "input_change"
    assert diff["attribution"]["input_changed"] is True
    assert diff["attribution"]["rule_changed"] is False
    assert diff["summary"]["changed"] == 1
    entry = diff["changed"][0]
    # Per-version evidence kept separate, never summed across versions.
    assert entry["from"]["support"]["alt_fragments"] == 4
    assert entry["to"]["support"]["alt_fragments"] == 6
    assert entry["changes"]["support.alt_fragments"] == {"from": 4, "to": 6}


def test_compare_attributes_rule_change(tmp_path):
    outdir = tmp_path / "out"
    config1 = _make_case(tmp_path, alt="A", name="c1")
    _call(config1, outdir)
    # Same reads, stricter min_alt_fraction -> rule change only.
    data = json.loads(config1.read_text(encoding="utf-8"))
    data["parameters"] = {"min_alt_fraction": 0.9}
    config2 = tmp_path / "c1" / "config2.json"
    config2.write_text(json.dumps(data), encoding="utf-8")
    _call(config2, outdir)
    diff = compare_versions(RunStore(outdir), "v0001", "v0002")
    assert diff["attribution"]["cause"] == "rule_change"
    assert diff["attribution"]["input_changed"] is False
    assert diff["attribution"]["rule_changed"] is True
    entry = diff["changed"][0]
    assert entry["changes"]["status"]["from"] == "PASS"
    assert entry["changes"]["status"]["to"] == "REJECTED"


def test_compare_added_and_removed(tmp_path):
    outdir = tmp_path / "out"
    _call(_make_case(tmp_path, alt="A", name="c1"), outdir)
    _call(_make_case(tmp_path, alt="T", name="c2"), outdir)
    diff = compare_versions(RunStore(outdir), "v0001", "v0002")
    assert diff["summary"]["added"] == 1
    assert diff["summary"]["removed"] == 1
    assert diff["added"][0]["variant_id"].endswith("G>T")
    assert diff["removed"][0]["variant_id"].endswith("G>A")


def test_export_after_restart_and_verify(tmp_path):
    outdir = tmp_path / "out"
    _call(_make_case(tmp_path, alt="A", name="c1"), outdir)
    _call(_make_case(tmp_path, alt="T", name="c2"), outdir)
    # Fresh store instance simulates a process restart.
    store = RunStore(outdir)
    dest = tmp_path / "exported"
    store.export("v0001", dest)
    for name in OUTPUT_FILENAMES:
        assert (dest / name).read_bytes() == \
            (outdir / "runs" / "v0001" / name).read_bytes()
    store.verify_all()


def test_tampered_version_fails_verify_and_export(tmp_path):
    outdir = tmp_path / "out"
    _call(_make_case(tmp_path, alt="A", name="c1"), outdir)
    target = outdir / "runs" / "v0001" / "candidates.json"
    target.write_text('{"forged": true}\n', encoding="utf-8")
    store = RunStore(outdir)
    with pytest.raises(VerificationError):
        store.verify("v0001")
    with pytest.raises(VerificationError):
        store.export("v0001", tmp_path / "dest")

