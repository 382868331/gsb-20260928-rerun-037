import json

from variant_pipeline.cli import main
from variant_pipeline.store import OUTPUT_FILENAMES, RunStore
from tests.helpers import RNAME, write_case


def test_export_contains_ledger_and_version_decision(tmp_path):
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    fasta, sam = write_case(case_dir, [])
    config = case_dir / "config.json"
    config.write_text(json.dumps({
        "sample": {"sample_id": "S"},
        "reference": {
            "fasta_path": str(fasta),
            "name": RNAME,
            "assembly": "x",
        },
        "reads_path": str(sam),
    }), encoding="utf-8")
    out = tmp_path / "out"
    assert main(["call", "--config", str(config), "--outdir", str(out)]) == 0
    exported = tmp_path / "exported"
    assert main(["export", "--outdir", str(out), "--version", "v0001",
                 "--dest", str(exported)]) == 0
    for name in OUTPUT_FILENAMES:
        assert (exported / name).is_file()
    manifest = json.loads((exported / "manifest.json").read_text("utf-8"))
    assert "evidence_ledger.json" in manifest["outputs"]
    ledger = json.loads((exported / "evidence_ledger.json").read_text("utf-8"))
    assert ledger["schema"] == "fragment-evidence-ledger-v1"
