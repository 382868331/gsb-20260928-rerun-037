"""Versioned, immutable batch storage.

Every analysis run becomes its own batch directory under the output root:
``batches/<timestamp>_<short-content-id>``. A batch keeps the resolved
analysis artifacts together with a snapshot of every input that produced
them (config JSON, FASTA and SAM bytes). Nothing in an existing batch is
ever rewritten: a repeat request with the same inputs creates a new batch
whose identity/duplicate status is recorded in ``index.json``.

Publication is atomic: artifacts are written to a hidden staging directory,
the manifest is stamped last, and the staging directory is renamed into
place. A failed or interrupted run therefore leaves no partial batch and
cannot damage earlier versions. The batch index is replaced atomically as
well.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import ANALYSIS_VERSION, __version__
from .config import RunConfig

INDEX_NAME = "index.json"
MANIFEST_NAME = "manifest.json"
BATCHES_DIRNAME = "batches"
ARTIFACT_NAMES = (
    "candidates.json",
    "candidates.tsv",
    "coverage.tsv",
    "report.txt",
)
INPUT_NAMES = ("config.json", "reference.fasta", "reads.sam")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


@dataclass(frozen=True)
class BatchInputSet:
    config_bytes: bytes
    reference_bytes: bytes
    reads_bytes: bytes
    config_path: str
    reference_path: str
    reads_path: str

    def content_id(self) -> str:
        h = hashlib.sha256()
        h.update(b"variant-pipeline-batch-v1\n")
        h.update(b"analysis_version=")
        h.update(ANALYSIS_VERSION.encode("utf-8"))
        h.update(b"\nconfig:")
        h.update(self.config_bytes)
        h.update(b"\nreference:")
        h.update(self.reference_bytes)
        h.update(b"\nreads:")
        h.update(self.reads_bytes)
        return h.hexdigest()

    @property
    def input_hashes(self) -> dict[str, str]:
        return {
            "config.json": sha256_bytes(self.config_bytes),
            "reference.fasta": sha256_bytes(self.reference_bytes),
            "reads.sam": sha256_bytes(self.reads_bytes),
        }

    @property
    def input_sources(self) -> dict[str, str]:
        return {
            "config.json": self.config_path,
            "reference.fasta": self.reference_path,
            "reads.sam": self.reads_path,
        }


def load_input_set(
    cfg: RunConfig,
    config_path: str | Path,
    reads_override: str | None = None,
) -> BatchInputSet:
    config_path = Path(config_path)
    reads_path = Path(reads_override) if reads_override else Path(cfg.reads_path)
    if cfg.reads_path is None and reads_override is None:
        raise ValueError("no reads path: set reads_path in config or pass --reads")
    reference_path = Path(cfg.reference.fasta_path)
    return BatchInputSet(
        config_bytes=config_path.read_bytes(),
        reference_bytes=reference_path.read_bytes(),
        reads_bytes=reads_path.read_bytes(),
        config_path=str(config_path),
        reference_path=str(reference_path),
        reads_path=str(reads_path),
    )


def load_index(root: str | Path) -> dict:
    path = Path(root) / INDEX_NAME
    if not path.exists():
        return {"schema": "variant-pipeline-index/1", "batches": []}
    return json.loads(path.read_text(encoding="utf-8"))


def write_index_atomic(root: str | Path, index: dict) -> None:
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    final = root / INDEX_NAME
    tmp = root / f".{INDEX_NAME}.{os.getpid()}.tmp"
    tmp.write_text(
        json.dumps(index, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    os.replace(tmp, final)


def _same_content_id(index: dict, content_id: str) -> list[str]:
    return [
        b["batch_id"]
        for b in index.get("batches", [])
        if b.get("content_id") == content_id
    ]


def find_duplicate_batches(root: str | Path, content_id: str) -> list[str]:
    """Batch ids previously produced from the same inputs under the same rules."""
    return _same_content_id(load_index(root), content_id)


def _unique_batch_dir(root: Path, content_id: str) -> Path:
    for _ in range(1000):
        batch_id = f"{utc_stamp()}_{content_id[:12]}"
        candidate = root / BATCHES_DIRNAME / batch_id
        if not candidate.exists():
            return candidate
    raise RuntimeError("could not allocate a unique batch id")


def publish_batch(
    root: str | Path,
    cfg: RunConfig,
    inputs: BatchInputSet,
    payload: dict,
    artifacts: dict[str, bytes],
) -> dict:
    """Stage and atomically publish one immutable batch, then update index."""
    root = Path(root)
    index = load_index(root)
    content_id = inputs.content_id()
    duplicates_of = _same_content_id(index, content_id)

    batch_dir = _unique_batch_dir(root, content_id)
    stage = batch_dir.parent / f".stage_{batch_dir.name}_{os.getpid()}"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    try:
        input_dir = stage / "input"
        input_dir.mkdir()
        input_bytes = {
            "config.json": inputs.config_bytes,
            "reference.fasta": inputs.reference_bytes,
            "reads.sam": inputs.reads_bytes,
        }
        for name, data in input_bytes.items():
            (input_dir / name).write_bytes(data)

        for name, data in artifacts.items():
            (stage / name).write_bytes(data)

        files = {}
        for name, data in input_bytes.items():
            files[f"input/{name}"] = sha256_bytes(data)
        for name, data in artifacts.items():
            files[name] = sha256_bytes(data)

        manifest = {
            "schema": "variant-pipeline-batch/1",
            "batch_id": batch_dir.name,
            "created_at_utc": utc_stamp(),
            "package_version": __version__,
            "analysis_version": ANALYSIS_VERSION,
            "content_id": content_id,
            "duplicate_of": duplicates_of,
            "sample_id": cfg.sample.sample_id,
            "reference": {
                "name": cfg.reference.name,
                "assembly": cfg.reference.assembly,
            },
            "input_sources": inputs.input_sources,
            "input_hashes": inputs.input_hashes,
            "parameters": cfg.parameters.__dict__,
            "files": files,
        }
        (stage / MANIFEST_NAME).write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        os.replace(stage, batch_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise

    entry = {
        "batch_id": batch_dir.name,
        "created_at_utc": manifest["created_at_utc"],
        "analysis_version": ANALYSIS_VERSION,
        "package_version": __version__,
        "sample_id": cfg.sample.sample_id,
        "reference_assembly": cfg.reference.assembly,
        "content_id": content_id,
        "input_hashes": inputs.input_hashes,
        "duplicate_of": duplicates_of,
        "n_candidates": len(payload.get("candidates", [])),
        "manifest": MANIFEST_NAME,
    }
    index.setdefault("batches", []).append(entry)
    write_index_atomic(root, index)
    return {"manifest": manifest, "entry": entry, "batch_dir": batch_dir}


def batch_dir(root: str | Path, batch_id: str) -> Path:
    path = Path(root) / BATCHES_DIRNAME / batch_id
    if not path.is_dir():
        raise KeyError(f"unknown batch id: {batch_id}")
    return path


def resolve_batch_id(root: str | Path, ref: str) -> str:
    """Resolve a full id, a unique id prefix, or 'latest' to a batch id."""
    index = load_index(root)
    batches = index.get("batches", [])
    if not batches:
        raise KeyError("no batches in this output root")
    if ref == "latest":
        return batches[-1]["batch_id"]
    exact = [b for b in batches if b["batch_id"] == ref]
    if exact:
        return exact[0]["batch_id"]
    prefix = [b for b in batches if b["batch_id"].startswith(ref)]
    if len(prefix) == 1:
        return prefix[0]["batch_id"]
    if len(prefix) > 1:
        raise KeyError(f"batch prefix {ref!r} is ambiguous")
    if (Path(root) / BATCHES_DIRNAME / ref).is_dir():
        return ref
    raise KeyError(f"unknown batch id: {ref}")


def load_manifest(root: str | Path, batch_id: str) -> dict:
    path = batch_dir(root, batch_id) / MANIFEST_NAME
    return json.loads(path.read_text(encoding="utf-8"))


def load_payload(root: str | Path, batch_id: str) -> dict:
    path = batch_dir(root, batch_id) / "candidates.json"
    return json.loads(path.read_text(encoding="utf-8"))


def verify_batch(root: str | Path, batch_id: str) -> dict:
    """Check stored hashes only; never re-run analysis on historical inputs."""
    bdir = batch_dir(root, batch_id)
    manifest = load_manifest(root, batch_id)
    checks = {}
    ok = True
    for rel, expected in manifest["files"].items():
        target = bdir / rel
        if not target.exists():
            checks[rel] = {"ok": False, "reason": "missing"}
            ok = False
            continue
        actual = sha256_file(target)
        good = actual == expected
        checks[rel] = {"ok": good, "sha256": actual}
        ok = ok and good
    return {"batch_id": batch_id, "ok": ok, "files": checks}


def export_batch(
    root: str | Path,
    batch_id: str,
    dest: str | Path,
    overwrite: bool = False,
) -> Path:
    """Copy the frozen batch artifacts; historical bytes are never regenerated."""
    src = batch_dir(root, batch_id)
    dest = Path(dest)
    if dest.exists():
        if not overwrite:
            raise FileExistsError(f"export destination already exists: {dest}")
        if not dest.is_dir():
            raise NotADirectoryError(f"export destination is not a directory: {dest}")
    dest.mkdir(parents=True, exist_ok=True)
    for name in (*ARTIFACT_NAMES, MANIFEST_NAME):
        shutil.copy2(src / name, dest / name)
    input_dest = dest / "input"
    input_dest.mkdir(exist_ok=True)
    for name in INPUT_NAMES:
        shutil.copy2(src / "input" / name, input_dest / name)
    return dest
