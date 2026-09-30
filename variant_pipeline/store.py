"""Versioned run store: every analysis publishes a new immutable version.

Each ``call`` writes its outputs into a fresh ``runs/vNNNN/`` directory
under the store root (the ``--outdir``). Existing versions are never
overwritten or recomputed in place: publishing goes through a temporary
staging directory plus an atomic rename, so a failed publish leaves all
previously published versions untouched. A repeated request with identical
inputs and interpretation rules returns the already published version
instead of re-analysing.

The manifest of each version records the input fingerprints (config,
reads, reference FASTA) and the analysis rule version together with the
output file hashes, so input / rules / output stay bound together and any
later modification of a historical version is detectable via ``verify``.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import ANALYSIS_VERSION, __version__

RUNS_DIRNAME = "runs"
INDEX_FILENAME = "index.json"
MANIFEST_FILENAME = "manifest.json"
VERSION_PATTERN = re.compile(r"^v(\d{4,})$")

OUTPUT_FILENAMES = (
    "candidates.json",
    "candidates.tsv",
    "coverage.tsv",
    "report.txt",
)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class VersionInfo:
    version: str
    path: Path
    manifest: dict


class StoreError(RuntimeError):
    pass


class VerificationError(StoreError):
    pass


class RunStore:
    """File-backed store of immutable analysis versions."""

    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.runs_dir = self.root / RUNS_DIRNAME

    # -- publishing ------------------------------------------------------

    def fingerprint(self, *, config_dict: dict, reads_path: str,
                    fasta_path: str) -> str:
        """Stable identity of one analysis request.

        Binds the exact inputs (reads + reference bytes), the full run
        configuration and the interpretation rule version, so a repeated
        request can be recognised without re-running the analysis.
        """
        material = json.dumps(
            {
                "analysis_version": ANALYSIS_VERSION,
                "package_version": __version__,
                "config": config_dict,
                "reads_sha256": sha256_file(reads_path),
                "fasta_sha256": sha256_file(fasta_path),
            },
            sort_keys=True,
            ensure_ascii=False,
        )
        return sha256_text(material)

    def find_duplicate(self, fingerprint: str) -> VersionInfo | None:
        for info in self.list_versions():
            if info.manifest.get("fingerprint") == fingerprint:
                return info
        return None

    def publish(self, outputs: dict[str, str | Path], *,
                config_dict: dict, reads_path: str, fasta_path: str,
                fingerprint: str) -> VersionInfo:
        """Publish a new version atomically; never touches old versions.

        ``outputs`` maps output filenames to the staged file paths to
        import. Files are copied into a temporary sibling directory, the
        manifest (with output hashes) is written there, and the directory
        is atomically renamed to ``runs/vNNNN``. If anything fails before
        the rename, previously published versions are unaffected and the
        partial staging directory is removed.
        """
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(
            prefix=".staging-", dir=str(self.runs_dir)))
        try:
            output_hashes = {}
            for name in OUTPUT_FILENAMES:
                src = Path(outputs[name])
                dst = staging / name
                shutil.copyfile(src, dst)
                output_hashes[name] = sha256_file(dst)
            manifest = {
                "version": None,  # filled in below, before rename
                "created_at": datetime.now(timezone.utc).isoformat(),
                "analysis_version": ANALYSIS_VERSION,
                "package_version": __version__,
                "fingerprint": fingerprint,
                "inputs": {
                    "config": config_dict,
                    "config_sha256": sha256_text(json.dumps(
                        config_dict, sort_keys=True, ensure_ascii=False)),
                    "reads_path": str(reads_path),
                    "reads_sha256": sha256_file(reads_path),
                    "fasta_path": str(fasta_path),
                    "fasta_sha256": sha256_file(fasta_path),
                },
                "outputs": output_hashes,
            }
            version = self._next_version()
            manifest["version"] = version
            (staging / MANIFEST_FILENAME).write_text(
                json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            final_dir = self.runs_dir / version
            # Atomic on the same filesystem; fails rather than overwriting
            # if the target somehow already exists.
            if final_dir.exists():
                raise StoreError(
                    f"version directory already exists: {final_dir}")
            os.rename(staging, final_dir)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        self._register(version)
        return VersionInfo(version=version, path=final_dir, manifest=manifest)

    # -- reading ---------------------------------------------------------

    def list_versions(self) -> list[VersionInfo]:
        infos = []
        if not self.runs_dir.is_dir():
            return infos
        for entry in sorted(self.runs_dir.iterdir()):
            if entry.is_dir() and VERSION_PATTERN.match(entry.name):
                infos.append(self.get(entry.name))
        return infos

    def get(self, version: str) -> VersionInfo:
        path = self.runs_dir / version
        manifest_path = path / MANIFEST_FILENAME
        if not manifest_path.is_file():
            raise StoreError(f"unknown version: {version}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return VersionInfo(version=version, path=path, manifest=manifest)

    def verify(self, version: str) -> VersionInfo:
        """Recompute output hashes; raise if the version was rewritten."""
        info = self.get(version)
        for name, recorded in info.manifest.get("outputs", {}).items():
            path = info.path / name
            if not path.is_file():
                raise VerificationError(
                    f"{version}: missing output file {name}")
            actual = sha256_file(path)
            if actual != recorded:
                raise VerificationError(
                    f"{version}: {name} hash mismatch "
                    f"(recorded {recorded[:12]}..., "
                    f"actual {actual[:12]}...); historical content "
                    "appears to have been rewritten")
        return info

    def verify_all(self) -> list[VersionInfo]:
        return [self.verify(info.version) for info in self.list_versions()]

    def export(self, version: str, dest: str | Path) -> Path:
        """Export a historical version after verifying its integrity."""
        info = self.verify(version)
        dest = Path(dest)
        if dest.exists() and any(dest.iterdir()):
            raise StoreError(f"export destination is not empty: {dest}")
        dest.mkdir(parents=True, exist_ok=True)
        for name in list(OUTPUT_FILENAMES) + [MANIFEST_FILENAME]:
            shutil.copyfile(info.path / name, dest / name)
        return dest

    # -- internals -------------------------------------------------------

    def _next_version(self) -> str:
        highest = 0
        if self.runs_dir.is_dir():
            for entry in self.runs_dir.iterdir():
                match = VERSION_PATTERN.match(entry.name)
                if match:
                    highest = max(highest, int(match.group(1)))
        return f"v{highest + 1:04d}"

    def _register(self, version: str) -> None:
        index_path = self.root / INDEX_FILENAME
        entries = []
        if index_path.is_file():
            entries = json.loads(index_path.read_text(encoding="utf-8"))
        if version not in entries:
            entries.append(version)
        fd, tmp = tempfile.mkstemp(prefix=".index-", dir=str(self.root))
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(entries, indent=2) + "\n")
        os.replace(tmp, index_path)

