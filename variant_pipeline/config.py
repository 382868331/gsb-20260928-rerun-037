from __future__ import annotations

import json
from pathlib import Path

from .models import Parameters, SampleConfig


def load_manifest(path: str | Path) -> SampleConfig:
    manifest_path = Path(path)
    with manifest_path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)

    try:
        sample_id = raw["sample_id"]
        reference = raw["reference"]
        reads = raw["reads"]
    except KeyError as exc:
        raise ValueError(f"Manifest is missing required field: {exc.args[0]}") from exc

    parameter_values = raw.get("parameters", {})
    allowed = Parameters.__dataclass_fields__
    unknown = sorted(set(parameter_values) - set(allowed))
    if unknown:
        raise ValueError(f"Unknown analysis parameter(s): {', '.join(unknown)}")

    base_dir = manifest_path.parent
    config = SampleConfig(
        sample_id=sample_id,
        reference_name=reference["name"],
        reference_version=reference["version"],
        reference_fasta=str(_resolve(base_dir, reference["fasta"])),
        reads_sam=str(_resolve(base_dir, reads["sam"])),
        parameters=Parameters(**parameter_values),
    )
    return config


def _resolve(base_dir: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base_dir / path
