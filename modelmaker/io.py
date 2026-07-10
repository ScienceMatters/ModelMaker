from __future__ import annotations

import csv
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path

from .models import MutationSpec, TargetSpec


def _open_csv(path: Path):
    return path.open(newline="", encoding="utf-8-sig")


def read_mutations_csv(path: Path) -> list[MutationSpec]:
    out: list[MutationSpec] = []
    with _open_csv(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out.append(
                MutationSpec(
                    gene=row["gene"].strip().upper(),
                    mutation=row["mutation"].strip().upper(),
                    position=int(row["position"]),
                    numbering_offset=int(row.get("numbering_offset", 0) or 0),
                )
            )
    return out


def read_targets_csv(path: Path) -> list[TargetSpec]:
    out: list[TargetSpec] = []
    with _open_csv(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out.append(
                TargetSpec(
                    gene=row["gene"].strip().upper(),
                    human=row["human"].strip(),
                    zfish=row["zfish"].strip(),
                )
            )
    return out


def write_json(path: Path, obj) -> None:
    if is_dataclass(obj):
        obj = asdict(obj)
    elif hasattr(obj, "__dict__") and not isinstance(obj, dict):
        obj = obj.__dict__
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
