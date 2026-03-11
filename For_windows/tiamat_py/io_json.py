from __future__ import annotations

import json
from pathlib import Path

from .model import Base, TiamatProject, _normalize_base


def load_json_project(path: str | Path) -> TiamatProject:
    file_path = Path(path)
    data = json.loads(file_path.read_text(encoding="utf-8"))
    project = _parse_json_project(data)
    project.metadata.setdefault("source_path", str(file_path))
    project.metadata.setdefault("source_format", "dnajson")
    return project


def save_json_project(project: TiamatProject, path: str | Path) -> None:
    file_path = Path(path)
    file_path.write_text(project.to_json() + "\n", encoding="utf-8")


def _parse_json_project(data: dict) -> TiamatProject:
    bases = data.get("bases")
    if not isinstance(bases, list):
        return TiamatProject.from_dict(data)
    if bases and isinstance(bases[0], dict) and "position" in bases[0]:
        parsed_bases = {}
        for item in bases:
            position = item.get("position") or [0.0, 0.0, 0.0]
            index = int(item["id"])
            parsed_bases[index] = Base(
                index=index,
                object_id=str(item.get("id", index)),
                x=float(position[0]),
                y=float(position[1]),
                z=float(position[2]),
                molecule=str(item.get("molecule", "DNA")).upper(),
                up=_as_int_or_none(item.get("up")),
                down=_as_int_or_none(item.get("down")),
                across=_as_int_or_none(item.get("across")),
                nucleotide=_normalize_base(item.get("type")),
            )
        return TiamatProject(
            bases=parsed_bases,
            metadata={"format": "external_dnajson"},
        )
    return TiamatProject.from_dict(data)


def _as_int_or_none(value):
    if value is None:
        return None
    return int(value)
