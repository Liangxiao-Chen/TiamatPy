from __future__ import annotations

import json
from pathlib import Path

from .model import Base, TiamatProject, _normalize_base, _as_int_or_none
from .persistence import atomic_write_text


def load_json_project(path: str | Path) -> TiamatProject:
    file_path = Path(path)
    data = json.loads(file_path.read_text(encoding="utf-8"))
    project = _parse_json_project(data)
    project.metadata["source_path"] = str(file_path.resolve())
    project.metadata["source_format"] = "dnajson"
    return project


def save_json_project(project: TiamatProject, path: str | Path) -> None:
    project.validate()
    data = project.to_dict()
    data["metadata"] = dict(project.metadata, source_path=str(Path(path).resolve()), source_format="dnajson")
    atomic_write_text(path, json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _parse_json_project(data: dict) -> TiamatProject:
    bases = data.get("bases")
    if not isinstance(bases, list):
        return TiamatProject.from_dict(data)
    if bases and isinstance(bases[0], dict) and "position" in bases[0]:
        parsed_bases = {}
        for item in bases:
            position = item["position"]
            if not isinstance(position, (list, tuple)) or len(position) != 3:
                raise ValueError("Each base must have exactly three position coordinates.")
            index = _as_int_or_none(item["id"])
            if index is None or index in parsed_bases:
                raise ValueError(f"Missing or duplicate base ID: {index!r}.")
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
