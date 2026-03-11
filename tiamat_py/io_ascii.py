from __future__ import annotations

from pathlib import Path

from .model import Base, TiamatProject


def load_ascii_dat(path: str | Path) -> TiamatProject:
    file_path = Path(path)
    text = file_path.read_text(encoding="utf-8", errors="ignore")
    project = parse_ascii_dat(text)
    project.metadata["source_path"] = str(file_path)
    project.metadata["source_format"] = "ascii_dat"
    return project


def parse_ascii_dat(text: str) -> TiamatProject:
    rows: list[list[str]] = []
    expected_count: int | None = None
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) == 1 and parts[0].isdigit():
            expected_count = int(parts[0])
            continue
        if len(parts) == 8 and parts[0].isdigit():
            rows.append(parts)
    if not rows:
        raise ValueError("no Tiamat ASCII rows found")
    if expected_count is not None and expected_count != len(rows):
        raise ValueError(f"expected {expected_count} bases but found {len(rows)} rows")

    object_id_to_index = {parts[1]: int(parts[0]) for parts in rows}
    bases = {}
    for parts in rows:
        index = int(parts[0])
        bases[index] = Base(
            index=index,
            object_id=parts[1],
            up=_resolve(parts[2], object_id_to_index),
            down=_resolve(parts[3], object_id_to_index),
            across=_resolve(parts[4], object_id_to_index),
            x=float(parts[5]),
            y=float(parts[6]),
            z=float(parts[7]),
        )
    return TiamatProject(
        bases=bases,
        metadata={
            "source_format": "ascii_dat",
            "base_count": len(bases),
        },
    )


def _resolve(token: str, object_id_to_index: dict[str, int]) -> int | None:
    if token == "00000000":
        return None
    return object_id_to_index[token]

