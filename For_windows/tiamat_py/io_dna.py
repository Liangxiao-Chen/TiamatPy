from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct

from .dna_archive import decode_archive, TYPE_CODES
from .model import Base, TiamatProject


@dataclass(frozen=True)
class DnaRecord:
    offset: int
    x: float
    y: float
    z: float
    nucleotide: str | None
    molecule: str | None
    flag: int


def _read(path):
    data = Path(path).read_bytes()
    if data.startswith(b"\xff\xff"):
        project, objects = decode_archive(data)
        records = [DnaRecord(obj["offset"], *obj["position"], obj["nucleotide"], obj["molecule"], obj["flag"])
                   for obj in objects]
        return project, records

    # Old TiamatPy PDB exports contained only positions, not topology.
    # Never invent connections or borrow a nearby JSON file for these records.
    if len(data) < 56 or len(data) % 56:
        raise ValueError("Unsupported or truncated .dna archive; no valid Nucleobase records.")
    end = len(data) - 56 if data.endswith(b"\x00" * 56) else len(data)
    records = []
    for offset in range(0, end, 56):
        payload = data[offset:offset + 56]
        x, y, z, code, flag = struct.unpack_from("<dddII", payload)
        if code not in TYPE_CODES or flag not in (0, 1) or payload[32:52] != b"\x00" * 20 or payload[52:] != b"\x01\x00\x00\x00":
            raise ValueError(f"Unsupported .dna record at byte {offset}.")
        records.append(DnaRecord(offset, x, y, z, TYPE_CODES[code], "RNA" if code == 2 else "DNA", flag))
    if not records:
        raise ValueError("No Nucleobase records in .dna file.")
    project = TiamatProject({i: Base(i, str(i), r.x, r.y, r.z, molecule=r.molecule, nucleotide=r.nucleotide)
                            for i, r in enumerate(records)}, metadata={
        "dna_topology": "unavailable",
        "import_warnings": ["This legacy position-only file contains no strand or pair references. "
                            "Bases were imported without bonds. Open its .dnajson instead to recover topology."]})
    return project, records


def load_dna_project(path: str | Path) -> TiamatProject:
    project, records = _read(path)
    project.metadata.update(source_path=str(Path(path).resolve()), source_format="dna", dna_record_count=len(records))
    return project


def parse_dna_records(path: str | Path) -> list[DnaRecord]:
    return _read(path)[1]
