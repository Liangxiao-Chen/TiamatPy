from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import struct

from .io_json import save_json_project
from .model import Base, TiamatProject, complement_for


PENTOSE_RING_ATOMS = {"C1'", "C2'", "C3'", "C4'", "O4'"}
RESIDUE_MAP = {
    "DA": ("DNA", "A"),
    "DT": ("DNA", "T"),
    "DC": ("DNA", "C"),
    "DG": ("DNA", "G"),
    "DU": ("DNA", "U"),
    "A": ("RNA", "A"),
    "U": ("RNA", "U"),
    "C": ("RNA", "C"),
    "G": ("RNA", "G"),
}
DNA_RECORD_TYPE_CODES = {
    "A": 0,
    "T": 1,
    "U": 2,
    "C": 3,
    "G": 4,
}
ANGSTROM_TO_NM = 0.1


def load_pdb_sugar_center_project(path: str | Path) -> TiamatProject:
    file_path = Path(path)
    residues_by_chain: dict[str, list[tuple[int, str, str, tuple[float, float, float]]]] = defaultdict(list)
    residue_atoms: dict[tuple[str, int, str], dict[str, tuple[float, float, float]]] = {}
    residue_names: dict[tuple[str, int, str], str] = {}

    for line in file_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        if atom_name not in PENTOSE_RING_ATOMS:
            continue
        residue_name = line[17:20].strip().upper()
        if residue_name not in RESIDUE_MAP:
            continue
        chain_id = line[21].strip() or "_"
        residue_number = int(line[22:26].strip())
        insertion_code = line[26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        key = (chain_id, residue_number, insertion_code)
        residue_atoms.setdefault(key, {})[atom_name] = (x, y, z)
        residue_names[key] = residue_name

    for key, atoms in residue_atoms.items():
        if not PENTOSE_RING_ATOMS.issubset(atoms):
            missing = sorted(PENTOSE_RING_ATOMS - set(atoms))
            raise ValueError(f"Residue {key[0]} {key[1]} is missing pentose atoms: {', '.join(missing)}")
        residue_name = residue_names[key]
        molecule, nucleotide = RESIDUE_MAP[residue_name]
        points = [atoms[atom_name] for atom_name in sorted(PENTOSE_RING_ATOMS)]
        center = tuple(
            (sum(point[index] for point in points) / len(points)) * ANGSTROM_TO_NM
            for index in range(3)
        )
        residues_by_chain[key[0]].append((key[1], key[2], residue_name, center))

    bases: dict[int, Base] = {}
    chain_indices: dict[str, list[int]] = {}
    next_index = 1
    for chain_id in sorted(residues_by_chain):
        ordered_residues = sorted(residues_by_chain[chain_id], key=lambda item: (item[0], item[1]))
        indices: list[int] = []
        for residue_number, insertion_code, residue_name, center in ordered_residues:
            molecule, nucleotide = RESIDUE_MAP[residue_name]
            object_id = f"{chain_id}:{residue_number}{insertion_code}".rstrip()
            bases[next_index] = Base(
                index=next_index,
                object_id=object_id,
                x=center[0],
                y=center[1],
                z=center[2],
                molecule=molecule,
                nucleotide=nucleotide,
            )
            indices.append(next_index)
            next_index += 1
        for left, right in zip(indices, indices[1:]):
            bases[left].down = right
            bases[right].up = left
        chain_indices[chain_id] = indices

    project = TiamatProject(
        bases=bases,
        metadata={
            "format": "pdb_sugar_centers",
            "source_path": str(file_path),
            "source_format": "pdb",
            "coordinate_unit": "nm",
            "pdb_input_unit": "angstrom",
        },
    )
    _pair_reverse_complement_chains(project, chain_indices)
    project.refresh_strands()
    return project


def save_minimal_dna_project(project: TiamatProject, path: str | Path, include_sidecar: bool = True) -> tuple[Path, Path | None]:
    file_path = Path(path)
    record_order = _dna_record_order(project)
    payload = bytearray()
    for index in record_order:
        base = project.bases[index]
        payload.extend(_pack_dna_record(base))
    payload.extend(b"\x00" * 56)
    file_path.write_bytes(bytes(payload))

    sidecar_path: Path | None = None
    if include_sidecar:
        sidecar_path = file_path.with_suffix(".dnajson")
        sidecar_project = TiamatProject(
            bases={
                index: Base.from_dict(base.to_dict())
                for index, base in project.bases.items()
            },
            metadata=dict(project.metadata),
        )
        sidecar_project.metadata["source_path"] = str(file_path)
        sidecar_project.metadata["source_format"] = "dna"
        save_json_project(sidecar_project, sidecar_path)
    return file_path, sidecar_path


def convert_pdb_to_dna_files(
    pdb_path: str | Path,
    dna_path: str | Path | None = None,
    include_sidecar: bool = True,
) -> tuple[TiamatProject, Path, Path | None]:
    pdb_file = Path(pdb_path)
    output_path = Path(dna_path) if dna_path is not None else pdb_file.with_suffix(".dna")
    project = load_pdb_sugar_center_project(pdb_file)
    dna_file, sidecar_file = save_minimal_dna_project(project, output_path, include_sidecar=include_sidecar)
    return project, dna_file, sidecar_file


def _pair_reverse_complement_chains(project: TiamatProject, chain_indices: dict[str, list[int]]) -> None:
    if len(chain_indices) != 2:
        return
    first_chain, second_chain = sorted(chain_indices)
    first_indices = chain_indices[first_chain]
    second_indices = chain_indices[second_chain]
    if len(first_indices) != len(second_indices):
        return
    for left_index, right_index in zip(first_indices, reversed(second_indices)):
        left_base = project.bases[left_index]
        right_base = project.bases[right_index]
        if left_base.nucleotide is None or right_base.nucleotide is None:
            return
        if complement_for(left_base.nucleotide, right_base.molecule) != right_base.nucleotide:
            return
        if complement_for(right_base.nucleotide, left_base.molecule) != left_base.nucleotide:
            return
    for left_index, right_index in zip(first_indices, reversed(second_indices)):
        project.bases[left_index].across = right_index
        project.bases[right_index].across = left_index


def _dna_record_order(project: TiamatProject) -> list[int]:
    if len(project.strands) == 2:
        first = project.strands[0].base_indices
        second = project.strands[1].base_indices
        if len(first) == len(second) and all(project.bases[first[index]].across == second[-index - 1] for index in range(len(first))):
            order: list[int] = []
            for left, right in zip(first, reversed(second)):
                order.extend((left, right))
            return order
    return [index for strand in project.strands for index in strand.base_indices]


def _pack_dna_record(base: Base) -> bytes:
    type_code = 5 if base.nucleotide is None else DNA_RECORD_TYPE_CODES[base.nucleotide]
    return (
        struct.pack("<dddII", base.x, base.y, base.z, type_code, 0)
        + (b"\x00" * 20)
        + struct.pack("<I", 1)
    )
