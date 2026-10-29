from __future__ import annotations

from collections import defaultdict
import math
from pathlib import Path

from .dna_archive import encode_archive
from .geometry import cross, dot, mean, norm, perpendicular, sub, unit
from .io_json import save_json_project
from .model import Base, TiamatProject, bases_can_pair

PENTOSE_RING_ATOMS = {"C1'", "C2'", "C3'", "C4'", "O4'"}
RESIDUE_MAP = {"DA": ("DNA", "A"), "DT": ("DNA", "T"), "DC": ("DNA", "C"),
               "DG": ("DNA", "G"), "DU": ("DNA", "U"),
               "A": ("RNA", "A"), "U": ("RNA", "U"), "C": ("RNA", "C"), "G": ("RNA", "G")}
ANGSTROM_TO_NM = 0.1
RING_ATOMS = ("N1", "C2", "N3", "C4", "C5", "C6")


def load_pdb_sugar_center_project(path: str | Path) -> TiamatProject:
    file_path = Path(path)
    residues = {}
    segment = 0
    model_seen = False
    for line in file_path.read_text(encoding="utf-8", errors="strict").splitlines():
        record = line[:6].strip()
        if record == "MODEL":
            if model_seen:
                break
            model_seen = True
            continue
        if record == "ENDMDL":
            break
        if record == "TER":
            segment += 1
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        residue_name = line[17:20].strip().upper()
        if residue_name not in RESIDUE_MAP:
            continue
        chain = line[21].strip() or "_"
        key = (segment, chain, int(line[22:26]), line[26].strip())
        residue = residues.setdefault(key, {"name": residue_name, "atoms": defaultdict(dict)})
        atom = line[12:16].strip().replace("*", "'")
        alternate = line[16].strip()
        position = tuple(float(line[start:start+8]) * ANGSTROM_TO_NM for start in (30, 38, 46))
        occupancy = float(line[54:60].strip() or "1")
        residue["atoms"][alternate][atom] = (position, occupancy)

    bases = {}
    atom_sets = {}
    previous = {}
    notes = []
    for key, residue in residues.items():
        conformers = residue["atoms"]
        alternatives = [alt for alt in conformers if alt]
        selected = min(alternatives, key=lambda alt: (-sum(item[1] for item in conformers[alt].values()) / len(conformers[alt]), alt)) if alternatives else ""
        # One conformer per residue; blank atoms are shared, not a competing model.
        atoms = {name: value[0] for name, value in conformers.get("", {}).items()}
        atoms.update({name: value[0] for name, value in conformers.get(selected, {}).items()})
        chain_key = key[:2]
        if not PENTOSE_RING_ATOMS.issubset(atoms):
            notes.append(f"Skipped residue {key[1]}:{key[2]}{key[3]}: incomplete pentose ring.")
            previous.pop(chain_key, None)
            continue
        molecule, nucleotide = RESIDUE_MAP[residue["name"]]
        index = len(bases) + 1
        center = mean([atoms[name] for name in sorted(PENTOSE_RING_ATOMS)])
        base = Base(index, f"{key[1]}:{key[2]}{key[3]}:segment{key[0]}", *center, molecule=molecule, nucleotide=nucleotide)
        if all(name in atoms for name in RING_ATOMS):
            try:
                normal = unit(cross(sub(atoms["C4"], atoms["C2"]), sub(atoms["C6"], atoms["C2"])))
                base.plane_u = unit(perpendicular(sub(mean([atoms[n] for n in RING_ATOMS]), center), normal))
                base.plane_v = cross(normal, base.plane_u)
            except ValueError:
                pass
        bases[index] = base
        atom_sets[index] = atoms
        if chain_key in previous:
            prev_index, prev_key = previous[chain_key]
            old_atoms = atom_sets[prev_index]
            if "O3'" in old_atoms and "P" in atoms:
                connected = 0.10 <= math.dist(old_atoms["O3'"], atoms["P"]) <= 0.23
            else:
                sequential = key[2] == prev_key[2] + 1 or (key[2] == prev_key[2] and key[3] > prev_key[3])
                connected = sequential and 0.2 <= math.dist(bases[prev_index].position, center) <= 1.0
                if connected:
                    notes.append("Some backbone links were inferred from residue order and sugar spacing because O3'/P atoms are missing.")
            if connected:
                bases[prev_index].down = index
                base.up = prev_index
        previous[chain_key] = (index, key)
    if not bases:
        raise ValueError("No supported nucleic acid residues with complete pentose rings in this PDB model.")
    project = TiamatProject(bases, metadata={"format": "pdb_sugar_centers", "source_path": str(file_path.resolve()),
        "source_format": "pdb", "coordinate_unit": "nm", "pdb_input_unit": "angstrom", "pdb_model": "first",
        "pairing_method": "base_planes_and_heavy_atom_contacts"})
    _infer_pairs(project, atom_sets)
    if notes:
        project.metadata["import_warnings"] = list(dict.fromkeys(notes))[:12]
    project.validate(check_pairs=True)
    return project


def _infer_pairs(project, atom_sets):
    """Conservative geometric candidates, not sequence-alignment pairing.

    Require coplanar base rings and at least two expected donor/acceptor heavy-atom
    contacts (0.22-0.36 nm). Ambiguous competing partners are left unassigned.
    """
    geometry = {}
    cells = defaultdict(list)
    for index, atoms in atom_sets.items():
        base = project.bases[index]
        if base.plane_u is None or base.plane_v is None:
            continue
        center = mean([atoms[n] for n in RING_ATOMS])
        normal = unit(cross(base.plane_u, base.plane_v))
        geometry[index] = (center, normal)
        cells[tuple(math.floor(v / 1.2) for v in center)].append(index)
    candidates = defaultdict(list)
    contacts = {("A", "T"): (("N6", "O4"), ("N1", "N3")),
                ("A", "U"): (("N6", "O4"), ("N1", "N3")),
                ("C", "G"): (("N4", "O6"), ("N3", "N1"), ("O2", "N2")),
                ("G", "U"): (("O6", "N3"), ("N1", "O2"))}
    for index, (center, normal) in geometry.items():
        cell = tuple(math.floor(v / 1.2) for v in center)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in cells.get((cell[0]+dx, cell[1]+dy, cell[2]+dz), ()):
                        if other <= index:
                            continue
                        left, right = project.bases[index], project.bases[other]
                        if other in (left.up, left.down) or not bases_can_pair(left, right):
                            continue
                        other_center, other_normal = geometry[other]
                        separation = sub(other_center, center)
                        if abs(dot(normal, other_normal)) < 0.85 or abs(dot(separation, normal)) > 0.18 or abs(dot(separation, other_normal)) > 0.18:
                            continue
                        expected = contacts.get((left.nucleotide, right.nucleotide))
                        if expected is None:
                            expected = tuple((b, a) for a, b in contacts.get((right.nucleotide, left.nucleotide), ()))
                        distances = [math.dist(atom_sets[index][a], atom_sets[other][b]) for a, b in expected
                                     if a in atom_sets[index] and b in atom_sets[other]]
                        valid = [d for d in distances if 0.22 <= d <= 0.36]
                        if len(valid) < 2:
                            continue
                        score = sum(abs(d - 0.29) for d in valid) / len(valid)
                        candidates[index].append((score, other))
                        candidates[other].append((score, index))
    best = {}
    for index, options in candidates.items():
        options.sort()
        if len(options) == 1 or options[1][0] - options[0][0] > 0.02:
            best[index] = options[0][1]
    for index, other in best.items():
        if best.get(other) == index:
            project.bases[index].across = other


def save_minimal_dna_project(project: TiamatProject, path: str | Path, include_sidecar: bool = True) -> tuple[Path, Path | None]:
    # Retain the old API name; the archive now contains complete object references.
    file_path = Path(path)
    file_path.write_bytes(encode_archive(project))
    sidecar = file_path.with_suffix(".dnajson") if include_sidecar else None
    if sidecar is not None:
        save_json_project(project, sidecar)
    return file_path, sidecar


def convert_pdb_to_dna_files(pdb_path: str | Path, dna_path: str | Path | None = None,
                             include_sidecar: bool = True) -> tuple[TiamatProject, Path, Path | None]:
    project = load_pdb_sugar_center_project(pdb_path)
    path, sidecar = save_minimal_dna_project(project, dna_path or Path(pdb_path).with_suffix(".dna"), include_sidecar)
    return project, path, sidecar
