"""Backbone-preserving coarse-grained export, independently implemented.

Tiamat positions (including PDB sugar centers) are treated as backbone sites,
NOT centers of mass. This is an initial configuration, not an energy minimizer.
Site offsets follow the documented oxDNA2 and oxRNA models; see docs/validation.md.
"""
from __future__ import annotations

from collections import deque
import math
from pathlib import Path
import random

from .geometry import add, sub, scale, dot, cross, norm, unit, perpendicular
from .model import TiamatProject, complement_for
from .persistence import atomic_write_files

NM_PER_OXDNA_UNIT = 0.8518
BACKBONE_OFFSETS = {"DNA2": (-0.34, 0.3408, 0.0), "RNA": (-0.4, 0.0, 0.2)}


def resolve_oxdna_output_paths(output_path: str | Path) -> tuple[Path, Path]:
    path = Path(output_path)
    return path.with_suffix(".top"), path.with_suffix(".oxdna")


def _tangent_to_five_prime(project, base):
    before = project.bases[base.up].position if base.up is not None else base.position
    after = project.bases[base.down].position if base.down is not None else base.position
    tangent = sub(before, after)
    if norm(tangent) > 1e-10:
        return tangent
    if base.up is not None:
        return sub(before, base.position)
    return None


def _orientation_frames(project):
    frames = {}
    for index, base in project.bases.items():
        tangent = _tangent_to_five_prime(project, base)
        normal = cross(base.plane_u, base.plane_v) if base.plane_u is not None and base.plane_v is not None else None
        try:
            if base.across is not None:
                a1 = unit(sub(project.bases[base.across].position, base.position))
                if normal is None:
                    # Pair midpoints cancel the winding of the two backbones.
                    midpoint = scale(add(base.position, project.bases[base.across].position), 0.5)
                    for neighbor_id, sign in ((base.up, 1), (base.down, -1)):
                        if neighbor_id is None:
                            continue
                        neighbor = project.bases[neighbor_id]
                        if neighbor.across is not None:
                            other_midpoint = scale(add(neighbor.position, project.bases[neighbor.across].position), 0.5)
                            normal = scale(sub(other_midpoint, midpoint), sign)
                            if norm(perpendicular(normal, a1)) > 1e-10:
                                break
                    if normal is None or norm(perpendicular(normal, a1)) < 1e-10:
                        normal = tangent
                if normal is None:
                    continue
                a3 = unit(perpendicular(normal, a1))
            elif normal is not None:
                a1 = unit(base.plane_u)
                a3 = unit(perpendicular(normal, a1))
            elif tangent is not None and base.up is not None and base.down is not None:
                a3 = unit(tangent)
                bend = sub(add(project.bases[base.up].position, project.bases[base.down].position), scale(base.position, 2))
                a1 = unit(perpendicular(bend, a3))
            else:
                continue
            if tangent is not None and dot(a3, tangent) < 0:
                a3 = scale(a3, -1)
            elif tangent is None and base.across is not None and index > base.across:
                a3 = scale(a3, -1)
            frames[index] = (a1, a3)
        except ValueError:
            continue

    # Parallel transport from established frames into unpaired tails. No fake
    # world-space nucleotide is introduced for a missing neighbor.
    queue = deque(sorted(frames))
    while queue:
        source = queue.popleft()
        base = project.bases[source]
        for neighbor in (base.up, base.down, base.across):
            if neighbor is None or neighbor in frames:
                continue
            target = project.bases[neighbor]
            tangent = _tangent_to_five_prime(project, target)
            try:
                a3 = unit(tangent) if tangent is not None else frames[source][1]
                a1 = unit(perpendicular(frames[source][0], a3))
                if target.across == source:
                    a1 = unit(sub(base.position, target.position))
                    a3 = unit(perpendicular(scale(frames[source][1], -1), a1))
                frames[neighbor] = (a1, a3)
                queue.append(neighbor)
            except ValueError:
                continue
    missing = sorted(set(project.bases) - set(frames))
    if missing:
        raise ValueError(f"Cannot determine base orientation for IDs {missing[:8]}. "
                         "An isolated point or straight unpaired strand has no defined base plane. "
                         "Use a helix/PDB with base planes or supply plane_u/plane_v in the project; "
                         "export will not invent world-axis orientations.")
    return frames


def build_oxdna_configuration(project: TiamatProject, *, model: str | None = None):
    project.validate(check_pairs=True)
    molecules = {base.molecule for base in project.bases.values()}
    if not molecules:
        raise ValueError("oxDNA export requires at least one base.")
    if len(molecules) != 1:
        raise ValueError("oxDNA export currently supports only all-DNA or all-RNA projects.")
    expected = "DNA2" if molecules == {"DNA"} else "RNA"
    model = model or expected
    if model != expected:
        raise ValueError(f"This project requires the {expected} export model, not {model}.")
    frames = _orientation_frames(project)
    offsets = BACKBONE_OFFSETS[model]
    configurations = {}
    for index, base in project.bases.items():
        a1, a3 = frames[index]
        a2 = cross(a3, a1)
        offset = add(add(scale(a1, offsets[0]), scale(a2, offsets[1])), scale(a3, offsets[2]))
        center = sub(scale(base.position, 1 / NM_PER_OXDNA_UNIT), offset)
        configurations[index] = (center, a1, a3)
    return configurations


def _export_sequence(project, default_base):
    assignments = {index: base.nucleotide for index, base in project.bases.items()}
    rng = random.Random(0)
    for index in sorted(project.bases):
        base = project.bases[index]
        alphabet = "AUCG" if base.molecule == "RNA" else "ATCG"
        if assignments[index] is None:
            if default_base is None:
                raise ValueError(f"Base {index} has no assigned nucleotide. Assign all bases before oxDNA export.")
            choice = default_base.upper()
            if choice not in (*alphabet, "R"):
                raise ValueError(f"Invalid default base {choice!r} for {base.molecule}.")
            if base.across is not None and assignments[base.across] is not None:
                assignments[index] = complement_for(assignments[base.across], base.molecule)
            else:
                assignments[index] = rng.choice(alphabet) if choice == "R" else choice
                if base.across is not None:
                    assignments[base.across] = complement_for(assignments[index], project.bases[base.across].molecule)
        if assignments[index] not in alphabet:
            raise ValueError(f"Base {index} has an invalid nucleotide for {base.molecule}.")
    return assignments


def export_oxdna_files(project: TiamatProject, output_path: str | Path, *, default_base: str | None = None,
                       tiamat_version: int = 2, model: str | None = None) -> tuple[Path, Path]:
    if tiamat_version != 2:
        raise ValueError("Coordinates must be in nm (Tiamat version 2). Legacy empirical scaling is no longer applied.")
    configurations = build_oxdna_configuration(project, model=model)
    assignments = _export_sequence(project, default_base)
    project.refresh_strands()
    # Classic topology rows are 3' -> 5'; up is 5', down is 3'.
    ordered = [(strand.strand_id, index) for strand in project.strands for index in reversed(strand.base_indices)]
    local_ids = {index: local for local, (_, index) in enumerate(ordered)}
    topology = [f"{len(project.bases)} {len(project.strands)}"]
    extent = max(abs(v) for center, _, _ in configurations.values() for v in center)
    box = max(250, math.ceil(2 * extent + 20))
    configuration = ["t = 0", f"b = {box} {box} {box}", "E = 0 0 0"]
    for strand_id, index in ordered:
        base = project.bases[index]
        topology.append(f"{strand_id} {assignments[index]} {local_ids.get(base.down, -1)} {local_ids.get(base.up, -1)}")
        center, a1, a3 = configurations[index]
        configuration.append(" ".join(f"{v:.12g}" for v in (*center, *a1, *a3)) + " 0 0 0 0 0 0")
    top, conf = resolve_oxdna_output_paths(output_path)
    atomic_write_files({top: "\n".join(topology) + "\n", conf: "\n".join(configuration) + "\n"})
    return top, conf
