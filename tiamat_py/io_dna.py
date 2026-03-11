from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import math
from pathlib import Path
import struct

from .io_json import load_json_project
from .model import Base, TiamatProject


NUCLEOBASE_CLASS_MARKER = b"\xff\xff\x05\x00\x0a\x00Nucleobase"
DNA_TYPE_TO_BASE = {
    0: "A",
    1: "T",
    2: "U",
    3: "C",
    4: "G",
}
GENERIC_BASE_TYPE = 5
SIDECAR_TOLERANCE = 1e-2
PARTIAL_SIDECAR_MATCH_MIN_FRACTION = 0.75
ACROSS_DISTANCE_MIN = 1.75
ACROSS_DISTANCE_MAX = 2.05
ORDER_ACROSS_DISTANCE_MAX = 1.95
BACKBONE_DISTANCE_MIN = 0.35
BACKBONE_DISTANCE_MAX = 0.95
ORDER_BACKBONE_DISTANCE_MAX = 6.0


@dataclass(frozen=True)
class DnaRecord:
    offset: int
    x: float
    y: float
    z: float
    nucleotide: str | None
    molecule: str | None
    flag: int


def load_dna_project(path: str | Path) -> TiamatProject:
    file_path = Path(path)
    records = parse_dna_records(file_path)
    sidecar_project, sidecar_path = _load_matching_sidecar(file_path, records)
    if sidecar_project is not None:
        sidecar_project.metadata["source_path"] = str(file_path)
        sidecar_project.metadata["source_format"] = "dna"
        sidecar_project.metadata["dna_sidecar_path"] = str(sidecar_path)
        sidecar_project.metadata["dna_record_count"] = len(records)
        return sidecar_project

    project = _build_fallback_project(records)
    project.metadata["source_path"] = str(file_path)
    project.metadata["source_format"] = "dna"
    project.metadata["dna_record_count"] = len(records)
    project.metadata["dna_topology"] = "record_order_fallback"
    return project


def parse_dna_records(path: str | Path) -> list[DnaRecord]:
    file_path = Path(path)
    data = file_path.read_bytes()
    candidate_records: list[list[DnaRecord]] = []
    for start in _scan_starts(data):
        records = _scan_records(data, start)
        if records:
            candidate_records.append(records)
    if not candidate_records:
        records = _scan_records(data, 0)
        if not records:
            raise ValueError("could not locate any Nucleobase records in the .dna archive")
        candidate_records.append(records)
    filtered = _filter_records_with_sidecar(file_path, candidate_records)
    if filtered is not None:
        return filtered
    return candidate_records[0]


def _scan_starts(data: bytes) -> list[int]:
    starts: list[int] = []
    position = data.find(NUCLEOBASE_CLASS_MARKER)
    if position != -1:
        starts.append(position + len(NUCLEOBASE_CLASS_MARKER))
    plain = data.find(b"Nucleobase")
    if plain != -1:
        starts.append(max(0, plain - 8))
        starts.append(plain)
    starts.append(0)
    unique: list[int] = []
    seen: set[int] = set()
    for start in starts:
        if start in seen:
            continue
        seen.add(start)
        unique.append(start)
    return unique


def _scan_records(data: bytes, start: int) -> list[DnaRecord]:
    records: list[DnaRecord] = []
    for offset in range(start, len(data) - 56):
        x, y, z = struct.unpack_from("<ddd", data, offset)
        if not _looks_like_coordinate_triplet(x, y, z):
            continue
        type_code = struct.unpack_from("<I", data, offset + 24)[0]
        flag = struct.unpack_from("<I", data, offset + 28)[0]
        parsed = _parse_record_fields(data, offset, type_code, flag)
        if parsed is None:
            continue
        nucleotide, molecule, parsed_flag = parsed
        records.append(
            DnaRecord(
                offset=offset,
                x=x,
                y=y,
                z=z,
                nucleotide=nucleotide,
                molecule=molecule,
                flag=parsed_flag,
            )
        )
    return records


def _parse_record_fields(
    data: bytes,
    offset: int,
    type_code: int,
    flag: int,
) -> tuple[str | None, str | None, int] | None:
    if type_code in DNA_TYPE_TO_BASE:
        if flag not in (0, 1):
            return None
        if not _matches_explicit_record_layout(data, offset):
            return None
        nucleotide = DNA_TYPE_TO_BASE[type_code]
        molecule = None
        if nucleotide == "U":
            molecule = "RNA"
        elif nucleotide == "T":
            molecule = "DNA"
        return nucleotide, molecule, flag

    if type_code == GENERIC_BASE_TYPE:
        if not _matches_generic_record_layout(data, offset):
            return None
        return None, None, 0

    return None


def _matches_explicit_record_layout(data: bytes, offset: int) -> bool:
    # Layout used by the first DNA sample:
    #   type @ +24, flag @ +28, 20 zero bytes, LE uint32(1) @ +52
    if data[offset + 32:offset + 52] == b"\x00" * 20 and struct.unpack_from("<I", data, offset + 52)[0] == 1:
        return True

    # Layout used by the RNA sample:
    #   type @ +24, flag @ +28, 16 zero bytes, BE uint32(1) @ +48, 4 zero bytes @ +52
    if (
        data[offset + 32:offset + 48] == b"\x00" * 16
        and struct.unpack_from(">I", data, offset + 48)[0] == 1
        and data[offset + 52:offset + 56] == b"\x00" * 4
    ):
        return True

    # Layout used by the standalone SSCCH test.dna sample:
    #   type @ +24, flag @ +28, 24 zero bytes through +55
    if data[offset + 32:offset + 56] == b"\x00" * 24:
        return True

    return False


def _matches_generic_record_layout(data: bytes, offset: int) -> bool:
    # Generic layout used by older DNA files:
    #   marker @ +24, 24 zero bytes, LE uint32(1) @ +52
    if data[offset + 28:offset + 52] == b"\x00" * 24 and struct.unpack_from("<I", data, offset + 52)[0] == 1:
        return True

    # Generic RNA layout used by current test.dna:
    #   marker @ +24, 20 zero bytes, BE uint32(1) @ +48, 4 zero bytes @ +52
    if (
        data[offset + 28:offset + 48] == b"\x00" * 20
        and struct.unpack_from(">I", data, offset + 48)[0] == 1
        and data[offset + 52:offset + 56] == b"\x00" * 4
    ):
        return True

    return False


def _looks_like_coordinate_triplet(x: float, y: float, z: float) -> bool:
    if not all(math.isfinite(value) and abs(value) < 1e6 for value in (x, y, z)):
        return False
    # Filter out denormal garbage that can accidentally satisfy the record-layout bytes.
    return max(abs(x), abs(y), abs(z)) >= 1e-4


def _load_matching_sidecar(
    dna_path: Path,
    records: list[DnaRecord],
) -> tuple[TiamatProject | None, Path | None]:
    best_exact: tuple[float, TiamatProject, Path] | None = None
    best_partial: tuple[float, float, TiamatProject, Path] | None = None
    for candidate in _sidecar_candidates(dna_path):
        try:
            project = load_json_project(candidate)
        except Exception:
            continue
        score = _match_sidecar_score(project, records)
        if score is not None:
            if best_exact is None or score < best_exact[0]:
                best_exact = (score, project, candidate)
            continue
        partial = _partial_sidecar_score(project, records)
        if partial is None:
            continue
        coverage, distance = partial
        if best_partial is None or coverage > best_partial[0] or (coverage == best_partial[0] and distance < best_partial[1]):
            best_partial = (coverage, distance, project, candidate)
    if best_exact is not None:
        return best_exact[1], best_exact[2]
    if best_partial is not None:
        return best_partial[2], best_partial[3]
    return None, None


def _filter_records_with_sidecar(
    dna_path: Path,
    candidate_records: list[list[DnaRecord]],
) -> list[DnaRecord] | None:
    for records in candidate_records:
        for candidate in _sidecar_candidates(dna_path):
            try:
                project = load_json_project(candidate)
            except Exception:
                continue
            if len(project.bases) >= len(records):
                continue
            filtered = _records_matching_project(project, records)
            if filtered is not None:
                return filtered
    return None


def _sidecar_candidates(dna_path: Path) -> list[Path]:
    parent = dna_path.parent
    candidates = [
        dna_path.with_suffix(".dnajson"),
        dna_path.with_suffix(".json"),
        parent / "output.dnajson",
        parent / "output.json",
    ]
    for extra in sorted(parent.glob("*.dnajson")):
        candidates.append(extra)
    for extra in sorted(parent.glob("*.json")):
        candidates.append(extra)

    seen: set[Path] = set()
    unique_candidates: list[Path] = []
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved in seen or not candidate.exists() or candidate == dna_path:
            continue
        seen.add(resolved)
        unique_candidates.append(candidate)
    return unique_candidates


def _match_sidecar_score(project: TiamatProject, records: list[DnaRecord]) -> float | None:
    if len(project.bases) != len(records):
        return None

    unmatched = {
        index: base
        for index, base in project.bases.items()
    }
    total_distance = 0.0
    sorted_records = sorted(records, key=lambda item: (item.x, item.y, item.z, item.offset))
    for record in sorted_records:
        best_index = None
        best_distance = None
        for index, base in unmatched.items():
            if (
                base.nucleotide is not None
                and record.nucleotide is not None
                and base.nucleotide != record.nucleotide
            ):
                continue
            if (
                base.molecule is not None
                and record.molecule is not None
                and base.molecule != record.molecule
            ):
                continue
            distance = math.dist((record.x, record.y, record.z), base.position)
            if distance > SIDECAR_TOLERANCE:
                continue
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_index = index
        if best_index is None or best_distance is None:
            return None
        total_distance += best_distance
        del unmatched[best_index]
    return total_distance


def _partial_sidecar_score(project: TiamatProject, records: list[DnaRecord]) -> tuple[float, float] | None:
    if not project.bases or not records:
        return None
    if len(records) >= len(project.bases):
        return None
    coverage = len(records) / len(project.bases)
    if coverage < PARTIAL_SIDECAR_MATCH_MIN_FRACTION:
        return None

    unmatched = {
        index: base
        for index, base in project.bases.items()
    }
    total_distance = 0.0
    sorted_records = sorted(records, key=lambda item: (item.x, item.y, item.z, item.offset))
    for record in sorted_records:
        best_index = None
        best_distance = None
        for index, base in unmatched.items():
            if (
                base.nucleotide is not None
                and record.nucleotide is not None
                and base.nucleotide != record.nucleotide
            ):
                continue
            if (
                base.molecule is not None
                and record.molecule is not None
                and base.molecule != record.molecule
            ):
                continue
            distance = math.dist((record.x, record.y, record.z), base.position)
            if distance > SIDECAR_TOLERANCE:
                continue
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_index = index
        if best_index is None or best_distance is None:
            return None
        total_distance += best_distance
        del unmatched[best_index]
    return coverage, total_distance


def _records_matching_project(project: TiamatProject, records: list[DnaRecord]) -> list[DnaRecord] | None:
    unmatched = list(records)
    matched_records: list[DnaRecord] = []
    for index in sorted(project.bases):
        base = project.bases[index]
        best_position = None
        best_distance = None
        for position, record in enumerate(unmatched):
            if (
                base.nucleotide is not None
                and record.nucleotide is not None
                and base.nucleotide != record.nucleotide
            ):
                continue
            if (
                base.molecule is not None
                and record.molecule is not None
                and base.molecule != record.molecule
            ):
                continue
            distance = math.dist((record.x, record.y, record.z), base.position)
            if distance > SIDECAR_TOLERANCE:
                continue
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_position = position
        if best_position is None:
            return None
        matched_records.append(unmatched.pop(best_position))
    return matched_records


def _build_fallback_project(records: list[DnaRecord]) -> TiamatProject:
    bases = {
        index: Base(
            index=index,
            object_id=f"dna:{record.offset}",
            x=record.x,
            y=record.y,
            z=record.z,
            molecule=record.molecule or "",
            nucleotide=record.nucleotide,
        )
        for index, record in enumerate(records)
    }
    _assign_record_order_topology(records, bases)
    project = TiamatProject(
        bases=bases,
        metadata={"format": "dna_fixed_records"},
    )
    _assign_strand_molecules(project)
    return project


def _assign_record_order_topology(records: list[DnaRecord], bases: dict[int, Base]) -> None:
    for chunk in _record_order_chunks(records):
        _assign_record_order_chunk(records, bases, chunk)


def _record_order_chunks(records: list[DnaRecord]) -> list[tuple[int, int]]:
    if not records:
        return []

    chunks: list[tuple[int, int]] = []
    start = 0
    for index in range(len(records) - 1):
        if _record_distance(records, index, index + 1) > ORDER_BACKBONE_DISTANCE_MAX:
            chunks.append((start, index + 1))
            start = index + 1
    chunks.append((start, len(records)))
    return chunks


def _assign_record_order_chunk(
    records: list[DnaRecord],
    bases: dict[int, Base],
    chunk: tuple[int, int],
) -> None:
    start, end = chunk
    positions = list(range(start, end))
    if not positions:
        return

    colors: dict[int, int] = {positions[0]: 0}
    for left, right in zip(positions, positions[1:]):
        relation = _classify_record_order_relation(records, left, right)
        if relation == "across" and bases[left].across is None and bases[right].across is None:
            bases[left].across = right
            bases[right].across = left
        if relation == "across":
            colors[right] = 1 - colors[left]
        else:
            colors[right] = colors[left]

    wrap_relation = None
    if len(positions) > 1:
        wrap_relation = _classify_record_order_relation(records, positions[-1], positions[0])
        if wrap_relation == "across" and bases[positions[-1]].across is None and bases[positions[0]].across is None:
            bases[positions[-1]].across = positions[0]
            bases[positions[0]].across = positions[-1]
        expected_color = colors[positions[-1]] if wrap_relation != "across" else 1 - colors[positions[-1]]
        if wrap_relation == "break" or expected_color != colors[positions[0]]:
            wrap_relation = None

    strands_by_color: dict[int, list[int]] = defaultdict(list)
    for position in positions:
        strands_by_color[colors[position]].append(position)

    directions = _infer_track_directions(records, strands_by_color)
    for color, strand_positions in strands_by_color.items():
        direction = directions.get(color, "forward")
        for left, right in zip(strand_positions, strand_positions[1:]):
            _attach_oriented_backbone_pair(bases, left, right, direction)
        if (
            wrap_relation == "backbone"
            and len(strand_positions) > 1
            and strand_positions[0] == positions[0]
            and strand_positions[-1] == positions[-1]
        ):
            _attach_oriented_backbone_pair(bases, strand_positions[-1], strand_positions[0], direction)


def _classify_record_order_relation(records: list[DnaRecord], left: int, right: int) -> str:
    distance = _record_distance(records, left, right)
    if ACROSS_DISTANCE_MIN <= distance <= ORDER_ACROSS_DISTANCE_MAX:
        return "across"
    if distance <= ORDER_BACKBONE_DISTANCE_MAX:
        return "backbone"
    return "break"


def _record_distance(records: list[DnaRecord], left: int, right: int) -> float:
    first = records[left]
    second = records[right]
    return math.dist((first.x, first.y, first.z), (second.x, second.y, second.z))


def _infer_track_directions(
    records: list[DnaRecord],
    strands_by_color: dict[int, list[int]],
) -> dict[int, str]:
    directions = {
        color: _infer_track_direction(records, strand_positions)
        for color, strand_positions in strands_by_color.items()
    }
    if len(directions) == 2:
        colors = sorted(directions)
        left_direction = directions[colors[0]]
        right_direction = directions[colors[1]]
        if left_direction is None and right_direction is not None:
            directions[colors[0]] = _opposite_direction(right_direction)
        elif right_direction is None and left_direction is not None:
            directions[colors[1]] = _opposite_direction(left_direction)
    return {
        color: direction or "forward"
        for color, direction in directions.items()
    }


def _infer_track_direction(
    records: list[DnaRecord],
    strand_positions: list[int],
) -> str | None:
    if len(strand_positions) < 2:
        return None
    same_track_gaps = [
        records[right].offset - (records[left].offset + 56)
        for left, right in zip(strand_positions, strand_positions[1:])
    ]
    for gap in same_track_gaps:
        if gap == 0:
            return "forward"
        if gap == 3:
            return "reverse"
    return None


def _opposite_direction(direction: str) -> str:
    return "reverse" if direction == "forward" else "forward"


def _assign_geometric_across_links(bases: dict[int, Base]) -> None:
    nearest_candidates: dict[int, tuple[float, int]] = {}
    indices = sorted(bases)
    for index in indices:
        base = bases[index]
        best: tuple[float, int] | None = None
        for other_index in indices:
            if other_index == index:
                continue
            other = bases[other_index]
            distance = math.dist(base.position, other.position)
            if not (ACROSS_DISTANCE_MIN <= distance <= ACROSS_DISTANCE_MAX):
                continue
            if best is None or distance < best[0]:
                best = (distance, other_index)
        if best is not None:
            nearest_candidates[index] = best

    for index, (_distance, other_index) in nearest_candidates.items():
        candidate = nearest_candidates.get(other_index)
        if candidate is None or candidate[1] != index:
            continue
        if bases[index].across is None and bases[other_index].across is None:
            bases[index].across = other_index
            bases[other_index].across = index


def _assign_local_backbone_links(bases: dict[int, Base]) -> None:
    candidates: list[tuple[float, int, int]] = []
    indices = sorted(bases)
    for position, index in enumerate(indices):
        base = bases[index]
        for other_index in indices[position + 1:]:
            other = bases[other_index]
            if base.across == other_index:
                continue
            distance = math.dist(base.position, other.position)
            if BACKBONE_DISTANCE_MIN <= distance <= BACKBONE_DISTANCE_MAX:
                candidates.append((distance, index, other_index))

    degree = {index: 0 for index in indices}
    for _distance, left, right in sorted(candidates):
        if degree[left] >= 2 or degree[right] >= 2:
            continue
        if _already_connected(bases[left], right) or _already_connected(bases[right], left):
            continue
        _attach_backbone_neighbor(bases[left], right)
        _attach_backbone_neighbor(bases[right], left)
        degree[left] += 1
        degree[right] += 1


def _assign_order_backbone_links(records: list[DnaRecord], bases: dict[int, Base]) -> None:
    for left, right in zip(range(len(records) - 1), range(1, len(records))):
        first = bases[left]
        second = bases[right]
        if first.across == right or second.across == left:
            continue
        if _already_connected(first, right) or _already_connected(second, left):
            continue
        if first.up is not None and first.down is not None:
            continue
        if second.up is not None and second.down is not None:
            continue
        distance = math.dist(first.position, second.position)
        if not (BACKBONE_DISTANCE_MAX < distance <= ORDER_BACKBONE_DISTANCE_MAX):
            continue
        _attach_backbone_neighbor(first, right)
        _attach_backbone_neighbor(second, left)


def _attach_backbone_neighbor(base: Base, neighbor: int) -> None:
    if base.up is None:
        base.up = neighbor
        return
    if base.down is None:
        base.down = neighbor


def _attach_backbone_pair(bases: dict[int, Base], left: int, right: int) -> None:
    if _already_connected(bases[left], right) or _already_connected(bases[right], left):
        return
    if bases[left].up is not None and bases[left].down is not None:
        return
    if bases[right].up is not None and bases[right].down is not None:
        return
    _attach_backbone_neighbor(bases[left], right)
    _attach_backbone_neighbor(bases[right], left)


def _attach_oriented_backbone_pair(
    bases: dict[int, Base],
    left: int,
    right: int,
    direction: str,
) -> None:
    if _already_connected(bases[left], right) or _already_connected(bases[right], left):
        return
    if direction == "reverse":
        _set_oriented_backbone_edge(bases, right, left)
        return
    _set_oriented_backbone_edge(bases, left, right)


def _set_oriented_backbone_edge(
    bases: dict[int, Base],
    up_index: int,
    down_index: int,
) -> None:
    up_base = bases[up_index]
    down_base = bases[down_index]
    if up_base.down is not None or down_base.up is not None:
        return
    up_base.down = down_index
    down_base.up = up_index


def _assign_strand_molecules(project: TiamatProject) -> None:
    for strand in project.strands:
        molecule = _strand_molecule(project, strand.base_indices)
        for index in strand.base_indices:
            project.bases[index].molecule = molecule


def _already_connected(base: Base, neighbor: int) -> bool:
    return neighbor in {base.up, base.down}


def _infer_default_molecule(records: list[DnaRecord]) -> str:
    nucleotides = {record.nucleotide for record in records}
    if "U" in nucleotides:
        return "RNA"
    return "DNA"


def _strand_molecule(project: TiamatProject, base_indices: list[int]) -> str:
    explicit = {project.bases[index].molecule for index in base_indices if project.bases[index].molecule}
    if "RNA" in explicit:
        return "RNA"
    if "DNA" in explicit:
        return "DNA"

    nucleotides = {project.bases[index].nucleotide for index in base_indices}
    if "U" in nucleotides:
        return "RNA"
    if "T" in nucleotides:
        return "DNA"
    return "DNA"
