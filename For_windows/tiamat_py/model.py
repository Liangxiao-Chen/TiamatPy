from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
import colorsys
import json
import math
import random
import re
import sys
import time
from typing import Any


NUCLEOTIDES = ("A", "C", "G", "T")
ALLOWED_BASES = ("A", "C", "G", "T", "U")
FULL_NAME_TO_BASE = {
    "Adenine": "A",
    "Cytosine": "C",
    "Guanine": "G",
    "Thymine": "T",
    "Uracil": "U",
}
BASE_COLORS = {
    "A": "#22c55e",
    "C": "#3b82f6",
    "G": "#facc15",
    "T": "#ef4444",
    "U": "#a855f7",
    "N": "#9ca3af",
    None: "#9ca3af",
}


@dataclass
class Base:
    index: int
    object_id: str
    x: float
    y: float
    z: float
    molecule: str = "DNA"
    up: int | None = None
    down: int | None = None
    across: int | None = None
    nucleotide: str | None = None
    strand_id: int | None = None
    backbone_color: str | None = None
    plane_u: tuple[float, float, float] | None = None
    plane_v: tuple[float, float, float] | None = None

    @property
    def position(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "object_id": self.object_id,
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "molecule": self.molecule,
            "up": self.up,
            "down": self.down,
            "across": self.across,
            "nucleotide": self.nucleotide,
            "strand_id": self.strand_id,
            "backbone_color": self.backbone_color,
            "plane_u": list(self.plane_u) if self.plane_u is not None else None,
            "plane_v": list(self.plane_v) if self.plane_v is not None else None,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Base":
        return cls(
            index=int(data["index"]),
            object_id=str(data["object_id"]),
            x=float(data["x"]),
            y=float(data["y"]),
            z=float(data["z"]),
            molecule=str(data.get("molecule", "DNA")).upper(),
            up=_as_int_or_none(data.get("up")),
            down=_as_int_or_none(data.get("down")),
            across=_as_int_or_none(data.get("across")),
            nucleotide=_normalize_base(data.get("nucleotide")),
            strand_id=_as_int_or_none(data.get("strand_id")),
            backbone_color=_normalize_color(data.get("backbone_color")),
            plane_u=_as_vector_or_none(data.get("plane_u")),
            plane_v=_as_vector_or_none(data.get("plane_v")),
        )


@dataclass
class Strand:
    strand_id: int
    base_indices: list[int]
    color: str

    def sequence(self, project: "TiamatProject") -> str:
        return "".join(project.bases[index].nucleotide or "N" for index in self.base_indices)


@dataclass
class TiamatProject:
    bases: dict[int, Base]
    metadata: dict[str, Any] = field(default_factory=dict)
    strands: list[Strand] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.strands:
            self.refresh_strands()

    def refresh_strands(self) -> None:
        self.strands = []
        components = self._strand_components()
        for strand_id, component in enumerate(components, start=1):
            ordered = self._order_component(component)
            color = _strand_color(strand_id)
            self.strands.append(Strand(strand_id=strand_id, base_indices=ordered, color=color))
            for base_index in ordered:
                self.bases[base_index].strand_id = strand_id

    def summary(self) -> dict[str, int]:
        paired = sum(1 for base in self.bases.values() if base.across is not None)
        assigned = sum(1 for base in self.bases.values() if base.nucleotide)
        return {
            "bases": len(self.bases),
            "strands": len(self.strands),
            "paired_bases": paired,
            "assigned_bases": assigned,
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "format": "tiamat_py_v1",
            "metadata": self.metadata,
            "bases": [self.bases[index].to_dict() for index in sorted(self.bases)],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "TiamatProject":
        bases = {item["index"]: Base.from_dict(item) for item in data["bases"]}
        return cls(bases=bases, metadata=dict(data.get("metadata", {})))

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    def along_neighbors(self, index: int) -> list[int]:
        base = self.bases[index]
        return sorted({neighbor for neighbor in (base.up, base.down) if neighbor is not None})

    def all_edges(self) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
        along: set[tuple[int, int]] = set()
        across: set[tuple[int, int]] = set()
        for index, base in self.bases.items():
            for neighbor in (base.up, base.down):
                if neighbor is not None:
                    along.add(tuple(sorted((index, neighbor))))
            if base.across is not None:
                across.add(tuple(sorted((index, base.across))))
        return sorted(along), sorted(across)

    def backbone_edges(self) -> list[tuple[int, int]]:
        along: set[tuple[int, int]] = set()
        for index, base in self.bases.items():
            if base.down is not None and base.down in self.bases:
                along.add((index, base.down))
            elif base.up is not None and base.up in self.bases and self.bases[base.up].down != index:
                along.add((base.up, index))
        return sorted(along)

    def delete_bases(self, indices: list[int] | set[int]) -> None:
        targets = {index for index in indices if index in self.bases}
        if not targets:
            return

        remaining = set(self.bases) - targets
        for index in sorted(remaining):
            base = self.bases[index]
            if base.up is not None and base.up not in remaining:
                base.up = None
            if base.down is not None and base.down not in remaining:
                base.down = None
                base.backbone_color = None
            if base.across is not None and base.across not in remaining:
                base.across = None

        for index in targets:
            self.bases.pop(index, None)

        self._set_sticky_end_groups(self.sticky_end_groups())

        if self.bases:
            self.refresh_strands()
        else:
            self.strands = []

    def reset_bases(self, indices: list[int] | None = None) -> None:
        targets = indices or list(self.bases)
        for index in targets:
            self.bases[index].nucleotide = None

    def set_backbone_color(self, indices: list[int] | set[int], color: str) -> int:
        normalized = _normalize_color(color)
        if normalized is None:
            return 0
        changed = 0
        for index in sorted(indices):
            base = self.bases.get(index)
            if base is None or base.down is None or base.down not in self.bases:
                continue
            if base.backbone_color == normalized:
                continue
            base.backbone_color = normalized
            changed += 1
        return changed

    def reset_backbone_color(self, indices: list[int] | set[int]) -> int:
        changed = 0
        for index in sorted(indices):
            base = self.bases.get(index)
            if base is None:
                continue
            if base.backbone_color is None:
                continue
            base.backbone_color = None
            changed += 1
        return changed

    def fill_unpaired_generic_bases(self, nucleotide: str) -> int:
        normalized = _normalize_fill_base(nucleotide)
        filled = 0
        for base in self.bases.values():
            if base.nucleotide is not None:
                continue
            if base.across is not None and base.across in self.bases:
                continue
            base.nucleotide = _fill_base_for_molecule(normalized, base.molecule)
            filled += 1
        return filled

    def can_create_complementary(self, indices: list[int] | set[int]) -> bool:
        selected = self._selected_pair(indices)
        if selected is None:
            return False
        left, right = (self.bases[index] for index in selected)
        if left.across is not None or right.across is not None:
            return False
        return bases_can_pair(left, right)

    def create_complementary(self, indices: list[int] | set[int]) -> bool:
        selected = self._selected_pair(indices)
        if selected is None or not self.can_create_complementary(selected):
            return False
        left_index, right_index = selected
        self.bases[left_index].across = right_index
        self.bases[right_index].across = left_index
        return True

    def can_delete_complementary(self, indices: list[int] | set[int]) -> bool:
        selected = self._selected_pair(indices)
        if selected is None:
            return False
        left_index, right_index = selected
        return self.bases[left_index].across == right_index and self.bases[right_index].across == left_index

    def delete_complementary(self, indices: list[int] | set[int]) -> bool:
        selected = self._selected_pair(indices)
        if selected is None or not self.can_delete_complementary(selected):
            return False
        left_index, right_index = selected
        self.bases[left_index].across = None
        self.bases[right_index].across = None
        return True

    def can_create_down(self, indices: list[int] | set[int]) -> bool:
        selected = self._selected_pair(indices)
        if selected is None:
            return False
        left_index, right_index = selected
        left = self.bases[left_index]
        right = self.bases[right_index]
        return (
            (left.down is None and right.up is None)
            or (right.down is None and left.up is None)
        )

    def create_down(self, indices: list[int] | set[int]) -> bool:
        selected = self._selected_pair(indices)
        if selected is None or not self.can_create_down(selected):
            return False
        left_index, right_index = selected
        left = self.bases[left_index]
        right = self.bases[right_index]
        if left.down is None and right.up is None:
            left.down = right_index
            right.up = left_index
        else:
            right.down = left_index
            left.up = right_index
        self.refresh_strands()
        return True

    def can_delete_down(self, indices: list[int] | set[int]) -> bool:
        selected = sorted(index for index in indices if index in self.bases)
        return len(selected) == 1 and self.bases[selected[0]].down is not None and self.bases[selected[0]].down in self.bases

    def delete_down(self, indices: list[int] | set[int]) -> bool:
        selected = sorted(index for index in indices if index in self.bases)
        if not self.can_delete_down(selected):
            return False
        index = selected[0]
        down_index = self.bases[index].down
        self.bases[index].down = None
        self.bases[index].backbone_color = None
        if down_index is not None and down_index in self.bases and self.bases[down_index].up == index:
            self.bases[down_index].up = None
        self.refresh_strands()
        return True

    def sticky_end_groups(self) -> list[dict[str, Any]]:
        groups = self.metadata.get("sticky_ends")
        if not isinstance(groups, list):
            return []
        normalized: list[dict[str, Any]] = []
        seen_ids: set[int] = set()
        next_id = 1
        for item in groups:
            if not isinstance(item, dict):
                continue
            left_indices = _normalize_index_list(item.get("left_indices"), self.bases)
            right_indices = _normalize_index_list(item.get("right_indices"), self.bases)
            if not left_indices or len(left_indices) != len(right_indices):
                continue
            sticky_id = _as_int_or_none(item.get("id"))
            while sticky_id is None or sticky_id in seen_ids:
                sticky_id = next_id
                next_id += 1
            seen_ids.add(sticky_id)
            color = _normalize_color(item.get("color")) or _sticky_end_color(sticky_id)
            normalized.append(
                {
                    "id": sticky_id,
                    "left_indices": left_indices,
                    "right_indices": right_indices,
                    "color": color,
                }
            )
        return normalized

    def can_create_sticky_end(self, indices: list[int] | set[int]) -> bool:
        segments = _sticky_end_segments_for_selection(self, indices)
        if segments is None:
            return False
        used_indices = self._sticky_end_index_map()
        return not any(index in used_indices for segment in segments for index in segment)

    def create_sticky_end(self, indices: list[int] | set[int]) -> tuple[bool, str | None]:
        segments = _sticky_end_segments_for_selection(self, indices)
        if segments is None:
            return False, "Select two continuous nucleotide segments of the same length."
        used_indices = self._sticky_end_index_map()
        if any(index in used_indices for segment in segments for index in segment):
            return False, "Selected bases are already part of a sticky end."
        left_indices, right_indices = segments
        if not _sticky_end_sequences_match(self, left_indices, right_indices):
            return False, "Selected sequences are not reverse complementary."
        groups = self.sticky_end_groups()
        sticky_id = max((group["id"] for group in groups), default=0) + 1
        groups.append(
            {
                "id": sticky_id,
                "left_indices": list(left_indices),
                "right_indices": list(right_indices),
                "color": _sticky_end_color(sticky_id),
            }
        )
        self._set_sticky_end_groups(groups)
        return True, None

    def can_delete_sticky_end(self, indices: list[int] | set[int]) -> bool:
        return self._selected_sticky_end_id(indices) is not None

    def delete_sticky_end(self, indices: list[int] | set[int]) -> bool:
        sticky_id = self._selected_sticky_end_id(indices)
        if sticky_id is None:
            return False
        groups = [group for group in self.sticky_end_groups() if group["id"] != sticky_id]
        self._set_sticky_end_groups(groups)
        return True

    def _selected_sticky_end_id(self, indices: list[int] | set[int]) -> int | None:
        sticky_index_map = self._sticky_end_index_map()
        selected_group_ids = {
            sticky_index_map[index]
            for index in indices
            if index in self.bases and index in sticky_index_map
        }
        selected_indices = {index for index in indices if index in self.bases}
        if not selected_indices or len(selected_group_ids) != 1:
            return None
        if len(selected_group_ids) == 1 and len(selected_indices) == len(
            {index for index in selected_indices if index in sticky_index_map}
        ):
            return next(iter(selected_group_ids))
        return None

    def _sticky_end_index_map(self) -> dict[int, int]:
        mapping: dict[int, int] = {}
        for group in self.sticky_end_groups():
            for index in group["left_indices"] + group["right_indices"]:
                mapping[index] = group["id"]
        return mapping

    def _set_sticky_end_groups(self, groups: list[dict[str, Any]]) -> None:
        if groups:
            self.metadata["sticky_ends"] = groups
        else:
            self.metadata.pop("sticky_ends", None)

    def set_nucleotide(self, index: int, nucleotide: str | None) -> None:
        base = self.bases[index]
        base.nucleotide = _normalize_base(nucleotide)
        if base.across is None or base.across not in self.bases:
            return
        across = self.bases[base.across]
        if base.nucleotide is None:
            across.nucleotide = None
            return
        if basepair_state(base, across) == "wobble":
            return
        across.nucleotide = complement_for(base.nucleotide, across.molecule)

    def generate_sequences(
        self,
        seed: int | None = None,
        unique_sequence_limit: int = 10,
        repetition_limit: int = 6,
        gc_percentage: float = 0.5,
        gc_tolerance: float = 0.05,
        timeout: float = 16.0,
    ) -> None:
        if not self.bases:
            return
        if unique_sequence_limit < 4:
            raise ValueError("Unique Sequence Limit must be at least 4.")
        if repetition_limit < 2:
            raise ValueError("Repetition Limit must be at least 2.")
        if not 0.0 <= gc_percentage <= 1.0:
            raise ValueError("G/C Percentage must be between 0 and 1.")
        if not 0.0 <= gc_tolerance <= 0.5:
            raise ValueError("G/C tolerance must be between 0 and 0.5.")
        if timeout <= 0.0:
            raise ValueError("Operation Timeout must be greater than 0.")

        previous = {index: base.nucleotide for index, base in self.bases.items()}
        rng = random.Random(seed)
        deadline = time.monotonic() + timeout
        min_gc = max(0, math.ceil((gc_percentage - gc_tolerance) * len(self.bases)))
        max_gc = min(len(self.bases), math.floor((gc_percentage + gc_tolerance) * len(self.bases)))
        if min_gc > max_gc:
            raise ValueError("G/C Percentage and tolerance leave no valid GC range.")

        try:
            assignments = _try_periodic_sequence_assignments(
                project=self,
                rng=rng,
                unique_limit=unique_sequence_limit,
                repetition_limit=repetition_limit,
                min_gc=min_gc,
                max_gc=max_gc,
            )
            if assignments is None:
                assignments = _generate_sequence_assignments(
                    project=self,
                    rng=rng,
                    deadline=deadline,
                    unique_limit=unique_sequence_limit,
                    repetition_limit=repetition_limit,
                    min_gc=min_gc,
                    max_gc=max_gc,
                )
        except TimeoutError:
            for index, nucleotide in previous.items():
                self.bases[index].nucleotide = nucleotide
            raise

        for index, base in self.bases.items():
            base.nucleotide = assignments[index]

    def export_sequences_text(self) -> str:
        lines = []
        for strand in self.strands:
            lines.append(f">strand_{strand.strand_id} length={len(strand.base_indices)}")
            lines.append(strand.sequence(self))
        return "\n".join(lines) + "\n"

    def describe_base(self, index: int) -> str:
        base = self.bases[index]
        return (
            f"Base {base.index}\n"
            f"  strand: {base.strand_id}\n"
            f"  nucleotide: {base.nucleotide or 'N'}\n"
            f"  up/down/across: {base.up or 0} / {base.down or 0} / {base.across or 0}\n"
            f"  position: ({base.x:.3f}, {base.y:.3f}, {base.z:.3f})"
        )

    def _strand_components(self) -> list[set[int]]:
        remaining = set(self.bases)
        components: list[set[int]] = []
        while remaining:
            start = min(remaining)
            stack = [start]
            component: set[int] = set()
            while stack:
                current = stack.pop()
                if current in component:
                    continue
                component.add(current)
                remaining.discard(current)
                stack.extend(self.along_neighbors(current))
            components.append(component)
        components.sort(key=lambda item: (min(item), len(item)))
        return components

    def _order_component(self, component: set[int]) -> list[int]:
        if len(component) == 1:
            return [next(iter(component))]
        ordered: list[int] = []
        remaining = set(component)
        while remaining:
            start_candidates = self._component_start_candidates(remaining)
            start = start_candidates[0] if start_candidates else min(remaining)
            segment = self._walk_down_chain(start, remaining)
            if not segment:
                ordered.extend(sorted(remaining))
                break
            ordered.extend(segment)
            remaining.difference_update(segment)
        return ordered

    def _component_start_candidates(self, component: set[int]) -> list[int]:
        candidates = []
        for index in component:
            up = self.bases[index].up
            if up is None or up not in component:
                candidates.append(index)
        return sorted(candidates)

    def _walk_down_chain(self, start: int, component: set[int]) -> list[int]:
        ordered: list[int] = []
        current = start
        visited: set[int] = set()
        while current in component and current not in visited:
            ordered.append(current)
            visited.add(current)
            down = self.bases[current].down
            if down is None or down not in component or down in visited:
                break
            current = down
        return ordered

    def _selected_pair(self, indices: list[int] | set[int]) -> tuple[int, int] | None:
        selected = sorted(index for index in indices if index in self.bases)
        if len(selected) != 2:
            return None
        return (selected[0], selected[1])


def _as_int_or_none(value: Any) -> int | None:
    if value in (None, "", 0, "0"):
        return None
    return int(value)


def _normalize_base(value: Any) -> str | None:
    if value in (None, "", "N"):
        return None
    base = str(value).upper()
    base = FULL_NAME_TO_BASE.get(base.capitalize(), base)
    if base not in ALLOWED_BASES:
        raise ValueError(f"unsupported nucleotide: {value}")
    return base


def _normalize_fill_base(value: Any) -> str:
    base = _normalize_base(value)
    if base is None:
        raise ValueError("fill base must be A, T/U, C, or G")
    if base not in {"A", "C", "G", "T", "U"}:
        raise ValueError(f"unsupported fill base: {value}")
    return base


def _fill_base_for_molecule(base: str, molecule: str) -> str:
    normalized = _normalize_fill_base(base)
    if normalized in {"T", "U"}:
        return "U" if str(molecule).upper() == "RNA" else "T"
    return normalized


def _normalize_color(value: Any) -> str | None:
    if value in (None, ""):
        return None
    text = str(value).strip()
    if not text:
        return None
    if re.fullmatch(r"#[0-9a-fA-F]{6}", text):
        return text.lower()
    return text


def _as_vector_or_none(value: Any) -> tuple[float, float, float] | None:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        return None
    try:
        return (float(value[0]), float(value[1]), float(value[2]))
    except (TypeError, ValueError):
        return None


def complement_for(base: str, target_molecule: str = "DNA") -> str:
    normalized = _normalize_base(base)
    if normalized is None:
        raise ValueError("cannot compute complement for undefined base")
    if normalized == "A":
        return "U" if str(target_molecule).upper() == "RNA" else "T"
    if normalized in {"T", "U"}:
        return "A"
    if normalized == "C":
        return "G"
    if normalized == "G":
        return "C"
    raise ValueError(f"unsupported nucleotide: {base}")


def basepair_state(left: Base, right: Base) -> str:
    if left.nucleotide is None or right.nucleotide is None:
        return "undefined"
    if (
        str(left.molecule).upper() == "RNA"
        and str(right.molecule).upper() == "RNA"
        and {left.nucleotide, right.nucleotide} == {"G", "U"}
    ):
        return "wobble"
    try:
        if (
            right.nucleotide == complement_for(left.nucleotide, right.molecule)
            and left.nucleotide == complement_for(right.nucleotide, left.molecule)
        ):
            return "canonical"
    except ValueError:
        return "invalid"
    return "invalid"


def bases_can_pair(left: Base, right: Base) -> bool:
    return basepair_state(left, right) in {"canonical", "wobble", "undefined"}


def _try_periodic_sequence_assignments(
    project: TiamatProject,
    rng: random.Random,
    unique_limit: int,
    repetition_limit: int,
    min_gc: int,
    max_gc: int,
) -> dict[int, str] | None:
    if 2 * len(project.bases) < 4 * min_gc or 2 * len(project.bases) > 4 * max_gc:
        return None

    base_to_location = {
        base_index: (strand.strand_id, position)
        for strand in project.strands
        for position, base_index in enumerate(strand.base_indices)
    }
    pair_constraints: dict[int, dict[int, list[tuple[int, int]]]] = {strand.strand_id: {} for strand in project.strands}
    for index, base in project.bases.items():
        across = base.across
        if across is None or across not in project.bases or index > across:
            continue
        left_strand, left_pos = base_to_location[index]
        right_strand, right_pos = base_to_location[across]
        pair_constraints[left_strand].setdefault(right_strand, []).append((left_pos, right_pos))
        pair_constraints[right_strand].setdefault(left_strand, []).append((right_pos, left_pos))

    cycle_families = [
        {
            "DNA": ("A", "C", "T", "G"),
            "RNA": ("A", "C", "U", "G"),
        },
        {
            "DNA": ("A", "G", "T", "C"),
            "RNA": ("A", "G", "U", "C"),
        },
    ]
    order = list(cycle_families)
    rng.shuffle(order)
    strand_ids = [strand.strand_id for strand in project.strands]
    components = _connected_strand_components(strand_ids, pair_constraints)

    for family in order:
        states: dict[int, tuple[int, int]] = {}
        valid = True
        randomized_components = list(components)
        rng.shuffle(randomized_components)
        for component in randomized_components:
            if not _solve_periodic_component(component, pair_constraints, states, rng):
                valid = False
                break
        if not valid:
            continue

        assignments: dict[int, str] = {}
        for strand in project.strands:
            sign, offset = states[strand.strand_id]
            for position, base_index in enumerate(strand.base_indices):
                molecule = str(project.bases[base_index].molecule).upper()
                cycle = family["RNA"] if molecule == "RNA" else family["DNA"]
                assignments[base_index] = cycle[(offset + sign * position) % 4]

        gc_count = sum(1 for nucleotide in assignments.values() if nucleotide in {"G", "C"})
        if not (min_gc <= gc_count <= max_gc):
            continue
        if not _assignments_respect_pairs(project, assignments):
            continue
        if not _assignments_respect_sequence_constraints(project, assignments, unique_limit, repetition_limit):
            continue
        return assignments
    return None


def _generate_sequence_assignments(
    project: TiamatProject,
    rng: random.Random,
    deadline: float,
    unique_limit: int,
    repetition_limit: int,
    min_gc: int,
    max_gc: int,
) -> dict[int, str]:
    strand_order = [strand.strand_id for strand in sorted(project.strands, key=lambda strand: (-len(strand.base_indices), strand.strand_id))]
    strands_by_id = {strand.strand_id: strand for strand in project.strands}
    strand_sequences = {strand.strand_id: [None] * len(strand.base_indices) for strand in project.strands}
    strand_required = {
        strand.strand_id: set(_sequence_alphabet(str(project.bases[strand.base_indices[0]].molecule).upper()))
        for strand in project.strands
        if strand.base_indices
    }
    strand_positions = {strand.strand_id: 0 for strand in project.strands}
    base_to_location = {
        base_index: (strand.strand_id, position)
        for strand in project.strands
        for position, base_index in enumerate(strand.base_indices)
    }
    assignments: dict[int, str] = {}
    total_bases = len(project.bases)

    sys.setrecursionlimit(max(sys.getrecursionlimit(), total_bases * 2 + 1000))

    def propagate_frontiers() -> bool:
        for strand_id in strand_order:
            sequence = strand_sequences[strand_id]
            frontier = strand_positions[strand_id]
            while frontier < len(sequence) and sequence[frontier] is not None:
                if not _position_constraints_ok(
                    sequence=sequence,
                    position=frontier,
                    required_bases=strand_required[strand_id],
                    unique_limit=unique_limit,
                    repetition_limit=repetition_limit,
                ):
                    return False
                frontier += 1
            strand_positions[strand_id] = frontier
        return True

    def choose_next_step(gc_count: int) -> tuple[int, list[tuple[tuple[float, float, float], dict[int, str], int]]] | None:
        chosen: tuple[int, list[tuple[tuple[float, float, float], dict[int, str], int]]] | None = None
        chosen_metric: tuple[int, int] | None = None
        assigned_count = len(assignments)
        for strand_id in strand_order:
            frontier = strand_positions[strand_id]
            sequence = strand_sequences[strand_id]
            if frontier >= len(sequence):
                continue
            strand = strands_by_id[strand_id]
            base_index = strand.base_indices[frontier]
            candidates = _candidate_infos_for_frontier(
                project=project,
                base_index=base_index,
                assignments=assignments,
                strand_sequences=strand_sequences,
                strand_required=strand_required,
                base_to_location=base_to_location,
                current_gc=gc_count,
                assigned_count=assigned_count,
                total_bases=total_bases,
                min_gc=min_gc,
                max_gc=max_gc,
                unique_limit=unique_limit,
                repetition_limit=repetition_limit,
                rng=rng,
            )
            if not candidates:
                return None
            metric = (len(candidates), -(len(sequence) - frontier))
            if chosen is None or metric < chosen_metric:
                chosen = (strand_id, candidates)
                chosen_metric = metric
                if len(candidates) == 1:
                    break
        return chosen

    def backtrack(gc_count: int) -> bool:
        if time.monotonic() >= deadline:
            raise TimeoutError("Could not generate sequences that satisfy the current limits within the timeout.")
        saved_frontiers = dict(strand_positions)
        if not propagate_frontiers():
            strand_positions.update(saved_frontiers)
            return False
        if all(strand_positions[strand_id] >= len(strand_sequences[strand_id]) for strand_id in strand_order):
            if min_gc <= gc_count <= max_gc and _assignments_respect_pairs(project, assignments):
                return True
            strand_positions.update(saved_frontiers)
            return False

        step = choose_next_step(gc_count)
        if step is None:
            strand_positions.update(saved_frontiers)
            return False
        _strand_id, candidates = step
        for _score, changes, new_gc in candidates:
            changed_locations: list[tuple[int, int, int]] = []
            for changed_index, nucleotide in changes.items():
                assignments[changed_index] = nucleotide
                target_strand_id, target_position = base_to_location[changed_index]
                strand_sequences[target_strand_id][target_position] = nucleotide
                changed_locations.append((changed_index, target_strand_id, target_position))
            if backtrack(new_gc):
                return True
            for changed_index, target_strand_id, target_position in reversed(changed_locations):
                assignments.pop(changed_index, None)
                strand_sequences[target_strand_id][target_position] = None
            strand_positions.update(saved_frontiers)
        strand_positions.update(saved_frontiers)
        return False

    if backtrack(0):
        return assignments
    raise TimeoutError("Could not generate sequences that satisfy the current limits within the timeout.")


def _candidate_infos_for_frontier(
    project: TiamatProject,
    base_index: int,
    assignments: dict[int, str],
    strand_sequences: dict[int, list[str | None]],
    strand_required: dict[int, set[str]],
    base_to_location: dict[int, tuple[int, int]],
    current_gc: int,
    assigned_count: int,
    total_bases: int,
    min_gc: int,
    max_gc: int,
    unique_limit: int,
    repetition_limit: int,
    rng: random.Random,
) -> list[tuple[tuple[float, float, float], dict[int, str], int]]:
    base = project.bases[base_index]
    candidates = list(_sequence_alphabet(base.molecule))
    rng.shuffle(candidates)
    results: list[tuple[tuple[float, float, float], dict[int, str], int]] = []
    strand_id, position = base_to_location[base_index]
    sequence = strand_sequences[strand_id]

    for candidate in candidates:
        changes = _assignment_changes_for_candidate(project, base_index, candidate, assignments)
        if changes is None:
            continue
        added_gc = sum(1 for nucleotide in changes.values() if nucleotide in {"G", "C"})
        new_assigned_count = assigned_count + len(changes)
        remaining = total_bases - new_assigned_count
        new_gc = current_gc + added_gc
        if new_gc > max_gc or new_gc + remaining < min_gc:
            continue

        original = sequence[position]
        sequence[position] = candidate
        current_ok = _position_constraints_ok(
            sequence=sequence,
            position=position,
            required_bases=strand_required[strand_id],
            unique_limit=unique_limit,
            repetition_limit=repetition_limit,
        )
        sequence[position] = original
        if not current_ok:
            continue

        gc_target = (min_gc + max_gc) / 2.0
        gc_error = abs(gc_target - new_gc)
        results.append(
            (
                (gc_error, float(len(changes) > 1), rng.random()),
                changes,
                new_gc,
            )
        )

    results.sort(key=lambda item: item[0])
    return results


def _assignment_changes_for_candidate(
    project: TiamatProject,
    base_index: int,
    candidate: str,
    assignments: dict[int, str],
) -> dict[int, str] | None:
    base = project.bases[base_index]
    changes: dict[int, str] = {}
    existing = assignments.get(base_index)
    if existing is not None and existing != candidate:
        return None
    if existing is None:
        changes[base_index] = candidate

    across_index = base.across
    if across_index is None or across_index not in project.bases:
        return changes

    across = project.bases[across_index]
    complement = complement_for(candidate, across.molecule)
    across_existing = assignments.get(across_index)
    if across_existing is not None and across_existing != complement:
        return None
    if across_existing is None:
        changes[across_index] = complement
    return changes


def _assignments_respect_pairs(project: TiamatProject, assignments: dict[int, str]) -> bool:
    for index, base in project.bases.items():
        nucleotide = assignments.get(index)
        if nucleotide is None:
            return False
        if base.across is None or base.across not in project.bases:
            continue
        across = project.bases[base.across]
        across_nucleotide = assignments.get(base.across)
        if across_nucleotide is None:
            return False
        if across_nucleotide != complement_for(nucleotide, across.molecule):
            return False
    return True


def _connected_strand_components(
    strand_ids: list[int],
    pair_constraints: dict[int, dict[int, list[tuple[int, int]]]],
) -> list[set[int]]:
    remaining = set(strand_ids)
    components: list[set[int]] = []
    while remaining:
        start = remaining.pop()
        component = {start}
        stack = [start]
        while stack:
            current = stack.pop()
            for neighbor in pair_constraints.get(current, {}):
                if neighbor in component:
                    continue
                component.add(neighbor)
                remaining.discard(neighbor)
                stack.append(neighbor)
        components.append(component)
    return components


def _solve_periodic_component(
    component: set[int],
    pair_constraints: dict[int, dict[int, list[tuple[int, int]]]],
    states: dict[int, tuple[int, int]],
    rng: random.Random,
) -> bool:
    local_states = {strand_id: states[strand_id] for strand_id in component if strand_id in states}
    if not local_states:
        root = next(iter(component))
        root_states = [(sign, offset) for sign in (1, 3) for offset in range(4)]
        rng.shuffle(root_states)
        for state in root_states:
            local_states[root] = state
            if _backtrack_periodic_component(component, pair_constraints, local_states, rng):
                states.update(local_states)
                return True
            local_states.pop(root, None)
        return False
    if _backtrack_periodic_component(component, pair_constraints, local_states, rng):
        states.update(local_states)
        return True
    return False


def _backtrack_periodic_component(
    component: set[int],
    pair_constraints: dict[int, dict[int, list[tuple[int, int]]]],
    states: dict[int, tuple[int, int]],
    rng: random.Random,
) -> bool:
    if len(states) == len(component):
        return True

    choice: tuple[int, list[tuple[int, int]]] | None = None
    choice_metric: tuple[int, int, int] | None = None
    for strand_id in component:
        if strand_id in states:
            continue
        options: set[tuple[int, int]] | None = None
        constrained_neighbors = 0
        for neighbor, constraints in pair_constraints.get(strand_id, {}).items():
            if neighbor not in states:
                continue
            constrained_neighbors += 1
            neighbor_options = _neighbor_state_options(states[neighbor], pair_constraints[neighbor][strand_id])
            options = neighbor_options if options is None else options & neighbor_options
        if options is None:
            options = {(sign, offset) for sign in (1, 3) for offset in range(4)}
        if not options:
            return False
        ordered = list(options)
        rng.shuffle(ordered)
        metric = (len(ordered), -constrained_neighbors, strand_id)
        if choice is None or choice_metric is None or metric < choice_metric:
            choice = (strand_id, ordered)
            choice_metric = metric

    if choice is None:
        return False
    strand_id, options = choice
    for option in options:
        states[strand_id] = option
        if _backtrack_periodic_component(component, pair_constraints, states, rng):
            return True
        states.pop(strand_id, None)
    return False


def _neighbor_state_options(
    source_state: tuple[int, int],
    constraints: list[tuple[int, int]],
) -> set[tuple[int, int]]:
    source_sign, source_offset = source_state
    options: set[tuple[int, int]] = set()
    first_source_pos, first_target_pos = constraints[0]
    for target_sign in (1, 3):
        target_offset = (source_offset + source_sign * first_source_pos + 2 - target_sign * first_target_pos) % 4
        if all(
            (target_offset + target_sign * target_pos - source_offset - source_sign * source_pos - 2) % 4 == 0
            for source_pos, target_pos in constraints
        ):
            options.add((target_sign, target_offset))
    return options


def _assignments_respect_sequence_constraints(
    project: TiamatProject,
    assignments: dict[int, str],
    unique_limit: int,
    repetition_limit: int,
) -> bool:
    for strand in project.strands:
        sequence = [assignments[index] for index in strand.base_indices]
        required_bases = set(_sequence_alphabet(str(project.bases[strand.base_indices[0]].molecule).upper()))
        if not _repetition_limit_ok(sequence, repetition_limit):
            return False
        if not _unique_window_ok(sequence, unique_limit, required_bases):
            return False
    return True


def _position_constraints_ok(
    sequence: list[str | None],
    position: int,
    required_bases: set[str],
    unique_limit: int,
    repetition_limit: int,
) -> bool:
    current = sequence[position]
    if current is None:
        return False
    run_length = 1
    cursor = position - 1
    while cursor >= 0 and sequence[cursor] == current:
        run_length += 1
        if run_length >= repetition_limit:
            return False
        cursor -= 1
    if unique_limit > 0 and position + 1 >= unique_limit:
        window = sequence[position - unique_limit + 1 : position + 1]
        if any(base is None for base in window):
            return False
        if not required_bases.issubset(set(window)):
            return False
    return True


def _repetition_limit_ok(sequence: list[str], repetition_limit: int) -> bool:
    if repetition_limit <= 1:
        return False
    run_length = 1
    for previous, current in zip(sequence, sequence[1:]):
        if current == previous:
            run_length += 1
            if run_length >= repetition_limit:
                return False
        else:
            run_length = 1
    return True


def _unique_window_ok(sequence: list[str], unique_limit: int, required_bases: set[str]) -> bool:
    if unique_limit <= 0 or len(sequence) < unique_limit:
        return True
    for start in range(0, len(sequence) - unique_limit + 1):
        if not required_bases.issubset(set(sequence[start : start + unique_limit])):
            return False
    return True


def _sequence_alphabet(molecule: str) -> tuple[str, ...]:
    return ("A", "U", "C", "G") if str(molecule).upper() == "RNA" else ("A", "T", "C", "G")


def _strand_position_lookup(strand: Strand) -> dict[int, int]:
    return {index: position for position, index in enumerate(strand.base_indices)}


def _strand_color(strand_id: int) -> str:
    hue = ((strand_id * 0.61803398875) % 1.0)
    r, g, b = colorsys.hsv_to_rgb(hue, 0.55, 0.92)
    return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))


def _sticky_end_color(sticky_id: int) -> str:
    hue = ((sticky_id * 0.38196601125 + 0.17) % 1.0)
    r, g, b = colorsys.hsv_to_rgb(hue, 0.72, 0.9)
    return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))


def _normalize_index_list(values: Any, bases: dict[int, Base]) -> list[int]:
    if not isinstance(values, list):
        return []
    normalized: list[int] = []
    seen: set[int] = set()
    for value in values:
        index = _as_int_or_none(value)
        if index is None or index not in bases or index in seen:
            continue
        normalized.append(index)
        seen.add(index)
    return normalized


def _sticky_end_segments_for_selection(
    project: TiamatProject,
    indices: list[int] | set[int],
) -> tuple[list[int], list[int]] | None:
    selected = {index for index in indices if index in project.bases}
    if not selected:
        return None
    by_strand: dict[int, list[int]] = {}
    strands_by_id = {strand.strand_id: strand for strand in project.strands}
    for index in selected:
        strand_id = project.bases[index].strand_id
        if strand_id is None or strand_id not in strands_by_id:
            return None
        by_strand.setdefault(strand_id, []).append(index)
    if len(by_strand) != 2:
        return None

    segments: list[list[int]] = []
    for strand_id, strand_indices in sorted(by_strand.items()):
        strand = strands_by_id[strand_id]
        positions = [position for position, base_index in enumerate(strand.base_indices) if base_index in strand_indices]
        if not positions:
            return None
        if positions[-1] - positions[0] + 1 != len(positions):
            return None
        segment = strand.base_indices[positions[0] : positions[-1] + 1]
        segments.append(segment)
    if len(segments[0]) != len(segments[1]):
        return None
    return (segments[0], segments[1])


def _sticky_end_sequences_match(
    project: TiamatProject,
    left_indices: list[int],
    right_indices: list[int],
) -> bool:
    for left_index, right_index in zip(left_indices, reversed(right_indices)):
        if not bases_can_pair(project.bases[left_index], project.bases[right_index]):
            return False
    return True


def _choose_nucleotide(history: list[str], usage: Counter[str], rng: random.Random, molecule: str = "DNA") -> str:
    candidates = ["A", "C", "G", "U"] if str(molecule).upper() == "RNA" else ["A", "C", "G", "T"]
    if len(history) >= 2 and history[-1] == history[-2] and history[-1] in candidates:
        candidates.remove(history[-1])
    weights = defaultdict(list)
    for candidate in candidates:
        weights[usage[candidate]].append(candidate)
    lowest = min(weights)
    return rng.choice(weights[lowest])
