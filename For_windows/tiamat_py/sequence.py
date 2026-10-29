"""Randomized local repair with paired bases treated as one variable.

This satisfies sequence constraints, not thermodynamic sequence design. It does
not sample uniformly from all valid sequences and may time out on hard inputs.
"""
from __future__ import annotations

import time


def generate_assignments(project, rng, deadline, unique_limit, repetition_limit, min_gc, max_gc):
    locations = {}
    groups = []
    for index in sorted(project.bases):
        if index in locations:
            continue
        group = [index]
        locations[index] = (len(groups), False)
        across = project.bases[index].across
        if across is not None:
            locations[across] = (len(groups), True)
            group.append(across)
        groups.append(group)

    # Encoded A,T/U,C,G: complement is XOR 1, GC membership is >= 2.
    constraints = []
    incidence = [set() for _ in groups]
    for strand in project.strands:
        indices = strand.base_indices
        circular = project.bases[indices[-1]].down == indices[0]
        for width, kind in ((unique_limit, "unique"), (repetition_limit, "repeat")):
            # A full turn suffices when a circular strand is shorter than a window.
            size = min(width, len(indices)) if circular else width
            starts = len(indices) if circular else max(0, len(indices) - size + 1)
            for start in range(starts):
                if time.monotonic() >= deadline:
                    raise TimeoutError("Sequence generation timed out; the original sequence was kept.")
                window = [locations[indices[(start + offset) % len(indices)]] for offset in range(size)]
                constraint_id = len(constraints)
                constraints.append((kind, window))
                for variable, _ in window:
                    incidence[variable].add(constraint_id)

    target = (min_gc + max_gc) / (2 * len(project.bases))
    weights = [len(group) for group in groups]

    def gc_error(count):
        return max(min_gc - count, 0, count - max_gc)

    while time.monotonic() < deadline:
        values = [rng.randrange(2) + (2 if rng.random() < target else 0) for _ in groups]
        gc_count = sum(weight for weight, value in zip(weights, values) if value >= 2)

        def violation(constraint_id):
            kind, window = constraints[constraint_id]
            seen = {values[variable] ^ int(paired) for variable, paired in window}
            return len(seen) < 4 if kind == "unique" else len(seen) == 1

        failed = {i for i in range(len(constraints)) if violation(i)}
        for _step in range(max(1000, 30 * len(groups))):
            if time.monotonic() >= deadline:
                break
            if not failed and gc_error(gc_count) == 0:
                return {index: ("AUCG" if project.bases[index].molecule == "RNA" else "ATCG")[values[var] ^ int(paired)]
                        for index, (var, paired) in locations.items()}
            if failed:
                constraint_id = rng.choice(sorted(failed))
                candidates = sorted({variable for variable, _ in constraints[constraint_id][1]})
                # Bound per-step work for large windows.
                if len(candidates) > 12:
                    candidates = rng.sample(candidates, 12)
            else:
                candidates = [var for var, value in enumerate(values)
                              if (value >= 2) == (gc_count > max_gc)]
                if not candidates:
                    break
                candidates = rng.sample(candidates, min(12, len(candidates)))
            moves = []
            for variable in candidates:
                old = values[variable]
                affected = incidence[variable]
                old_failures = len(failed.intersection(affected))
                for value in range(4):
                    if value == old:
                        continue
                    values[variable] = value
                    next_gc = gc_count + weights[variable] * (int(value >= 2) - int(old >= 2))
                    delta = sum(violation(i) for i in affected) - old_failures + gc_error(next_gc) - gc_error(gc_count)
                    moves.append((delta, variable, value, next_gc))
                values[variable] = old
            if not moves:
                break
            best = min(move[0] for move in moves)
            choices = moves if rng.random() < 0.08 else [move for move in moves if move[0] == best]
            _, variable, value, gc_count = rng.choice(choices)
            values[variable] = value
            for constraint_id in incidence[variable]:
                if violation(constraint_id):
                    failed.add(constraint_id)
                else:
                    failed.discard(constraint_id)
    raise TimeoutError("Could not satisfy the sequence limits within the timeout; the original sequence was kept.")
