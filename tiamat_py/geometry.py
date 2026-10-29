"""Small, dependency-free vector operations shared by structural import/export."""
from __future__ import annotations

import math


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(a, factor):
    return tuple(x * factor for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    length = norm(a)
    if not math.isfinite(length) or length < 1e-10:
        raise ValueError("Cannot determine orientation from coincident or collinear coordinates.")
    return scale(a, 1 / length)


def perpendicular(a, normal):
    return sub(a, scale(normal, dot(a, normal)))


def mean(points):
    return tuple(sum(values) / len(points) for values in zip(*points))
