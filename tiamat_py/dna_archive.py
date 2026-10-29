"""Tiamat schema-5 MFC archive reader.

Objects are serialized recursively, with across/down/up references BEFORE their
coordinate payload. MFC object IDs are not nucleotide IDs: the final root list
defines the latter. An explicit stack avoids recursion limits on large designs.
Unknown schemas are rejected rather than reconstructed from geometry.
"""
from __future__ import annotations

import struct

from .model import Base, TiamatProject

TYPE_CODES = {0: "A", 1: "T", 2: "U", 3: "C", 4: "G", 5: None}


def decode_archive(data: bytes):
    position = 0

    def take(size):
        nonlocal position
        if position + size > len(data):
            raise ValueError(f"Truncated .dna archive at byte {position}.")
        result = data[position:position + size]
        position += size
        return result

    def integer(size=2):
        return int.from_bytes(take(size), "little")

    def new_class(expected):
        if integer() != 0xffff or integer() != 5:
            raise ValueError("Unsupported .dna archive schema (expected schema 5).")
        name = take(integer())
        if name != expected:
            raise ValueError(f"Unsupported .dna class {name!r}; expected {expected!r}.")

    for name, size in ((b"CTop", 12), (b"CPerspective", 92), (b"CSide", 12), (b"CFront", 12)):
        new_class(name)
        take(size)
    count = integer(4)
    if count == 0:
        if position != len(data):
            raise ValueError("Unexpected data in empty .dna archive.")
        return TiamatProject({}), []
    if count > len(data) // 59:
        raise ValueError("Invalid nucleotide count in .dna archive.")
    new_class(b"Nucleobase")
    objects = {}
    roots = []
    next_id = 10  # Four viewport classes + objects, then Nucleobase class 9.

    def new_object():
        nonlocal next_id
        if len(objects) >= count:
            raise ValueError("More nucleotide objects than declared in .dna archive.")
        object_id = next_id
        next_id += 1
        objects[object_id] = {"refs": []}
        return object_id

    def reference():
        tag = integer()
        extended = tag == 0x7fff
        if extended:
            tag = integer(4)
        if tag == (0x80000009 if extended else 0x8009):
            return new_object(), True
        if tag not in objects:
            raise ValueError(f"Unknown object reference {tag} at byte {position}.")
        return tag, False

    first = new_object()
    for root_index in range(count):
        root, fresh = (first, True) if root_index == 0 else reference()
        roots.append(root)
        stack = [root] if fresh else []
        while stack:
            obj = objects[stack[-1]]
            if len(obj["refs"]) < 3:
                present = integer(1)
                if present == 0:
                    obj["refs"].append(None)
                elif present == 1:
                    linked, new = reference()
                    obj["refs"].append(linked)
                    if new:
                        stack.append(linked)
                else:
                    raise ValueError(f"Invalid connection flag at byte {position - 1}.")
            else:
                obj["offset"] = position
                payload = take(56)
                x, y, z, type_code, flag = struct.unpack_from("<dddII", payload)
                if type_code not in TYPE_CODES:
                    raise ValueError(f"Unsupported nucleotide type {type_code} at byte {position - 32}.")
                molecule = "RNA" if payload[51] == 1 or type_code == 2 else "DNA"
                obj.update(position=(x, y, z), nucleotide=TYPE_CODES[type_code], molecule=molecule, flag=flag)
                stack.pop()
    if len(set(roots)) != count or len(objects) != count:
        raise ValueError("Duplicate or missing nucleotide in .dna root list.")
    if position != len(data):
        raise ValueError(f"Unsupported trailing data in .dna archive at byte {position}.")
    indices = {object_id: index for index, object_id in enumerate(roots)}
    bases = {}
    for object_id, index in indices.items():
        obj = objects[object_id]
        across, down, up = (indices.get(ref) for ref in obj["refs"])
        bases[index] = Base(index, str(index), *obj["position"], molecule=obj["molecule"],
                            nucleotide=obj["nucleotide"], across=across, down=down, up=up)
    project = TiamatProject(bases, metadata={"format": "dna_mfc_schema5", "dna_topology": "archive_references"})
    return project, sorted(objects.values(), key=lambda item: item["offset"])


def encode_archive(project: TiamatProject) -> bytes:
    """Write the supported schema with explicit topology, never a sidecar dependency."""
    project.validate()
    result = bytearray()

    def new_class(name):
        result.extend(struct.pack("<HHH", 0xffff, 5, len(name)) + name)

    new_class(b"CTop")
    result.extend(struct.pack("<fff", 10, 0, 0))
    new_class(b"CPerspective")
    result.extend(struct.pack("<f11d", 10, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1))
    for name in (b"CSide", b"CFront"):
        new_class(name)
        result.extend(struct.pack("<fff", 5, 0, 0))
    result.extend(struct.pack("<I", len(project.bases)))
    seen = {}

    def reference(index):
        if index in seen:
            tag = seen[index]
            result.extend(struct.pack("<H", tag) if tag < 0x7fff else struct.pack("<HI", 0x7fff, tag))
            return False
        if not seen:
            new_class(b"Nucleobase")
        else:
            result.extend(struct.pack("<H", 0x8009))
        seen[index] = 10 + len(seen)
        return True

    type_codes = {value: key for key, value in TYPE_CODES.items()}
    for index in sorted(project.bases):
        stack = [(index, 0)] if reference(index) else []
        while stack:
            index, slot = stack[-1]
            base = project.bases[index]
            if slot < 3:
                neighbor = (base.across, base.down, base.up)[slot]
                stack[-1] = (index, slot + 1)
                result.append(int(neighbor is not None))
                if neighbor is not None and reference(neighbor):
                    stack.append((neighbor, 0))
            else:
                payload = bytearray(struct.pack("<dddII", *base.position, type_codes[base.nucleotide], 0) + b"\x00" * 24)
                payload[51 if base.molecule == "RNA" else 52] = 1
                result.extend(payload)
                stack.pop()
    return bytes(result)
