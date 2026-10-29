# Structural Correctness and Regression Checks

## October 2026 Audit Fixes

1. **ID 0:** zero is a valid nucleotide ID in JSON, native archives, undo and redo.
   Null references are `null`, not zero. ASCII null sentinels stay in the ASCII importer.
2. **Save As:** the current document path is runtime state, separate from history.
   Reopening a renamed/copied JSON uses its actual path, not an embedded old path.
3. **Data protection:** New, Open, recent-file Open, Exit and window close prompt
   to save/discard/cancel. Failed/cancelled saving does not discard edits. JSON is
   staged in the destination directory, flushed and atomically replaced. This is
   not autosave or crash recovery. Two-file exports roll back ordinary I/O failures;
   a power loss during replacement of two files is not an atomic transaction.
4. **oxDNA frames:** the empirical fake-neighbor conversion was replaced by an
   independent local-frame reconstruction. Frames are orthonormal, and rigid
   transformations commute with conversion. Missing partners are never replaced
   with world-space dummy coordinates. Drag rotation also uses original frames
   rather than accumulating rotations on already-rotated frames.
5. **Native archives:** schema-5 MFC object references are decoded with an explicit
   stack. Across/down/up are read from references; the root array determines IDs.
   Tests cover DNA, RNA, generic RNA, circular strands, origin coordinates, deep
   graphs, moved files and unrelated JSON neighbors. Old TiamatPy position-only
   exports are explicitly reported as lacking topology, not silently connected.
6. **PDB pairing:** geometric, mutually unambiguous Watson-Crick/G-U candidates
   replace sequence-only chain matching. First-model selection, residue-level
   alternate conformers, TER boundaries and O3'-P continuity are respected.
7. **Sequences:** random initialization and local repair replace the period-four
   shortcut. A paired pair is a single variable, so complementarity is always
   maintained. Sliding diversity and repetition windows include circular joins.
   GC fraction is global across all bases. A fixed seed is reproducible; timeout
   leaves the input sequence unchanged. Sampling is not uniform and does not
   optimize melting temperature, secondary structure or off-target interactions.
8. **Graph validation:** duplicate IDs, dangling or nonreciprocal links, self-links,
   unknown molecule types and non-finite coordinates are rejected before traversal.
   Export additionally checks complementary pairs. Invalid inputs produce errors,
   not silent repair or unbounded traversal.

## oxDNA Coordinate Contract

Editor nucleotide positions, including imported PDB pentose centers, are treated
as representative **backbone sites**, not particle centers of mass. This is a
coarse-graining assumption, not an atomistic parameter fit.

One oxDNA length unit is **0.8518 nm**. With orthonormal `a1`, `a3` and
`a2 = a3 x a1`, the exported center is calculated by inverting the site offsets:

| Model | Backbone site relative to center of mass |
| --- | --- |
| DNA2 | `-0.34 a1 + 0.3408 a2` |
| RNA | `-0.4 a1 + 0.2 a3` |

Thus reconstructing backbone sites from the export recovers the input positions,
including their asymmetric groove geometry. **DNA1 uses different offsets** and
will not display the same result. A drawing with incorrect grooves cannot be
fixed simply by exporting it. Pair planes are estimated from stored frames or
local pair/backbone geometry; unpaired tails inherit neighboring orientations.
An isolated point or straight unpaired strand without a stored frame has an
underdetermined orientation and is rejected, rather than assigned a world axis.

The classic topology format is emitted in **3' to 5'** row order, with each row's
neighbors ordered **3', 5'** (`down`, `up`). Unknown bases require an explicit
fallback in the CLI; random fallback respects pairing without modifying the project.
The box grows with the geometry instead of being capped at 250 units.

Orthogonality and backbone reconstruction are regression-tested, but neither is
proof of low energy. Custom helix dimensions, stretched connections, clashes,
PDB coarse-graining and unusual pair geometry can require substantial relaxation.
Production molecular dynamics and a new Windows executable were not validated
by this code audit.

Model references: [oxDNA format and units](https://lorenzo-rovigatti.github.io/oxDNA/configurations.html),
[DNA site definitions](https://github.com/lorenzo-rovigatti/oxDNA/blob/master/src/Particles/DNANucleotide.cpp),
[RNA site definitions](https://github.com/lorenzo-rovigatti/oxDNA/blob/master/src/Particles/RNANucleotide.h),
[RNA parameters](https://github.com/lorenzo-rovigatti/oxDNA/blob/master/src/Interactions/rna_model.h).
The replacement exporter is independently implemented and does not contain the
previous tacoxDNA empirical conversion routines.

## Tests

`python -m unittest discover -s tests -v` runs without creating a graphical window.
Portable tests use the already-published demo and synthetic structures. Optional
local PDB/DNA fixtures add coverage without requiring or publishing private files.
The same source and tests are mirrored into `For_windows`; CI checks both trees.

Mac GUI startup and document editing/save/undo smoke checks are separate from the
headless tests. GUI prompts are unit-tested with controlled dialog responses.
Unknown archive versions, modified PDB residues, and unsupported mixed-molecule
simulation exports remain explicit limitations.
