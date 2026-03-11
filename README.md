# TiamatPy

`TiamatPy` is a pragmatic Python reimplementation of the core behavior inferred from `Tiamat_2.exe`.

What is implemented:

- Tiamat ASCII `.dat` import
- JSON save/load (`.dnajson`)
- Strand reconstruction from the base graph
- Heuristic sequence generation with Watson-Crick pairing across links
- Sequence export as FASTA-like text
- Static SVG rendering
- A four-panel Tk viewer with XY/YZ/XZ/3D views, zoom, pan, rotation in 3D, selection, base editing, and export actions
- NCBI genome search and FASTA download

What is not implemented yet:

- Native `.dna` binary parsing
- AVI rendering
- Full constraint editing
- Exact parity with the original MFC/OpenGL UI

Demo project:

```bash
python3 -m tiamat_py gui /Users/wyssuser/Documents/Tiamat_py/demo/demo_double_helix.dnajson
```

Viewer controls:

- `Command+A` or `Ctrl+A`: select all bases
- Left-drag in `3D`: rotate
- Right-drag in any view: pan
- Mouse wheel in any view: zoom
- Click: select nearest base
- `A` / `C` / `G` / `T`: set selected bases
- `Delete` / `Backspace`: clear selected bases

Typical commands:

```bash
python3 -m tiamat_py inspect /path/to/dna_ascii.dat
python3 -m tiamat_py convert /path/to/dna_ascii.dat /path/to/project.dnajson
python3 -m tiamat_py export-sequences /path/to/project.dnajson /path/to/sequences.txt --generate
python3 -m tiamat_py export-svg /path/to/project.dnajson /path/to/view.svg
python3 -m tiamat_py gui /path/to/project.dnajson
```

You can also start the viewer without an input file:

```bash
python3 -m tiamat_py
```
