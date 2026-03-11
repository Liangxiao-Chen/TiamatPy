# TiamatPy

`TiamatPy` is a pragmatic Python reimplementation of the core behavior inferred from `Tiamat_2.exe`.

## Download

Clone the repository:

```bash
git clone https://github.com/Liangxiao-Chen/TiamatPy.git
cd TiamatPy
```

Or download the ZIP from GitHub and extract it.

## Installation

`TiamatPy` requires Python 3 with `tkinter`.

Check whether `tkinter` is available:

```bash
python3 -m tkinter
```

If a small Tk window opens, `tkinter` is installed correctly.

### macOS

The easiest option is to install Python from [python.org](https://www.python.org/downloads/macos/). That installer usually includes `tkinter`.

Run:

```bash
python3 -m tiamat_py
```

### Linux

Install Python 3 and `tkinter`:

```bash
sudo apt install python3 python3-tk
```

Then run:

```bash
python3 -m tiamat_py
```

### Windows

Install Python 3 from [python.org](https://www.python.org/downloads/windows/).

During installation:

- enable `Add python.exe to PATH`
- keep the standard `tcl/tk and IDLE` component enabled

Check `tkinter`:

```powershell
python -m tkinter
```

If `python` is not available, try:

```powershell
py -m tkinter
```

Run the app:

```powershell
python -m tiamat_py
```

or

```powershell
py -m tiamat_py
```

### Windows `.exe`

If you want a double-clickable Windows build, use the self-contained build folder:

```powershell
cd For_windows
py -m PyInstaller --clean TiamatPy.spec
```

or double-click:

- `For_windows\build_exe.bat`

The packaged application will be created at:

```powershell
For_windows\dist\TiamatPy\TiamatPy.exe
```

To share it with other Windows users:

- send the whole `For_windows\dist\TiamatPy\` folder
- do not send only `TiamatPy.exe`
- the other user can run `TiamatPy.exe` by double-clicking it
- Python is not required on the target Windows machine

What is implemented:

- `.dna`, `.dnajson`, and nucleic-acid `.pdb` loading
- JSON save/load (`.dnajson`)
- Strand reconstruction from the base graph
- Heuristic sequence generation with Watson-Crick pairing across links
- Sequence export as FASTA-like text
- Static SVG rendering
- A four-panel Tk viewer with XY/YZ/XZ/3D views, zoom, pan, rotation in 3D, selection, base editing, and export actions
- NCBI genome search and FASTA download

What is not implemented yet:

- AVI rendering
- Exact parity with the original MFC/OpenGL UI

Demo files:

- `demo/test.dna`
- `demo/test.dnajson`

The two demo files represent the same structure in native and JSON forms.

Open the demo:

```bash
python3 -m tiamat_py gui demo/test.dna
python3 -m tiamat_py gui demo/test.dnajson
```

On Windows:

```powershell
python -m tiamat_py gui demo\test.dna
python -m tiamat_py gui demo\test.dnajson
```

If you built the packaged Windows app:

```powershell
For_windows\dist\TiamatPy\TiamatPy.exe
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
python3 -m tiamat_py gui /path/to/project.dna
python3 -m tiamat_py gui /path/to/project.dnajson
python3 -m tiamat_py gui /path/to/project.pdb
python3 -m tiamat_py export-sequences /path/to/project.dnajson /path/to/sequences.txt --generate
```

On Windows, use `python` or `py` and Windows-style paths, for example:

```powershell
python -m tiamat_py gui C:\path\to\project.dna
python -m tiamat_py gui C:\path\to\project.dnajson
python -m tiamat_py gui C:\path\to\project.pdb
```

You can also start the viewer without an input file:

```bash
python3 -m tiamat_py
```
