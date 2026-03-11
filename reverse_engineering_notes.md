# Reverse engineering notes for `Tiamat_2.exe`

Static findings from `/Users/wyssuser/Documents/Tiamat/Tiamat_2.exe`:

- SHA-256: `a25c3b85e8e94e8dbf175854460825a41af7c74848f61e5029d3b6d6e6373094`
- Size: `2483712` bytes
- Type: `PE32 executable (GUI) Intel 80386, for MS Windows`
- Timestamp in PE header: `2015-08-17 12:43:38`
- Native runtime: Visual C++ / MFC
- No CLR header, so it is not a `.NET` executable
- Graphics stack: `OPENGL32.dll`, `GLU32.dll`, `gdiplus.dll`
- Movie export stack: `AVIFIL32.dll`
- Network stack: `WS2_32.dll`

High-signal UI/resource strings:

- `Tiamat DNA Editor Version Alpha`
- `DNA Files (*.dna)`
- `Tiamat JSON (*.dnajson)`
- `Render Video`
- `Generate Sequences`
- `Create Freeform Strand`
- `View/Edit Constraints`
- `Genome Database`
- `GET /entrez/eutils/esearch.fcgi?tool=Tiamat&db=genome&term=`
- `GET /entrez/eutils/efetch.fcgi?tool=Tiamat&db=genome&retmode=text&rettype=fasta&seq_stop=1&id=`

What this implies:

- The original application is a native MFC/OpenGL CAD-style DNA editor, not a packaged scripting app.
- Exact source recovery is unrealistic from this binary alone.
- A Python port is best approached as behavioral reimplementation.
- The public `dna_ascii.dat` export format is enough to reconstruct a useful base graph:
  - column 1: row index
  - column 2: base object id
  - columns 3-5: symmetric neighbor ids (`up`, `down`, `across`)
  - columns 6-8: `x y z`

The `tiamat_py` package in this repository implements that reimplementation path.

