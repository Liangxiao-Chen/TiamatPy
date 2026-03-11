# Windows Build

This folder is a self-contained Windows build copy of the project.

From this folder, build the `.exe` with either:

```powershell
py -m PyInstaller --clean TiamatPy.spec
```

or double-click:

- `build_exe.bat`

The built application should appear at:

```powershell
dist\TiamatPy\TiamatPy.exe
```

For distribution, share the whole `dist\TiamatPy\` folder, not only `TiamatPy.exe`.
