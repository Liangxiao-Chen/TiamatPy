# Build on Windows with:
#   py -m pip install pyinstaller
#   py -m PyInstaller TiamatPy.spec

from pathlib import Path


project_root = Path(__file__).resolve().parent

datas = [
    (str(project_root / "Toolbar_icons"), "Toolbar_icons"),
    (str(project_root / "demo"), "demo"),
]


a = Analysis(
    [str(project_root / "tiamat_py" / "__main__.py")],
    pathex=[str(project_root)],
    binaries=[],
    datas=datas,
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="TiamatPy",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
