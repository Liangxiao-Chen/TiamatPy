@echo off
setlocal
py -m pip install -U pyinstaller
py -m PyInstaller --clean TiamatPy.spec
echo.
echo Build finished. If successful, check dist\TiamatPy.exe
pause
