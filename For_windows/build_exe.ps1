py -m pip install -U pyinstaller
py -m PyInstaller --clean TiamatPy.spec
Write-Host ""
Write-Host "Build finished. If successful, check dist\\TiamatPy.exe"
