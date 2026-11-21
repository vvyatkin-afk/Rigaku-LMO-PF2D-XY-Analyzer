@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"
py -3.14 -m PyInstaller --noconfirm --clean --onefile --console --hide-console hide-early --name LMO-PF2D-Analyzer --add-data "analysis;analysis" --hidden-import scipy.optimize --hidden-import scipy.special --hidden-import matplotlib.pyplot --hidden-import matplotlib.backends.backend_pdf --hidden-import numpy --hidden-import pypdf --hidden-import PIL lmo_pf2d_gui.py
if errorlevel 1 (
  echo Build FAILED.
  exit /b 1
)
echo.
echo Built: "%CD%\dist\LMO-PF2D-Analyzer.exe"
for %%F in ("dist\LMO-PF2D-Analyzer.exe") do echo Size: %%~zF bytes
endlocal
