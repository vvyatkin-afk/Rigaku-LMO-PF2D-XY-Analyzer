# LMO PF2D XY Analyzer

Windows GUI for **Rigaku SmartLabXE / HyPix-3000(H)** 2D pole-figure datasets. Processes one measured position or a batch of XY positions; produces per-point 2D pole figures, four independent peak fits, observed azimuthal FWHM, tilt and beta-cut FWHM, plus overall CSV/JSON and PDF maps.

**No instrument control.** This program only reads RASX files and writes reports to a new folder. It cannot move the XY stage or operate the X-ray source.

## Download / Installation

Download **LMO-PF2D-Analyzer.exe** from the repository's [Releases](../../releases). No separate Python installation is required. Windows 10/11 x64 supported by the build (built on Windows 11).

1. Open the EXE.
2. Select the folder containing Rigaku RASX frames; the program recursively finds matching files.
3. Click **Проверить данные** (Check data) to list series and missing frames.
4. Select **Одна выбранная точка** (one point) or **Все точки** (all points).
5. Select a parent results directory and click **Построить отчёт**.
6. Open results. Every run creates a *new timestamped folder*; originals are never changed.

A successful single-point run creates `LMO011_point_report.pdf`, a measured stereographic pole figure, independent phi fits, peak CSV and geometry logs. A batch additionally creates `LMO011_XY_Overall_Report.pdf`, `XY_summary.csv` and `XY_summary.json`.

## Supported raw acquisition

This initial release deliberately validates a **specific, real, verified acquisition geometry**, rather than silently misinterpreting arbitrary 2D XRD images:

| Parameter | Required |
|---|---|
| Reflection | LMO (011) |
| 2θ | 32.6° |
| ω | 16.3° |
| Detector | HyPix-3000(H), native 775×385 geometry |
| χ settings | 38°, 42° |
| φ | 0°, 2°, 4°, …, 358° for each χ |
| Images per point | 360 = 2×180 |
| Exposure | 5 seconds per frame |
| Output | Observed FWHM, with no instrument-resolution subtraction |

For each point there must be 360 files named as follows (different point folders **or** different prefixes):

```
P001_LMO011_image_1_38.000_001.rasx
P001_LMO011_image_1_38.000_002.rasx
...
P001_LMO011_image_2_42.000_180.rasx
P002_LMO011_image_1_38.000_001.rasx
...
```

The program checks every frame's angle, detector metadata and pixel calibration. A different measurement configuration is rejected, not automatically remapped. Four peaks are independently fitted; missing polar/azimuthal coverage remains unfilled.

## Build from source

Requires Windows x64, Python 3.14 and dependencies in `requirements.txt`.

```powershell
py -3.14 -m pip install -r requirements.txt
build_windows.bat
```

The command creates `dist\LMO-PF2D-Analyzer.exe` using PyInstaller `--onefile --console --hide-console hide-early`. The console is hidden for normal double-click launches but retained if started from a terminal, so automated tests and diagnostic commands can read stdout.

## CLI and validation

The same EXE supports non-interactive checks:

```powershell
LMO-PF2D-Analyzer.exe --inspect --source "D:\Data\LMO_011"
LMO-PF2D-Analyzer.exe --analyze --source "D:\Data\LMO_011" --output "D:\Data\Results_NEW_001"
LMO-PF2D-Analyzer.exe --analyze --source "D:\Data\LMO_011" --group-key ".__011_2D_LMO_011_PF2D_TEST" --output "D:\Data\Results_SINGLE"
LMO-PF2D-Analyzer.exe --smoke-gui "D:\Data\gui_smoke.json"
```

`--output` **must not already exist**. For every output point the program keeps validated geometry, CSV numeric grids, 2D peak fits, processing provenance and the final PDF.

## Validation

Tested with a 360-frame, actual LMO (011) measurement. Peak centers in φ (deg) and observed φ FWHM (deg) were 82.49 / 9.92, 176.55 / 13.58, 271.22 / 9.87, and 356.72 / 8.03, consistent with independently reviewed analysis of that recording. The batch-mode stress test may use two copies of the same measured set; this validates the batching mechanism, **not** two independently measured XY positions.

## Data protection and limitations

- Raw source folders are opened read-only. The app writes to a newly created report folder outside the raw directory.
- No raw laboratory data, cloud credentials, sample photos, internal network paths or proprietary XRD files are included in this repository/release.
- Results depend on actual XRD geometry and calibration. β=0 is the instrument's azimuth reference; it is not automatically tape longitudinal direction.
- Observed φ FWHM includes instrumental, tilt and sampling contributions. This is not an instrument-resolution-corrected in-plane mosaic spread.
- `LMO (200)` and arbitrary χ grids are **not yet supported** by this verified processing pipeline.

## Source map

`lmo_pf2d_gui.py` – Tkinter GUI, quick raw-frame preview and frozen application dispatcher.

`analysis/analyze_many_xy.py` – group validator, per-point PDF assembly, summary plots.

`analysis/validated_single_point_scripts` – validated detector geometry and independent 2D/φ fits; invoked in isolated child processes to keep the GUI responsive.

## v1.1.1 – Flat-folder 39-position Rigaku batches (2026-11-23)

- Handles exported names such as `011_2D_LMO_011_batch_image_1_38.000_001_01.rasx`
  where the *last* number is a separate acquisition/position suffix. Source files
  are never renamed or modified. Grouping uses both this suffix and actual XY
  coordinates recorded in the RASX motor metadata.
- Works with all position series in a **single input folder**. Uses metadata to
  verify the phi grid of each chi series; incomplete positions are clearly marked
  FAILED instead of mixed or interpolated.
- In Auto batch mode, creates `LMO011_XY_Overall_Report.pdf` and
  `XY_summary.csv` with XY coordinates and per-pole beta, tilt and observed
  phi/tilt widths, in addition to the individual pole figure PDFs.
- In the GUI, choose the folder containing the exported RASX files, click
  **Проверить данные**, select **Все точки (пакет)** and **Авто**, specify a
  separate results directory, then click **Построить отчёт**.
- Tested against synthesized naming/position metadata for 39 points; actual
  2026-11-22 sample data were not remotely accessible for validation.

## v1.1.0 – Auto acquisition (2026-11-22)

The new **Auto (step/continuous)** mode reads **actual RASX motor metadata** rather
than assuming a 2° step, 360 frames or only two specified chi values. The
original fixed-grid algorithm remains accessible in the GUI as **Legacy 2°**
and on the CLI as \`--analyze-legacy\`.

\`\`\`powershell
LMO-PF2D-Analyzer.exe --inspect --source "D:\\Data\\LMO011"
LMO-PF2D-Analyzer.exe --analyze --source "D:\\Data\\LMO011" --output "D:\\Data\\LMO011_NEW_RESULTS"
\`\`\`

Adaptive reconstruction uses the recorded Phi, Chi, X, Y, TwoTheta, Omega,
HyPix geometry and actual frame scan windows. It detects 4° and 5°
stationary acquisitions as well as continuous φ sweeps and checks for missing
φ positions without silently inventing frames. Continuous data are integrated
over the actual recorded angle window; stationary images are assigned their
recorded angle. Outputs include acquisition.json, map.npz, peaks.csv/.json
and a measured stereographic PF PDF for each accepted XY group.
The old rigorous 360-frame calibration pipeline is preserved for direct
comparison and independently verified legacy figures.

**Calibration boundary**: This version supports only the validated
HyPix-3000(H) horizontal 775×385 pixel geometry near LMO (011),
2θ≈32.6°, ω≈16.3°, detector distance≈150.6 mm. It rejects inconsistent
or incomplete input sets. Arbitrary detector geometry/phase, changing
optics, in-plane LMO (200), 3D angular inversion, and missing sectors
are not silently extrapolated.

**Offline Flow proposal:** 4° step-and-shoot is **not instrument validated**.
The native vendor Pole Figure(2D) activity must be opened in installed
SmartLab Studio II to check whether 2D(Single exposure) is accepted
and actually holds Phi fixed during every frame. Merely changing the XML
scan-mode text is not proof of executable support. The Flow drafts are
delivered with all hardware actions Skip=True. SaveMappingData=True is
set only in the *offline draft*; on the instrument the operator must
verify its output pathname and no overwrite behavior before enabling.
Do not run unattended; inspect the optics/stage/pinhole collision constraints.
