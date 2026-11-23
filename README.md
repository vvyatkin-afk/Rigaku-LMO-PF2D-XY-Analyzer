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

## v1.3.0 – One-click analysis, live point progress, English GUI (2026-11-23)

All application controls, dialogs, progress messages, and logs are in English.
For normal operation, **there is no separate metadata-check button**.

**One-click processing:** choose a folder containing RASX files (or its
unique parent), set the output parent, and press **Analyze & Generate Reports**.
The program discovers the RASX files **once**, parses the acquisition metadata,
validates each position's geometry, reconstructs the measured 2D pole figures,
fits all four poles, and produces the individual and batch PDF/CSV outputs.
The table of positions populates automatically from metadata during this
same job. No prior **Check data** pass is needed.

Progress is shown every few seconds: discovered and accepted frame counts,
which position is being processed, reconstructed frame counts within each
position, most recently reconstructed RASX filename, peak fit/report stage,
completion status, and elapsed time. Long network reads do not freeze the UI.

The actual image binary for each frame must still be read during
reconstruction; acquisition metadata and image intensity are distinct data.
The optimization removes the *extra preliminary metadata pass*, not the
scientifically necessary reading and fitting operations.

**One selected point** can now be requested directly by entering its numeric
acquisition index (e.g., 37) without first scanning every point via a
separate button. A metadata discovery is still performed inside the same
job to determine the proper XY/chi grouping.

**Fast report rebuild:** button **Rebuild PDFs from Saved Maps** uses the
previously generated `inspection.json` and `point_XXX/map.npz` arrays.
It refits/report-exports the existing reconstructed maps without reading
the raw network RASX data again.

**Disk safety:** all outputs are created in a new timestamped results
directory; the user-specified output parent may not yet exist and will be
created automatically. Old batch results and source RASX are never changed.

The 2D color map, four-pole fits, observed width definitions, and
conditional tilt model are unchanged from v1.2.1. This is a GUI/pipeline
orchestration release with regression tests for one metadata discovery,
English-language UI, and live worker events.

## v1.2.1 – RASX folder selection, output creation, reports-only CLI fix (2026-11-23)

- **Critical bug fix**: the v1.2.0 command-line worker accessed
  `args.reports_only` without registering the flag in `argparse`.
  That caused `AttributeError: Namespace has no attribute reports_only`
  in **both** normal and cached-map report generation, before any input
  could be processed. The option now exists and has integration tests that
  run the actual command-line dispatcher to produce four-page PDF reports.
- **Input folder auto-selection**: Choose a parent folder (up to 3 nested
  levels) containing a single directory with RASX files. The GUI resolves
  this to the exact RASX-containing folder and displays the selected path.
  If multiple subfolders contain different RASX scans, the program
  refuses to guess and requests a more specific folder.
- **Cached reports folder auto-selection**: Button
  `Новые PDF из готового расчёта` searches up to 2 nested levels for
  a unique prior `inspection.json` and saved `point_XXX/map.npz`.
  No RASX read/reconstruction is necessary; prior folders are preserved.
- **Nonexistent output folder accepted**: You may type a new output parent
  path manually; the GUI creates it recursively and adds a separately
  timestamped `LMO_PF2D_BATCH_...` or `LMO_PF2D_REPORTS_...` child.
  If the path cannot be created, it displays the actual permission or
  filesystem error, rather than failing silently.
- **Version in EXE filename**: `LMO-PF2D-Analyzer-v1.2.1.exe`
  and in the main window. CI also smoke-tests the built EXE
  with `--analyze --help` before publishing it.

### Processing existing 39-point data without rereading RASX

1. Run the **v1.2.1 EXE**.
2. Set the output parent (it need not exist; it will be created).
3. Click **Новые PDF из готового расчёта**.
4. Select the already calculated
   `LMO_PF2D_BATCH_20261123_034719` folder, or its parent if it contains
   exactly one saved batch directory.
5. Review the new `LMO_PF2D_REPORTS_...` directory, point PDFs,
   `XY_summary.csv`, and the comprehensive batch PDF.

This version keeps the same 2D reconstruction and color palette as v1.2.0.
It does not overwrite raw RASX files or earlier results.

## GitHub storage policy (2026-11-23)

Windows builds run only for a new `v*` release tag or manual dispatch,
not automatically on each `main` update. GitHub Actions does **not** upload
an EXE artifact or save a pip cache. The downloadable Windows EXE is attached
to the latest GitHub Release, with old binary Releases removed when newer
versions are published. Source files and scientific tests stay versioned.

## v1.2.0 – Detailed four-pole reports and cached-map report rebuild (2026-11-23)

This version addresses batch results from Ivan's 39-point LMO (011) measurements.
The original `magma` pole figure palette is preserved.

**Per-point PDF (4 pages)**: measured stereographic pole figure with fitted
P1-P4 centers; detailed numeric parameter table including four azimuth and
radial centers (`beta`, `t`), stereographic radius `r=tan(t/2)`,
four observed FWHM(beta) and FWHM(t), and fit R2; eight measured/model cut
charts (beta and radial t for each pole) using exactly the fitted 2D Gaussian
parameters; four local 2D intensity maps with half-max fitted ellipses.
The derived `peak_profile_cuts.csv`, `peaks.csv`, `peaks.json` and
`structure_tilt.json` are saved beside each PDF.

**Batch PDF and CSV**: XY maps for all four radial centers, all four
azimuthal and radial widths and stereographic radii; all four width trends
by measured position, a conditional structure-axis tilt map/direction;
full-position numeric tables for beta centers, radial centers, both widths,
stereographic radius, R2 and cone-fit tilt/direction/RMS. `XY_summary.csv`
has **numeric** position IDs (`1`, `2`, …, `39`) rather than leading-zero
text (`01`, `02`, …). All four values are exported separately.

**Fast re-reporting without rescanning a gigabyte of RASX**:

1. Open the updated Windows EXE (close the previous version first).
2. Click **Новые PDF из готового расчёта**.
3. Select the earlier `LMO_PF2D_BATCH_...` folder that contains
   `inspection.json` and `point_001/map.npz` etc.
4. Set a separate parent folder for outputs, if necessary. The program
   creates a **new** `LMO_PF2D_REPORTS_...` folder.
5. Open `LMO011_XY_Overall_Report.pdf`; the per-point PDFs are under
   `point_001/ ... point_039/`.

No raw RASX are read or modified in this path. The 2D Gaussian parameters
are refitted on the **stored measured maps** in the prior batch. This means
each fitted cut shown in the PDF corresponds to the currently exported
widths, not to a different independent 1D fit. The source batch result
folder also remains read-only.

CLI:

```powershell
LMO-PF2D-Analyzer.exe --analyze --reports-only --source "D:\Results\LMO_PF2D_BATCH_20261123_034719" --output "D:\Results\LMO_PF2D_REPORTS_NEW"
```

**Scientific interpretation**: The optional cone-axis tilt follows the
best-fit normal to a plane of four scattering-vector directions computed
from the four measured beta and t centers. It assumes the four peaks
belong to an equivalent cone. Its direction is in *instrument beta*, not
the tape longitudinal coordinate system; it is **not** an instrument
alignment-corrected tape tilt. RMS of the cone model is reported for QC.
All observed FWHM include detector/motor broadening; no instrumental
deconvolution is claimed. Incomplete data are explicitly marked FAILED.

## v1.1.2 – Responsive 5-second scan progress (2026-11-23)

The Windows GUI's **Проверить данные** button now reports progress at roughly
five-second intervals while checking large remote RASX folders, including:
- Stage: initial setup, finding files, sorting, reading RASX motor/header
  metadata, grouping by XY and filename suffix, validating chi/phi series
- Discovered and processed file counts (e.g., `1450/7020`)
- Header files successfully read, failures, the file being read, and the
  last successfully read filename
- Completed/valid group counts for the final per-position validation
- Elapsed time even when a single slow SMB file read temporarily stalls;
  all feedback is generated on the GUI thread rather than a network worker

The source RASX files remain read-only. While checking, the **Построить отчёт**
button is disabled until the scan finishes. This update does not change
the scientific reconstruction/calibration algorithm or accelerate remote SMB
storage; it makes slow scans observable.

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
