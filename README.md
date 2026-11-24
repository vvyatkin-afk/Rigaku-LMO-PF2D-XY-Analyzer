# v1.5.1 - Universal measured radial diffraction and lattice strain

One Windows application for Rigaku SmartLabXE / horizontal HyPix pole-figure RASX data. This release retains the published v1.4.1 pipeline and incorporates the v1.5.0 Deniz coverage changes. The old fixed-geometry worker and old report rebuild remain available. All new UI and reports are English; raw acquisitions are read-only.

## Start here

1. Select the folder containing raw 2D Phi-scan RASX files (including subfolders).
2. Set an output parent (it may not exist yet).
3. Optionally choose **Reference JSON** or use **Reference & Calibration...** to edit and save a JSON profile. The default LMO pseudocubic 3.880 A and wavelength 1.5406 A are only nominal, unconfirmed assumptions.
4. Select one group or all measured groups; click **Analyze & Generate Reports**.
5. See the individual PDF and CSV in each sample folder and the overall PDF, batch CSV and combined-point JSON in the new timestamped result directory.

The same program accepts 011, 111, 222 and other explicitly named/assigned hkl, other reference cells, arbitrary recorded Phi steps and supported detector motor settings. A filename or angle match is never a phase identification. Other mounting conventions, 1D scans, and non-Phi acquisition modes still fail explicitly. In-plane/axial or incomplete pole coverage may allow aggregate d-spacing but not a full indexed cell. Geometry outside the previously verified 011 setup requires experimental calibration.

## Radial analysis

The per-pixel absolute scattering angle is reconstructed from the recorded frame geometry. The original exposure-normalized, sideband-subtracted pole map is retained. In that same image read, counts and pixel coverage are accumulated into a (tilt, beta, 2theta) cube, with 0.025 deg radial bins over nominal 2theta +/-2.4 deg (independent of the original +/-1.2 deg pole-map ROI). Motor angles alone are not fitted peak positions. Peaks outside this window cannot be recovered by cached-map re-reporting.

For each resolved pole, a central angular aperture selects the radial profile; the observed intensity is averaged using actual sampled-pixel coverage. A Gaussian with a free linear background is fitted. Optional fixed-ratio equal-width wavelength doublet fitting is available if the optics justify it. Missing bins remain missing. Unresolved, clipped, gapped or poor profiles are rejected; the aggregate fallback is clearly marked rank 0 and must never be used as an indexed reciprocal vector.

The measured angle and configured or header-derived wavelength determine d=lambda/[2 sin(2theta/2)]. Each pole has a fitted direction and reciprocal vector g=(s_out-s_in)/lambda, in inverse Angstrom, without 2*pi. The vector azimuth corrects a historical half-bin offset in the stored pole-map convention; original beta reports remain unchanged for comparison. The direction still comes from the fitted integrated pole map rather than a full joint 3D peak fit. Phi sweeps use the recorded angular interval. Angular and radial sampling, overlapping phases, spatial displacement and detector calibration can bias the result.

## Indexed cell and deformation

The solver uses signed, correctly indexed reciprocal vectors H and G in the same single domain. In row notation, G=H B^T; the fitted reciprocal basis B gives direct basis A=B^-T. Reject rank-deficient, ill-conditioned, left-handed or inconsistent solutions. A rank-three solution using exactly three vectors has no redundancy. With a reference direct basis A0, F=A A0^-1. SVD gives the rotation and symmetric left stretch V, with strain represented as V-I in specimen coordinates. Pure rigid rotation does not create strain. Cell lengths, angles, volume, relative length/volume changes, principal strains/directions and specimen-normal/tape-axis projections are exported when determinable. Tape directions require a separately calibrated azimuth. This is a lattice-deformation model, not a mechanical stress or elastic-modulus measurement.

There are three indexing routes:

- Explicit signed hkl_by_rank entries for a known exact sample/acquisition. Verify ranks against that run; intensity rank is not a stable crystallographic identity across points.
- For a confirmed single-domain cubic texture, set cube_texture_azimuth_deg to the crystal a-axis azimuth and assume_single_domain=true. The code checks nearest candidate uniqueness and pole coverage.
- Default exploratory mode may infer a **hypothetical** four-pole cubic texture model for {011}, {111}, {222}, only if four resolved and covered vectors match the expected symmetry. The result is explicitly hypothetical and a/b labels are symmetry-ambiguous; this is not an independently indexed, calibrated cell. Disable allow_hypothetical_four_pole_cell to suppress it. Hypothetical indexing is never used in a cross-reflection combined solution.

If any required information is unavailable, the program reports N/A and a status, rather than inventing cell parameters. A single 002 or 200 d-spacing does not determine the general 3D cell.

The reference d-spacing is computed from the configurable reference cell metric. For noncubic cells, a family hint alone does not specify which permutation/sign member was measured; explicit indexing is needed for directional reference strain. Reference, wavelength and calibration confirmation flags are user declarations. Fit sigma is model-fit uncertainty, not absolute instrument accuracy. Optionally supply two_theta_systematic_sigma_deg and angular_sigma_deg for a conditional Monte Carlo uncertainty calculation; phase, wavelength, reference and sample-position systematic errors are not included.

See [example reference JSON](examples/reference_lmo.json). If use_header_wavelength is true and the manual wavelength is not confirmed, an explicitly recorded Ka/Ka1/Ka2/Kb line supplies its nominal wavelength. Generic Ka does not verify monochromatization. All lengths are Angstrom unless explicitly in mm, all angles are degrees, and detector center offsets are pixels. two_theta_offset_deg is added to the raw measured radial center; detector_center_offset_px and detector_distance_offset_mm change the actual reconstruction and require reprocessing raw files. Do not use these offsets to force agreement with a nominal reference.

For cross-reflection cell reconstruction, explicit indexing rules can assign the same physical_point_id to the same actual specimen point and domain. Exact source-folder/sample/position suffix matching prevents accidental merging. Stage XY coordinates must also agree; identical XY values alone never justify joining different samples. Example rule:

~~~json
{
  "sample": "LMO_011",
  "position_suffix": "01",
  "physical_point_id": "verified-specimen-A-position-01",
  "family_hkl": [0,1,1],
  "hkl_by_rank": {"1":[1,0,1], "2":[0,1,1], "3":[-1,0,1], "4":[0,-1,1]}
}
~~~

The example ranks are illustrative only, not a verified assignment for any raw measurement.

## Reports and outputs

Existing measured PF, acquisition coverage, peak tables, local maps and batch comparison pages remain. New pages include fitted pole cuts, conditional cone tilt, radial profiles, d-spacing and fit QC, cell/tensor tables, 3D cell diagrams, XY scatter maps (zero-centered signed scales), trends and histograms. No interpolation of unmeasured areas. All peak and cell data are machine-readable. Unassigned pole ranks are not crystallographically tracked between XY points.

Per-point: Radial_Cube.npz, Radial_Profiles.npz, Lattice_Analysis.json, Reflection_dSpacing.csv, Lattice_Parameters.csv, Lattice_Strain.csv, Strain_Tensor.json, Lattice_Settings.json, peak_profile_cuts.csv, structure_tilt.json, plus previous files. Batch: Batch_Lattice_Summary.csv, Reflection_dSpacing.csv, Lattice_Parameters.csv, Lattice_Strain.csv, Strain_Tensor.json, Combined_Physical_Points.json, Lattice_Settings.json, and the original summary/export files.

Cached **Rebuild PDFs from Saved Maps** can refit profiles and use a changed wavelength/reference/2theta offset from Radial_Cube.npz without rereading RASX. Old integrated maps have no radial axis; their lattice status is RAW_RADIAL_DATA_REQUIRED. A changed detector center/distance must be reconstructed from raw files. Legacy v1.3 integrated-map reports still use the original formatter; process raw RASX for new strain information.

~~~powershell
LMO-PF2D-Analyzer-v1.5.1.exe --analyze --source "D:\Raw\PF2D" --output "D:\Results\NEW" --reference-config "D:\Reference.json"
LMO-PF2D-Analyzer-v1.5.1.exe --analyze --reports-only --source "D:\Results\NEW" --output "D:\Results\REFIT" --reference-config "D:\Reference2.json"
~~~

## Scientific references and validation boundary

- IUCr, [The reciprocal lattice](https://www.iucr.org/what-we-do/education/pamphlets/reciprocal-lattice): scattering vector and Bragg relation.
- [Improved calibration of area detectors using multiple placements](https://journals.iucr.org/s/issues/2019/06/00/co5127/): beam center, distance, pixel dimensions, detector tilt and wavelength are needed for reliable angular calibration.
- [pyFAI calibration documentation](https://www.silx.org/doc/pyFAI/dev/calibration.html): calibrant-based geometry refinement and independent visual verification.

Synthetic tests check exact lattice inversion including pure rotation, isotropic strain down to 0.01% in the mathematical model, anisotropic and shear strain, rank failure, radial profiles, doublets, clipped/gapped data, cached re-reporting and actual RASX decoding. A mathematical test at 0.01% does **not** establish experimental sensitivity to 0.01%. Validate against a calibrated standard and independent theta-2theta measurements before using these outputs as acceptance criteria. Actual raw laboratory data and internal paths are not published.

## Historical documentation follows

The older version descriptions below remain for reproducibility. Their fixed geometry limitations and screenshots refer to those releases, not the v1.5.1 universal worker.

# v1.5.0 - LMO 011 / 111 / 222

Universal RASX analysis based on actual motor/detector metadata. New explicit PF111/PF222 recognition, nominal LMO angle candidates, per-chi angular coverage, and fit coverage QC. Existing v1.4.0/1.4.1 maps can be rebuilt. Filename and Bragg-angle hints are not phase identification. Real 111/222 measurement validation is pending; raw files are read-only. Nominal reference assumes pseudocubic LMO a=3.880 A, Cu K-alpha 1.5406 A.

# LMO PF2D XY Analyzer

Windows GUI for **Rigaku SmartLabXE / HyPix-3000(H)** 2D pole-figure datasets. Processes one measured position or a batch of XY positions; produces per-point 2D pole figures, four independent peak fits, observed azimuthal FWHM, tilt and beta-cut FWHM, plus overall CSV/JSON and PDF maps.

**No instrument control.** This program only reads RASX files and writes reports to a new folder. It cannot move the XY stage or operate the X-ray source.

## v1.4.1 - Visual report quality patch (2026-11-23)

Supersedes v1.4.0 for new work. Scientific reconstruction unchanged.
First-page map/colorbar/table geometry is revised to prevent overlap, and
all report strings use portable ASCII angle units ("deg") to eliminate
Windows character encoding artifacts. The peak resolution-quality guard
and the multi-sample report pipeline are retained. Cached v1.4.0 maps
can be regenerated into v1.4.1 PDFs without re-reading raw frames.
The Windows executable and GitHub tag include v1.4.1 in their filenames.
The GitHub build workflow retains the newest two versioned releases.

## v1.4.0 - Universal acquisition settings and sample reports (2026-11-23)

The default **Analyze & Generate Reports** workflow now uses analysis/universal_pf2d.py.
Older fixed LMO (011) processing is still available from the command line with --analyze-legacy.

**Automatic scanning:** recursively reads RASX images from the selected folder AND its
subfolders. Headers provide two-theta, omega, Phi and Chi positions, motor sweep
start/stop, exposure (milliseconds), XY coordinates, detector dimensions and
physical calibration. It does not assume two chi values, 2-degree Phi steps or
four equally spaced poles. Compatible frames group by source folder, sample
name, XY, acquisition suffix, reflection-specific 2theta and detector settings.
Unsupported 1D scans and detector mountings are skipped with reasons in
Rejected_Files.csv; files are never deleted, renamed or changed.

**Scientific mapping:** known Rigaku SmartLabXE / HyPix-3000(H) horizontal
orientation only. It uses the actual frame-specific 2theta for the signal
ROI, adjacent background bands, and the measured detector-pixel geometry;
integrates recorded continuous Phi windows; normalizes counts to the actual
exposure time. Maps contain **only sampled directions**; unknown sectors are
grey and are never interpolated or completed by crystallographic symmetry.
Peaks are found from measured intensity rather than forced at 90-degree
intervals. Observed FWHM values are not instrumental-resolution-corrected.

**Limitations:** automatically reading arbitrary angles does not calibrate
a physically different detector orientation. Those acquisitions fail
explicitly. A single Chi band is not a full stereographic hemisphere,
and instrument beta zero is NOT calibrated to tape longitudinal direction.
Reflections (e.g. 222) are *filename hints*, not independently confirmed
phase identifications. Tilt coordinates are not an independent measurement
of tape-axis tilt. Compare with SmartLab-calibrated measurements before
using unfamiliar setups for a quantitative acceptance criterion.

**Results:** creates an independent timestamped folder outside the input
tree, with sample_###_*/Pole_Figure_Report.pdf, map.npz,
acquisition.json, chi_series.csv, peaks.json/peaks.csv, plus
Overall_Pole_Figure_Report.pdf, Summary.csv, All_Peaks.csv,
Rejected_Files.csv and inspection.json. Original inputs untouched.
The re-report button reuses saved maps for v1.4 reports and auto-routes
old v1.3 caches to the old formatter.

CLI examples:

~~~powershell
LMO-PF2D-Analyzer-v1.4.0.exe --inspect --source "D:\Raw\Measurements"
LMO-PF2D-Analyzer-v1.4.0.exe --analyze --source "D:\Raw\Measurements" --output "D:\PF2D\Run_NEW"
LMO-PF2D-Analyzer-v1.4.0.exe --analyze --reports-only --source "D:\PF2D\Run_OLD" --output "D:\PF2D\New_PDF"
LMO-PF2D-Analyzer-v1.4.0.exe --analyze-legacy --source "D:\LMO011" --output "D:\Legacy_NEW"
~~~

## Download / Installation

Download **LMO-PF2D-Analyzer-v1.4.1.exe** from the repository's [Releases](../../releases). No separate Python installation is required. Windows 10/11 x64 supported by the build (built on Windows 11).

1. Open the EXE.
2. Select the folder containing Rigaku RASX frames; the program recursively finds matching files.
3. Click **Проверить данные** (Check data) to list series and missing frames.
4. Select **Одна выбранная точка** (one point) or **Все точки** (all points).
5. Select a parent results directory and click **Построить отчёт**.
6. Open results. Every run creates a *new timestamped folder*; originals are never changed.

A successful single-point run creates `LMO011_point_report.pdf`, a measured stereographic pole figure, independent phi fits, peak CSV and geometry logs. A batch additionally creates `LMO011_XY_Overall_Report.pdf`, `XY_summary.csv` and `XY_summary.json`.

## Supported raw acquisition

The historical v1.0 pipeline deliberately validated a **specific, real, verified acquisition geometry**, rather than silently misinterpreting arbitrary 2D XRD images:

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
