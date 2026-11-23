"""Universal acquisition-aware Rigaku SmartLabXE / HyPix 2D pole figure reader.

All calculations use measured RASX metadata. Never modify acquisition files.
No extrapolation into unmeasured sectors and no instrument-resolution correction.
Only the explicitly supported detector mounting convention is reconstructed.
"""
from __future__ import annotations
import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import math
import pathlib
import re
import sys
import time
import traceback
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import dataclass

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import least_squares

VERSION = "1.4.0"
FRAME_RE = re.compile(
    r"^(?P<sample>.+)_image_(?P<seq>\d+)_(?P<chi>-?\d+(?:\.\d+)?)_(?P<frame>\d+)(?:_(?P<point>\d+))?$",
    re.IGNORECASE)
EXPECTED_VECTOR = "0 -1 1 0"
NOTE = ("Measured instrument-coordinate pole figure; uncovered sectors are not interpolated. "
        "Widths include instrumental and sampling broadening; azimuth zero is not "
        "calibrated to tape length.")
SIG = lambda x: round(float(x), 4)

def emit(stage, **items):
    print("PROGRESS_JSON:" + json.dumps(dict(stage=stage, **items), ensure_ascii=False), flush=True)

def num(value, default=None):
    try:
        f = float(value)
        return f if math.isfinite(f) else default
    except (TypeError, ValueError, OverflowError):
        return default

def vals(text):
    return [float(v) for v in str(text or "").split()]

def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("._")[:65] or "sample"

def parse_frame(path):
    """Open only ZIP metadata, not image binary. Reject incompatible geometries."""
    p = pathlib.Path(path)
    with zipfile.ZipFile(p, "r") as z:
        if "Data0/Image0.bin" not in z.namelist():
            raise ValueError("not a two-dimensional image RASX")
        root = ET.fromstring(z.read("Data0/MesurementConditions0.xml"))
    ax = {x.get("Name"): num(x.get("Position")) for x in root.findall("Axes/Axis")}
    pairs = root.findall("ImageInformation/DTHeader/Pair")
    hdr = {p[0].text: p[1].text for p in pairs if len(p) >= 2}
    scan = root.find("ScanInformation")
    if scan is None:
        raise ValueError("missing ScanInformation")
    axis = scan.findtext("AxisName") or ""
    if axis.casefold() != "phi":
        raise ValueError("not a Phi-axis 2D pole-figure acquisition")
    phi, chi, tt, omega = [ax.get(n) for n in ("Phi", "Chi", "TwoTheta", "Omega")]
    if any(v is None for v in (phi, chi, tt, omega)):
        raise ValueError("missing Phi/Chi/TwoTheta/Omega motor positions")
    shape = [int(x) for x in vals(hdr.get("PXD_DETECTOR_DIMENSIONS"))]
    if len(shape) != 2 or min(shape) <= 0 or max(shape) > 8192:
        raise ValueError(f"invalid pixel dimensions: {shape}")
    origin = vals(hdr.get("PXD_SPATIAL_DISTORTION_INFO"))
    gonio = vals(hdr.get("PXD_GONIO_VALUES"))
    if len(origin) < 4 or len(gonio) < 6 or gonio[5] <= 0:
        raise ValueError("missing physical detector origin, pixel pitch or distance")
    detector = hdr.get("PXD_DETECTOR_IDENTIFICATION") or hdr.get("PXD_DETECTOR_DESCRIPTION") or "unknown"
    # This matrix was established for horizontal HyPix geometry (SmartLabXE).
    # Never silently apply it to a rotated detector or a different mounting.
    if "hypix" not in detector.casefold() or "h)" not in detector.casefold():
        raise ValueError(f"unsupported detector orientation/model: {detector}")
    if (hdr.get("PXD_SPATIAL_DISTORTION_VECTORS", "").strip() != EXPECTED_VECTOR
            or abs(gonio[0] - 270) > 0.1
            or abs(gonio[2]) > 0.1
            or hdr.get("PXD_DETECTOR_IMAGE_DIRECTION", "0") != "0"):
        raise ValueError("unsupported detector mounting/orientation; calibration required")
    if origin[2] <= 0 or origin[3] <= 0:
        raise ValueError("non-positive pixel size")
    data_type = (hdr.get("Data_type") or "").casefold()
    if data_type not in ("unsigned long int", "unsigned short int"):
        raise ValueError(f"unsupported detector pixel datatype: {data_type}")
    byte_order = (hdr.get("BYTE_ORDER") or "").casefold()
    if byte_order not in ("little_endian", "big_endian"):
        raise ValueError(f"unknown detector byte order: {byte_order}")
    exposure = num(hdr.get("PXD_DETECTOR_ACQUISITION_EXPTIME"))
    if exposure is None or exposure <= 0:
        raise ValueError("no valid detector exposure time")
    start, stop = num(scan.findtext("Start")), num(scan.findtext("Stop"))
    mode = (scan.findtext("Mode") or "UNKNOWN").upper()
    sweep = abs(stop - start) if start is not None and stop is not None and mode == "CONTINUOUS" else 0.
    m = FRAME_RE.match(p.stem)
    sample = m.group("sample") if m else p.stem
    point = m.group("point") if m else None
    hint = re.search(r"(?<!\d)(\d{3})(?:_2D|\b)", sample, flags=re.IGNORECASE)
    return dict(path=str(p), source_name=p.name, sample=sample,
        reflection_hint=(hint.group(1) if hint else "not specified"),
        acquisition_suffix=point, seq=(int(m.group("seq")) if m else None),
        frame_no=(int(m.group("frame")) if m else None), 
        phi=phi % 360., chi=chi, x=ax.get("X"), y=ax.get("Y"),
        two_theta=tt, omega=omega, start=start, stop=stop,
        mode=mode, sweep_deg=sweep, exposure_ms=exposure,
        detector=detector, dimensions=shape, detector_origin=origin[:4],
        detector_distance=gonio[5], binary_type=data_type, byte_order=byte_order,
        scan_start_time=scan.findtext("StartTime"), 
        detector_description=hdr.get("PXD_DETECTOR_DESCRIPTION"),
        geometry_basis="HyPix horizontal / SmartLabXE / recorded physical calibration")

def files_under(source):
    root = pathlib.Path(source)
    if not root.is_dir():
        raise FileNotFoundError(f"Source folder not found: {root}")
    # Always recurse, including when parent contains its own raw files.
    return sorted((p for p in root.rglob("*") if p.is_file() and p.suffix.lower() == ".rasx"),
                  key=lambda p: str(p).lower())

def group_frames(rows):
    groups = collections.defaultdict(list)
    for r in rows:
        p = pathlib.Path(r["path"])
        xy = tuple(round(r[k], 3) if r[k] is not None else None for k in ("x", "y"))
        # Different reflections/omega and detectors must never share one map.
        signature = (round(r["two_theta"], 2), round(r["omega"], 2),
                     r["detector"], tuple(r["dimensions"]))
        key = (str(p.parent), r["sample"], xy, r["acquisition_suffix"], signature)
        groups[key].append(r)
    return dict(groups)

def discover(source):
    found = files_under(source)
    emit("file_scan", scan_phase="enumerate", found=len(found), total=len(found))
    if not found:
        raise ValueError("No .rasx files found recursively")
    rows, rejected = [], []
    last = time.monotonic()
    for i, p in enumerate(found, 1):
        try:
            rows.append(parse_frame(p))
        except Exception as e:
            rejected.append(dict(file=str(p), error=str(e)))
        if (i == len(found)) or (time.monotonic() - last > 4):
            emit("file_scan", scan_phase="headers", found=len(found), total=len(found),
                 checked=i, accepted=len(rows), rejected=len(rejected), current=p.name)
            last = time.monotonic()
    groups = group_frames(rows)
    return groups, rejected

def series_summary(rows):
    angles = np.array(sorted({round(float(r["phi"]) % 360., 6) for r in rows}), dtype=float)
    if not len(angles):
        raise ValueError("empty chi/phi group")
    gaps = np.diff(np.r_[angles, angles[0] + 360.]) if len(angles) > 1 else np.array([360.])
    candidates = gaps[gaps > 1e-5]
    step = float(np.median(np.sort(candidates)[:max(1, int(len(candidates) * .8))])) if len(candidates) else 360.
    # A circular wrap gap and gaps >1.5 nominal steps identify missed sectors.
    big = [(float(angles[i]), round(float(g), 4)) for i, g in enumerate(gaps)
           if g > step * 1.5 + 1e-5]
    unique = len(angles)
    modes = sorted(set(r["mode"] for r in rows))
    continuous = modes == ["CONTINUOUS"]
    exp = [r["exposure_ms"] for r in rows]
    actual_sweeps = [r["sweep_deg"] for r in rows if r["sweep_deg"]]
    dup = len(rows) - unique
    warnings = []
    if big: warnings.append(f"{len(big)} Phi gaps exceed nominal step (up to {max(g for _, g in big):.3f} deg)")
    if dup: warnings.append(f"{dup} duplicate Phi centers: repeated exposures retained in intensity averaging")
    if len(modes) != 1: warnings.append("mixed scan modes")
    if continuous and actual_sweeps and abs(float(np.median(actual_sweeps)) - step) > max(.05, step*.15):
        warnings.append("Phi sweep and motor step differ; inspect angular overlaps/gaps")
    full_turn = unique > 2 and not big
    return dict(n_frames=len(rows), unique_phi=unique, phi_step_deg=step, 
        phi_min_deg=float(angles[0]), phi_max_deg=float(angles[-1]),
        mode="/".join(modes), full_turn=bool(full_turn), phi_gaps=big[:100],
        duplicate_phi=dup, exposure_ms_min=float(min(exp)), exposure_ms_max=float(max(exp)),
        sweep_deg=float(np.median(actual_sweeps)) if actual_sweeps else 0.,
        warnings=warnings)

def validate_group(rows):
    series = collections.defaultdict(list)
    for r in rows:
        series[round(r["chi"], 4)].append(r)
    if not series:
        raise ValueError("no valid Chi series")
    first = rows[0]
    for r in rows:
        for k, tolerance in (("two_theta", .03), ("omega", .03), ("detector_distance", 1.0)):
            if abs(r[k] - first[k]) > tolerance:
                raise ValueError(f"motor/detector {k} varies within one group")
        if r["dimensions"] != first["dimensions"] or r["detector"] != first["detector"]:
            raise ValueError("detector changed within one group")
    summaries = {str(ch): series_summary(seq) for ch, seq in sorted(series.items())}
    warnings = [f"Chi={ch}: {w}" for ch, summary in summaries.items() for w in summary["warnings"]]
    if not all(s["full_turn"] for s in summaries.values()):
        warnings.append("Incomplete Phi coverage: generate partial measured figure only")
    if len(series) == 1:
        warnings.append("Only one Chi setting: limited radial pole coverage")
    if not (abs(first["two_theta"]-32.6)<.3 and abs(first["omega"]-16.3)<.3):
        warnings.append("Detector coordinate mapping extended beyond validated LMO(011); verify against SmartLab")
    return dict(series=series, per_chi=summaries, warnings=warnings, full_turn=all(s["full_turn"] for s in summaries.values()))

def detector_map(row):
    """Map detector pixels to specimen tilt/azimuth using the known HyPix H mounting."""
    width, height = row["dimensions"]
    u, v = np.meshgrid(np.arange(width), np.arange(height))
    cx, cy, pitch_x, pitch_y = row["detector_origin"]
    distance = float(row["detector_distance"])
    tt = math.radians(row["two_theta"])
    om = math.radians(row["omega"])
    du = (u - cx) * pitch_x
    dv = (v - cy) * pitch_y
    rr = np.stack((-du, -dv * math.cos(tt) - distance * math.sin(tt),
                   dv * math.sin(tt) - distance * math.cos(tt)), axis=-1)
    s = rr / np.linalg.norm(rr, axis=-1, keepdims=True)
    tth = np.degrees(np.arccos(np.clip(-s[..., 2], -1., 1.)))
    q = s + np.array([0., 0., 1.])
    q /= np.linalg.norm(q, axis=-1, keepdims=True)
    q0 = np.stack((q[..., 0], math.cos(om)*q[..., 1]-math.sin(om)*q[..., 2],
                   math.sin(om)*q[..., 1]+math.cos(om)*q[..., 2]), axis=-1)
    alpha = math.radians(90.-row["chi"])
    qs = np.stack((q0[..., 0]*math.cos(alpha)+q0[..., 1]*math.sin(alpha),
                  -q0[..., 0]*math.sin(alpha)+q0[..., 1]*math.cos(alpha), q0[..., 2]), axis=-1)
    qs = np.where((qs[..., 0] > 0)[..., None], -qs, qs)
    tilt = np.degrees(np.arccos(np.clip(-qs[..., 0], -1., 1.)))
    azimuth = np.degrees(np.arctan2(qs[..., 2], -qs[..., 1])) % 360.
    # Absolute 2theta window is anchored to this scan's motor position, never 32.6 deg.
    diff = tth - row["two_theta"]
    signal = (np.abs(diff) < 1.20) & (tilt >= 0.) & (tilt < 90.)
    bg_left = (diff > -3.0) & (diff < -1.8) & (tilt >= 0.) & (tilt < 90.)
    bg_right = (diff > 1.8) & (diff < 3.0) & (tilt >= 0.) & (tilt < 90.)
    if signal.sum() < 120:
        raise ValueError(f"Detector has insufficient two-theta ROI around {row['two_theta']} degrees")
    if bg_left.sum() < 50 or bg_right.sum() < 50:
        raise ValueError("Detector lacks sidebands for safe background subtraction")
    return tilt, azimuth, signal, bg_left, bg_right

def frame_image(row):
    width, height = row["dimensions"]
    dtype = ("<" if row["byte_order"] == "little_endian" else ">") + (
        "u4" if row["binary_type"] == "unsigned long int" else "u2")
    expected = width * height * np.dtype(dtype).itemsize
    with zipfile.ZipFile(row["path"], "r") as z:
        member = z.getinfo("Data0/Image0.bin")
        if member.file_size != expected:
            raise ValueError(f"Image length {member.file_size} inconsistent with {expected} at {row['source_name']}")
        raw = z.read(member)
    return np.frombuffer(raw, dtype=dtype).reshape(height, width)

def reconstruct(rows, dest, callback=None, phi_bins=180, tilt_bins=90):
    valid = validate_group(rows)
    grid = np.zeros((tilt_bins, phi_bins), dtype=float)
    exposure_weights = np.zeros_like(grid)
    per_chi = []
    total = len(rows)
    done = 0
    last_update = 0.
    bin_angle = 360. / phi_bins
    for chi, seq in sorted(valid["series"].items()):
        ref = seq[0]
        tilt, az, signal, bg1, bg2 = detector_map(ref)
        inds = np.flatnonzero(signal)
        ti = np.floor(tilt.flat[inds]).astype(np.intp)
        ai = az.flat[inds]
        base = ti * phi_bins
        accepted = 0
        sum_raw = 0.
        for row in sorted(seq, key=lambda r:(r["phi"],r["source_name"])):
            image = frame_image(row)
            bkg = .5 * (float(image[bg1].mean()) + float(image[bg2].mean()))
            rate = (image.flat[inds].astype(float) - bkg) / (row["exposure_ms"] / 1000.)
            sum_raw += float(image[signal].sum())
            if row["mode"] == "CONTINUOUS" and row["sweep_deg"] > 0:
                # 5-point quadrature for the actual Phi window, including the wrap.
                phis = np.linspace(row["start"], row["stop"], 5)
            else:
                phis = [row["phi"]]
            for phi in phis:
                idx = base + (np.floor((ai + phi + bin_angle/2.) % 360. / bin_angle)
                              .astype(np.intp) % phi_bins)
                grid += np.bincount(idx, weights=rate/len(phis),
                                    minlength=phi_bins*tilt_bins).reshape(tilt_bins,phi_bins)
                exposure_weights += np.bincount(
                    idx, weights=np.full(len(idx), 1./len(phis)),
                    minlength=phi_bins*tilt_bins).reshape(tilt_bins,phi_bins)
            done += 1
            accepted += 1
            now = time.monotonic()
            if callback and (done == total or now-last_update >= 3.):
                callback(done,total,row["source_name"])
                last_update = now
        per_chi.append(dict(chi_deg=chi, frames=accepted, sum_raw_counts=sum_raw,
                            roi_pixel_count=len(inds), bg_pixel_count=int(bg1.sum()+bg2.sum())))
    intensity = np.divide(grid, exposure_weights, out=np.full_like(grid,np.nan), where=exposure_weights > 0)
    np.savez_compressed(dest/"map.npz", intensity=intensity, coverage=exposure_weights,
        tilt=np.arange(tilt_bins)+.5, azimuth=np.arange(phi_bins)*bin_angle+bin_angle/2.)
    return intensity, exposure_weights, dict(valid, chi_stats=per_chi)

def fit_peaks(data, coverage, n_max=10):
    """Locate measured local maxima; do not assume four peaks or perfect 90-degree spacing."""
    measurable = (coverage > 0) & np.isfinite(data)
    if np.count_nonzero(measurable) < 15:
        return []
    weight = gaussian_filter(measurable.astype(float),(1.4,1.2), mode=("constant","wrap"))
    smooth = gaussian_filter(np.where(measurable,data,0.),(1.4,1.2),
                             mode=("constant","wrap")) / np.maximum(weight,.001)
    reliable = measurable & (weight > .25)
    # Invalidate non-measured bins BEFORE local maxima search: otherwise
    # Gaussian division by tiny edge weights creates phantom neighboring maxima.
    smooth[~reliable] = -np.inf
    z = smooth[reliable]
    if len(z) < 15:
        return []
    med = float(np.median(z))
    p99 = float(np.percentile(z,99))
    mad = float(np.median(np.abs(z-med)))
    threshold = med + max(2.0*1.4826*mad, .22 * max(p99-med,0.), .015)
    local = maximum_filter(smooth,size=(7,7),mode=("constant","wrap"))
    candidates = np.argwhere((smooth >= local - 1e-10) & reliable & (smooth >= threshold))
    if not len(candidates):
        return []
    candidates = sorted(candidates,key=lambda ix:smooth[ix[0],ix[1]],reverse=True)
    peaks=[]
    n_t,n_a=data.shape
    az_centers=(np.arange(n_a)+.5)*360./n_a
    tilts=np.arange(n_t)+.5
    tt,aa=np.meshgrid(tilts,az_centers,indexing="ij")
    for i,j in candidates:
        t0=float(tilts[i]); a0=float(az_centers[j])
        if any(abs(t0-p["tilt_deg"])<7. and abs(((a0-p["beta_deg"]+180)%360)-180)<13.
               for p in peaks):
            continue
        delta=((aa-a0+180)%360)-180
        mask = measurable & (np.abs(delta)<=24.) & (np.abs(tt-t0)<=14.)
        if mask.sum() < 18:
            continue
        xx=delta[mask]
        yy=tt[mask]
        zz=data[mask]
        local_base=float(np.percentile(zz,25))
        amplitude=max(float(smooth[i,j]-local_base),.001)
        def res(params):
            amp,b0,ct,sb,st,bg=params
            pred=amp*np.exp(-.5*((xx-b0)/sb)**2-.5*((yy-ct)/st)**2)+bg
            return pred-zz
        try:
            fit=least_squares(res,[amplitude,0.,t0,5.,3.,local_base],
                bounds=([0.,-12.,max(0.,t0-12.),.4,.4,-np.inf],
                        [np.inf,12.,min(89.,t0+12.),28.,17.,np.inf]),max_nfev=250)
            amp,b0,ct,sb,st,bg=map(float,fit.x)
            # FWHM below two map bins is not physically resolved; this is
            # common for isolated edge bins and sparse detector statistics.
            # Do not output such spots as measured crystallographic peaks.
            radial_fwhm = st * (2*math.sqrt(2*math.log(2)))
            phi_fwhm = sb * (2*math.sqrt(2*math.log(2)))
            if radial_fwhm < 2.0 or phi_fwhm < 4.0:
                continue
            err=res(fit.x)
            sse=float(np.dot(err,err))
            sst=float(np.sum((zz-zz.mean())**2))
            r2=(1.-sse/sst) if sst>1e-12 else None
            if r2 is None or r2 < 0.65:
                continue
            fwhm=2*math.sqrt(2*math.log(2))
            peaks.append(dict(rank=0, status="OK" if fit.success else "FIT_WARNING",
                beta_deg=float((a0+b0)%360),tilt_deg=ct,
                stereo_radius=float(math.tan(math.radians(ct/2))),
                phi_fwhm_observed_deg=sb*fwhm, tilt_fwhm_observed_deg=st*fwhm,
                sigma_beta_deg=sb, sigma_tilt_deg=st, amplitude=amp,
                background=bg, fit_r2=r2, fit_rmse=float(math.sqrt(sse/len(zz))),
                n_fit_bins=int(len(zz)), seed_beta_deg=a0, seed_tilt_deg=t0))
        except Exception as exc:
            peaks.append(dict(rank=0,status="FIT_FAILED", beta_deg=a0,tilt_deg=t0, error=str(exc)))
        if len(peaks)>=n_max:
            break
    for j,p in enumerate(peaks,1):
        p["rank"]=j
    return peaks

def export_csv(path,rows,keys=None):
    if keys is None:
        keys=sorted(set().union(*(r.keys() for r in rows))) if rows else ["status"]
    with pathlib.Path(path).open("w",newline="",encoding="utf-8-sig") as f:
        writer=csv.DictWriter(f,fieldnames=keys,extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k:(json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v)
                            for k,v in row.items()})

def render_map(ax,data,coverage,peaks=None,title=None, vmax=None):
    """Stereographic pole projection with no visual interpolation over missing sectors."""
    n_t,n_a=data.shape
    theta=np.linspace(0.,2*math.pi,n_a+1)
    radius=np.tan(np.radians(np.arange(n_t+1)/2.))
    image=np.ma.masked_where((coverage<=0)|~np.isfinite(data),data)
    ax.set_facecolor("#e8ebed")
    mesh=ax.pcolormesh(theta,radius,image,cmap="magma",shading="flat",vmin=0.,vmax=vmax,
                       rasterized=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0.,1.)
    ax.set_yticks(np.tan(np.radians(np.array([15,30,45,60,75,90])/2)))
    ax.set_yticklabels(["15Â°","30Â°","45Â°","60Â°","75Â°","90Â°"],fontsize=7)
    ax.tick_params(labelsize=8)
    if title:ax.set_title(title,pad=14,fontsize=11)
    for peak in peaks or []:
        if peak.get("status","").startswith("FIT_") or peak.get("status")=="OK":
            if "beta_deg" not in peak:continue
            b=math.radians(peak["beta_deg"])
            rad=math.tan(math.radians(peak["tilt_deg"]/2))
            ax.plot([b],[rad],"o",ms=7,mfc="none",mec="#38ffff",mew=1.3)
            ax.annotate(str(peak["rank"]),(b,rad),xytext=(3,4),
                        textcoords="offset points",fontsize=8,color="#38ffff")
    return mesh

def _page_header(fig, title, subtitle=""):
    fig.text(.055,.965,title,fontsize=17,fontweight="bold",color="#253952",va="top")
    if subtitle: fig.text(.055,.928,subtitle,fontsize=9,color="#58677b",va="top")
    fig.text(.055,.025,NOTE,fontsize=7,color="#6d7580",va="bottom",wrap=True)

def _table(ax,headers,rows,col_widths=None,font=8):
    ax.axis("off")
    if not rows: rows=[["No entries"]+[""]*(len(headers)-1)]
    table=ax.table(cellText=[[str(x) for x in row] for row in rows],
                   colLabels=headers,cellLoc="left",colLoc="left",loc="upper center",
                   colWidths=col_widths)
    table.auto_set_font_size(False)
    table.set_fontsize(font)
    table.scale(1,1.55)
    for (r,c),cell in table.get_celld().items():
        cell.set_edgecolor("#d5dce3")
        cell.set_linewidth(.35)
        if r==0:
            cell.set_facecolor("#284761")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        elif r%2==0:
            cell.set_facecolor("#eef3f7")

def write_point_pdf(dest, item, data, cov, peaks, acquisition):
    sample=item["sample"]
    report_file=dest/"Pole_Figure_Report.pdf"
    nvalid=int(np.count_nonzero(cov>0))
    ntot=int(cov.size)
    pct=100*nvalid/max(ntot,1)
    cap = f"{sample} | 2theta={item['two_theta_deg']:.3f}Â° | {item['frames']} 2D frames"
    finite=data[np.isfinite(data) & (cov>0)]
    vmax=float(np.percentile(finite,99.7)) if finite.size else 1.
    vmax=max(vmax, .0001)
    with PdfPages(report_file) as pdf:
        fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
        _page_header(fig,"Measured 2D pole figure",cap)
        ax=fig.add_axes([.08,.13,.52,.73],projection="polar")
        mesh=render_map(ax,data,cov,peaks,vmax=vmax)
        cax=fig.add_axes([.60,.19,.016,.55])
        fig.colorbar(mesh,cax=cax,label="net counts / pixel / second")
        info=[
            f"Sample/group: {sample}",
            f"Reflection label: {item['reflection_hint']} (filename only)",
            f"Detector: {item['detector']}",
            f"2theta: {item['two_theta_deg']:.3f}Â°; omega: {item['omega_deg']:.3f}Â°",
            f"Chi series: {', '.join(f'{v:g}Â°' for v in item['chi_angles_deg'])}",
            f"Frames: {item['frames']}; peaks fitted: {len(peaks)}",
            f"Measured map bins: {nvalid}/{ntot} ({pct:.1f}%)",
            f"Phi coverage: {item['phi_completeness']}",
            f"Status: {item['status']}",
            f"XY: X={item['x_mm']} mm, Y={item['y_mm']} mm",
            "",
            "Unmeasured bins are grey and not inferred.",
            "Rotation reference is instrument beta.",
        ]
        fig.text(.65,.82,"\n".join(info),va="top",fontsize=9,linespacing=1.7)
        pdf.savefig(fig);plt.close(fig)
        fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
        _page_header(fig,"Peak locations and widths",cap)
        ax=fig.add_axes([.06,.36,.90,.47])
        denom=(cov>0).sum(axis=0); profile=np.divide(np.where(cov>0,data,0.).sum(axis=0),denom,out=np.full(data.shape[1],np.nan),where=denom>0)
        x=(np.arange(len(profile))+.5)*360./len(profile)
        ax.plot(x,profile,lw=1.35,color="#254f7c",label="Measured intensity (covered tilt bins)")
        for p in peaks:
            if "beta_deg" in p:
                ax.axvline(p["beta_deg"],color="#d77b2c",ls="--",lw=.8,alpha=.8)
                ax.text(p["beta_deg"],.96,p["rank"],transform=ax.get_xaxis_transform(),
                        ha="center",va="top",fontsize=8,color="#a55416")
        ax.set_xlim(0,360);ax.set_xticks(np.arange(0,361,30));ax.grid(alpha=.2)
        ax.set_xlabel("Instrument azimuth beta / degrees")
        ax.set_ylabel("Net counts / pixel / second")
        ax.set_title("Azimuthal profile (descriptive; not an independent FWHM fit)",fontsize=10)
        peak_rows=[]
        for p in peaks:
            def fmt(name):
                v=p.get(name)
                return f"{v:.2f}" if isinstance(v,(int,float)) else "-"
            peak_rows.append([p["rank"],fmt("beta_deg"),fmt("tilt_deg"),
                              fmt("phi_fwhm_observed_deg"),fmt("tilt_fwhm_observed_deg"),
                              fmt("fit_r2"),p.get("status","?")])
        tab=fig.add_axes([.06,.08,.90,.23])
        _table(tab,["#","Beta (deg)","Tilt (deg)","FWHM beta","FWHM tilt","R2","Fit QC"],
               peak_rows,font=8)
        pdf.savefig(fig);plt.close(fig)
        fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
        _page_header(fig,"Acquisition settings and coverage",cap)
        ax=fig.add_axes([.07,.41,.86,.39])
        coverage_profile=np.sum(cov>0,axis=0)
        xx=(np.arange(len(coverage_profile))+.5)*360./len(coverage_profile)
        ax.bar(xx,coverage_profile,width=360./len(xx),color="#4d83a5")
        ax.set_xlim(0,360);ax.set_xticks(np.arange(0,361,30))
        ax.set_xlabel("Instrument azimuth beta / degrees")
        ax.set_ylabel("Number of sampled tilt bins")
        ax.set_title("Actual detector coverage - no values substituted for gaps",fontsize=10)
        ax.grid(alpha=.18,axis="y")
        rows=[]
        for chi,meta in sorted(acquisition["per_chi"].items(),key=lambda x:float(x[0])):
            rows.append([f"{float(chi):g}", meta["n_frames"], f"{meta['phi_step_deg']:g}",
                         f"{meta['sweep_deg']:g}",meta["mode"],
                         "Yes" if meta["full_turn"] else "NO",
                         f"{meta['exposure_ms_min']:.0f}-{meta['exposure_ms_max']:.0f}"])
        tab=fig.add_axes([.06,.20,.90,.17])
        _table(tab,["Chi (deg)","Frames","Phi step","Sweep","Mode","360Â°","Exposure (ms)"],rows,font=8)
        warn=acquisition["warnings"]
        fig.text(.07,.16,"QUALITY / GEOMETRY NOTES",fontsize=9,fontweight="bold",color="#425670")
        message="\n".join("â€¢ "+x for x in warn[:6]) if warn else "â€¢ No acquisition gaps or metadata inconsistencies detected."
        fig.text(.07,.14,message,va="top",fontsize=8,color="#72552a",linespacing=1.25)
        pdf.savefig(fig);plt.close(fig)
        # Local maps and radial/azimuthal cuts around each actual fitted maximum.
        for start in range(0,len(peaks),4):
            fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
            _page_header(fig,"Measured peak neighborhoods",cap)
            for n,p in enumerate(peaks[start:start+4]):
                ax=fig.add_subplot(2,2,n+1)
                if "beta_deg" not in p:
                    ax.text(.5,.5,str(p.get("error","Fit failed")),ha="center",transform=ax.transAxes)
                    continue
                beta=p["beta_deg"]; tilt=p["tilt_deg"]
                nt,na=data.shape
                at=(np.arange(na)+.5)*360./na
                tv=np.arange(nt)+.5
                d=((at-beta+180)%360)-180
                cols=np.where(np.abs(d)<17)[0]
                cols=cols[np.argsort(d[cols])]
                rws=np.where(np.abs(tv-tilt)<12)[0]
                sub=data[np.ix_(rws,cols)]
                mask=cov[np.ix_(rws,cols)]<=0
                ax.pcolormesh(d[cols],tv[rws],np.ma.masked_where(mask,sub),
                              cmap="magma",shading="nearest",vmin=0,vmax=vmax)
                ax.plot(0,tilt,"c+",markersize=10,mew=1.6)
                ax.set_xlabel("Beta - peak center (deg)",fontsize=8)
                ax.set_ylabel("Tilt (deg)",fontsize=8)
                ax.set_title(f"Peak {p['rank']}: {beta:.1f}Â° / {tilt:.1f}Â° | "
                             f"FWHM {p.get('phi_fwhm_observed_deg',float('nan')):.1f}Â°",fontsize=9)
            fig.tight_layout(rect=[.04,.06,.98,.87])
            pdf.savefig(fig);plt.close(fig)
    return str(report_file)

def summarize_row(item):
    subset={k:item.get(k) for k in ("sample","reflection_hint","position_suffix",
            "x_mm","y_mm","two_theta_deg","omega_deg","detector",
            "frames","chi_angles_deg","phi_completeness","n_peaks","status","reason","warnings")}
    subset["chi_angles_deg"]=";".join(str(x) for x in item.get("chi_angles_deg",[]))
    subset["warnings"]="; ".join(item.get("warnings",[]))
    return subset

def write_overall_report(output,records):
    """Comparative QC across all discovered groups, including partial/failed ones."""
    report=output/"Overall_Pole_Figure_Report.pdf"
    ok=[r for r in records if r.get("status") in ("OK","PARTIAL")]
    with PdfPages(report) as pdf:
        fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
        _page_header(fig,"Batch pole-figure report",f"Rigaku HyPix | Universal PF2D Analyzer v{VERSION}")
        counts=collections.Counter(r.get("status","UNKNOWN") for r in records)
        txt=(f"Total sample/acquisition groups: {len(records)}\n"
             f"Complete: {counts.get('OK',0)}    Partial: {counts.get('PARTIAL',0)}"
             f"    Failed: {counts.get('FAILED',0)}\n"
             "Different 2theta reflections and different detector geometries are NEVER combined.")
        fig.text(.07,.84,txt,fontsize=11,linespacing=1.9,color="#26425e",va="top")
        tab=fig.add_axes([.055,.18,.91,.55])
        table_rows=[[str(r.get("point","")),str(r.get("sample",""))[:27],
            r.get("reflection_hint",""),str(r.get("frames","")),
            "; ".join(f"{ch:g}" for ch in r.get("chi_angles_deg",[])),
            str(r.get("n_peaks",0)),r.get("status","")] for r in records[:20]]
        _table(tab,["#","Sample/group","(hkl)*","2D frames","Chi settings (Â°)","Peaks","QC"],
            table_rows,col_widths=[.04,.32,.08,.09,.23,.08,.12],font=8)
        fig.text(.06,.08,"* Reflection index is a filename hint, not a measured phase identification.",fontsize=9)
        if len(records)>20:
            fig.text(.06,.06,f"First 20 of {len(records)} shown; complete machine-readable table in Summary.csv.")
        pdf.savefig(fig);plt.close(fig)
        # Large batches: one extra table page for every 28 records beyond cover.
        for start in range(20,len(records),28):
            fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
            _page_header(fig,"Group inventory",f"Groups {start+1} to {min(start+28,len(records))}")
            tab=fig.add_axes([.055,.07,.91,.78])
            rows=[[str(r.get("point","")),str(r.get("sample",""))[:28],
                f"{r.get('two_theta_deg',0):.2f}" if r.get("two_theta_deg") is not None else "-",
                str(r.get("frames","")),str(r.get("n_peaks",0)),r.get("status","")] for r in records[start:start+28]]
            _table(tab,["#","Sample/group","2theta (Â°)","Frames","Peaks","QC"],rows,
                   col_widths=[.06,.44,.13,.12,.10,.12],font=8)
            pdf.savefig(fig);plt.close(fig)
        # Representative measured maps with fixed scale per sample (not across different exposures).
        for start in range(0,len(ok),6):
            fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
            _page_header(fig,"Measured pole figures: comparison",f"Groups {start+1}â€“{min(start+6,len(ok))} of {len(ok)}")
            for k,row in enumerate(ok[start:start+6]):
                fp=output/row["folder"]/"map.npz"
                try:
                    with np.load(fp,allow_pickle=False) as z:
                        data=z["intensity"];cov=z["coverage"]
                    subplot=fig.add_subplot(2,3,k+1,projection="polar")
                    values=data[np.isfinite(data) & (cov>0)]
                    vmax=float(np.percentile(values,99.7)) if values.size else 1.
                    render_map(subplot,data,cov,title=(row["sample"][:24]+" ("+row["status"]+")"),vmax=max(vmax,.0001))
                except Exception as exc:
                    subplot=fig.add_subplot(2,3,k+1)
                    subplot.text(.1,.5,f"Cannot load map: {exc}",wrap=True)
            fig.tight_layout(rect=[.02,.05,.98,.87],h_pad=3.2,w_pad=2.)
            pdf.savefig(fig);plt.close(fig)
        if ok:
            fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
            _page_header(fig,"Sample comparison: observed peak metrics","First fitted peak in each sample; do not compare unlike reflections as texture equivalence")
            idx=np.arange(len(ok))
            top=[next((p for p in r.get("peaks",[]) if p.get("status")=="OK"),{}) for r in ok]
            beta=np.array([p.get("beta_deg",np.nan) for p in top])
            tilt=np.array([p.get("tilt_deg",np.nan) for p in top])
            fphi=np.array([p.get("phi_fwhm_observed_deg",np.nan) for p in top])
            frad=np.array([p.get("tilt_fwhm_observed_deg",np.nan) for p in top])
            axs=fig.subplots(2,2)
            for ax,vals,title in zip(axs.flat,(beta,tilt,fphi,frad),
                                     ("Beta center (deg)","Tilt coordinate (deg)",
                                      "Observed beta FWHM (deg)","Observed radial FWHM (deg)")):
                ax.plot(idx+1,vals,"o-",lw=1.2,ms=3)
                ax.set_title(title,fontsize=10)
                ax.set_xlabel("Sample/acquisition index",fontsize=8)
                ax.grid(alpha=.2)
            fig.tight_layout(rect=[.04,.07,.98,.87])
            pdf.savefig(fig);plt.close(fig)
        fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
        _page_header(fig,"Interpretation and quality limits","Do not interpret unmeasured sectors as zero intensity")
        notes=[
            "Acquisition settings are extracted separately from every RASX frame (scan motor metadata and detector calibration).",
            "2D counts are background-corrected using two-theta sidebands and normalized by individual frame exposure time.",
            "Every Phi sweep is integrated over its recorded angular interval, not inferred from a filename sequence.",
            "Different reflections, detector configurations, XY positions and export series are separated.",
            "No symmetry completion, pole-equivalent averaging, instrumental-resolution subtraction or texture ODF inversion is applied.",
            "Peak counts and observed FWHM depend on angular coverage; poor coverage can cause warnings or missing fits.",
            "True sample-axis tilt requires independent tape-axis and instrument-zero calibration, NOT provided by a pole map alone.",
            "The pixel-to-pole projection uses the supported HyPix horizontal SmartLabXE mounting convention.",
            "Other physical detector mountings are explicitly rejected until independently calibrated; numeric output is not guessed.",
        ]
        fig.text(.08,.85,"\n\n".join("â€¢ "+x for x in notes),fontsize=11,color="#3f5063",va="top",linespacing=1.6)
        pdf.savefig(fig);plt.close(fig)
    return str(report)

def save_report_package(out,records,rejected):
    payload=dict(version=VERSION,created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                 groups=records,rejected_files=rejected,notes=NOTE)
    (out/"inspection.json").write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding="utf-8")
    export_csv(out/"Summary.csv",[summarize_row(x) for x in records])
    export_csv(out/"Rejected_Files.csv",rejected,["file","error"])
    peak_table=[]
    for row in records:
        for p in row.get("peaks",[]):
            peak_table.append(dict(point=row.get("point"),sample=row.get("sample"),
                two_theta_deg=row.get("two_theta_deg"),sample_status=row.get("status"),**p))
    export_csv(out/"All_Peaks.csv",peak_table)
    return write_overall_report(out,records)

def _safe_out(source,destination):
    root=pathlib.Path(source).resolve()
    output=pathlib.Path(destination).resolve()
    if output.exists():
        raise FileExistsError(f"Output folder already exists: {output}")
    try:
        output.relative_to(root)
        raise ValueError("Output folder must be OUTSIDE the source tree (raw data read-only)")
    except ValueError as exc:
        if "OUTSIDE" in str(exc): raise
    output.parent.mkdir(parents=True,exist_ok=True)
    output.mkdir()
    return output

def _select_jobs(groups,group_key=None,position=None):
    ordered=sorted(groups.items(),key=lambda kv: (
        int(kv[0][3]) if kv[0][3] is not None else 10**9,
        str(kv[0][0]).lower(),str(kv[0][1]).lower(),str(kv[0][2]),str(kv[0][4])))
    jobs=[]
    for n,(key,rows) in enumerate(ordered,1):
        idx=int(key[3]) if key[3] is not None else n
        if group_key is not None and group_key not in (str(key),str(n),key[1]):
            continue
        if position is not None and position not in (n,idx):
            continue
        jobs.append((n,key,rows))
    if not jobs:
        raise ValueError("No sample/acquisition group matches selected position")
    if position is not None and len(jobs)!=1:
        raise ValueError("Position selection ambiguous: use --group-key for an exact group")
    return ordered,jobs

def run(source,output,inspect=False,group_key=None,position=None):
    emit("start",source=str(source),version=VERSION)
    groups,rejected=discover(source)
    if not groups:
        raise ValueError(f"No compatible HyPix 2D Phi RASX frames ({len(rejected)} rejected). "
                         "See rejection details or select a folder with 2D frames.")
    ordered,jobs=_select_jobs(groups,group_key,position)
    emit("metadata_done",groups=len(ordered),selected=len(jobs),
         valid_files=sum(len(v) for v in groups.values()),errors=len(rejected))
    for n,key,rows in jobs:
        emit("point_found",point=n,frames=len(rows),x=key[2][0],y=key[2][1],
             sample=key[1],reflection=rows[0]["reflection_hint"],
             chi_series=", ".join(str(x) for x in sorted({r["chi"] for r in rows})),
             two_theta=rows[0]["two_theta"])
    if inspect:
        for n,key,rows in jobs:
            qc=validate_group(rows)
            print(json.dumps(dict(group=n,sample=key[1],n_frames=len(rows),
                 chi=qc["per_chi"],warnings=qc["warnings"]),ensure_ascii=False),flush=True)
        return 0
    out=_safe_out(source,output)
    records=[]
    for j,(n,key,rows) in enumerate(jobs,1):
        first=rows[0]
        label=key[1]
        folder=f"sample_{n:03d}_{safe_name(label)}"
        item=dict(point=n,folder=folder,sample=label,
                  reflection_hint=first["reflection_hint"],position_suffix=key[3],
                  source_folder=key[0],source_group=str(key),x_mm=key[2][0],y_mm=key[2][1],
                  frames=len(rows),two_theta_deg=first["two_theta"],omega_deg=first["omega"],
                  detector=first["detector"],chi_angles_deg=sorted({r["chi"] for r in rows}),
                  n_peaks=0,peaks=[],status="FAILED",warnings=[])
        emit("point_start",point=n,index=j,total=len(jobs),frames=len(rows),sample=label)
        try:
            qc=validate_group(rows)
            item["warnings"]=list(qc["warnings"])
            item["phi_completeness"]="complete" if qc["full_turn"] else "partial"
            item["status"]="OK" if qc["full_turn"] else "PARTIAL"
            dest=out/folder
            dest.mkdir()
            emit("reconstruct",point=n,index=j,total=len(jobs),frames=len(rows))
            data,cov,extra=reconstruct(rows,dest,
                callback=lambda done,total,name:emit("frame_progress",
                    point=n,index=j,total=len(jobs),frames_done=done,
                    frames_total=total,current_file=name))
            emit("fit_and_report",point=n,index=j,total=len(jobs))
            peaks=fit_peaks(data,cov)
            if not peaks: item["warnings"].append("No reliable 2D peak found at the measured coverage")
            item.update(peaks=peaks,n_peaks=len(peaks))
            acquisition=dict(version=VERSION,sample=label,n_frames=len(rows),
                per_chi=qc["per_chi"],warnings=item["warnings"],
                chi_stats=extra["chi_stats"],motor_reference=first,
                method="Measured 2D HyPix projection + exposure-normalized counts / second",
                map_grid=dict(tilt_bin_deg=1.,azimuth_bin_deg=2.),notes=NOTE)
            (dest/"acquisition.json").write_text(
                json.dumps(acquisition,indent=2,ensure_ascii=False),encoding="utf-8")
            (dest/"peaks.json").write_text(json.dumps(peaks,indent=2),encoding="utf-8")
            export_csv(dest/"peaks.csv",peaks)
            export_csv(dest/"chi_series.csv",extra["chi_stats"])
            write_point_pdf(dest,item,data,cov,peaks,acquisition)
        except Exception as exc:
            item.update(status="FAILED",reason=str(exc))
            item["warnings"].append(f"Analysis error: {exc}")
            (out/f"sample_{n:03d}_FAILED.txt").write_text(traceback.format_exc(),encoding="utf-8")
        records.append(item)
        emit("point_done",point=n,index=j,total=len(jobs),status=item["status"],
             reason=item.get("reason",""),sample=label)
    emit("summary",total=len(records))
    overall=save_report_package(out,records,rejected)
    failed=sum(x["status"]=="FAILED" for x in records)
    partial=sum(x["status"]=="PARTIAL" for x in records)
    emit("complete",total=len(records),successful=len(records)-failed,failed=failed,
         partial=partial,output=str(out),report=overall)
    return 2 if failed else 0

def rebuild(source,output):
    src=pathlib.Path(source)
    meta=json.loads((src/"inspection.json").read_text(encoding="utf-8"))
    if meta.get("version")!=VERSION:
        raise ValueError("Saved maps have a different schema/version; use v1.3 reports-only for legacy maps")
    out=_safe_out(source,output)
    rows=[]
    for item in meta["groups"]:
        row=dict(item)
        if row["status"] in ("OK","PARTIAL"):
            old=src/row["folder"]
            new=out/row["folder"]
            new.mkdir()
            try:
                with np.load(old/"map.npz",allow_pickle=False) as z:
                    data=z["intensity"];cov=z["coverage"]
                acquisition=json.loads((old/"acquisition.json").read_text(encoding="utf-8"))
                np.savez_compressed(new/"map.npz",intensity=data,coverage=cov)
                if (old/"chi_series.csv").is_file():
                    import shutil
                    shutil.copy2(old/"chi_series.csv",new/"chi_series.csv")
                (new/"acquisition.json").write_text(json.dumps(acquisition,indent=2),encoding="utf-8")
                peaks=fit_peaks(data,cov)
                row["peaks"]=peaks
                row["n_peaks"]=len(peaks)
                (new/"peaks.json").write_text(json.dumps(peaks,indent=2),encoding="utf-8")
                export_csv(new/"peaks.csv",peaks)
                write_point_pdf(new,row,data,cov,peaks,acquisition)
            except Exception as exc:
                row.update(status="FAILED",reason="Reports rebuild: "+str(exc))
        rows.append(row)
    save_report_package(out,rows,meta.get("rejected_files",[]))
    return 2 if any(r["status"]=="FAILED" for r in rows) else 0

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True)
    p.add_argument("--output")
    p.add_argument("--inspect",action="store_true")
    p.add_argument("--reports-only",action="store_true")
    p.add_argument("--group-key")
    p.add_argument("--position",type=int)
    args=p.parse_args()
    if not args.inspect and not args.output:
        p.error("--output required for analysis or report rebuild")
    if args.reports_only:
        sys.exit(rebuild(args.source,args.output))
    sys.exit(run(args.source,args.output,args.inspect,args.group_key,args.position))

if __name__=="__main__":
    main()
