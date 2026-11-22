"""Acquisition-aware read-only Rigaku HyPix-3000(H) 2D pole figure processor.
Supports stationary phi step and continuous phi sweep with measured per-frame angles.
Never controls Rigaku hardware or writes into input folders.
"""
from __future__ import annotations
import argparse
import collections
import csv
import datetime as dt
import json
import math
import pathlib
import runpy
import shutil
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import least_squares

def normalize_phi(x):
    return float(round(float(x)%360.,8))%360.

def _f(x, default=None):
    try: return float(x)
    except (ValueError, TypeError): return default

def parse_rasx_header(path):
    """Read observed motor positions and true per-frame scan window from RASX ZIP."""
    path=pathlib.Path(path)
    with zipfile.ZipFile(path) as archive:
        root=ET.fromstring(archive.read("Data0/MesurementConditions0.xml"))
    axes={x.get("Name"): _f(x.get("Position")) for x in root.findall("Axes/Axis")}
    header={x[0].text: x[1].text for x in root.findall("ImageInformation/DTHeader/Pair")}
    scan=root.find("ScanInformation")
    mode=(scan.findtext("Mode") if scan is not None else "") or ""
    start=_f(scan.findtext("Start") if scan is not None else None)
    stop=_f(scan.findtext("Stop") if scan is not None else None)
    sweep=abs(stop-start) if start is not None and stop is not None and mode.upper()=="CONTINUOUS" else 0.
    phi=axes.get("Phi")
    chi=axes.get("Chi")
    if phi is None or chi is None:raise ValueError(f"Missing Phi/Chi motor metadata in {path.name}")
    dims=[int(float(v)) for v in header.get("PXD_DETECTOR_DIMENSIONS","").split()]
    if dims!=[775,385]:raise ValueError(f"Unverified detector dimensions {dims} in {path.name}; cannot reconstruct safely")
    if header.get("BYTE_ORDER")!="little_endian" or header.get("Data_type")!="unsigned long int":
        raise ValueError(f"Unsupported detector binary format in {path.name}")
    origin=[float(v) for v in header.get("PXD_SPATIAL_DISTORTION_INFO","").split()]
    gonio=[float(v) for v in header.get("PXD_GONIO_VALUES","").split()]
    if len(origin)<4 or len(gonio)<6:raise ValueError(f"Missing detector geometry in {path.name}")
    return dict(path=str(path),phi=normalize_phi(phi),chi=chi,x=axes.get("X"),y=axes.get("Y"),
      omega=axes.get("Omega"),two_theta=axes.get("TwoTheta"),mode=mode.upper(),
      sweep_deg=sweep,scan_start=start,scan_stop=stop,exposure_ms=_f(header.get("PXD_DETECTOR_ACQUISITION_EXPTIME")),
      detector_origin=origin[:4],detector_distance=gonio[5],detector_geometry=gonio,source_name=path.name)

def _median_step(phi):
    ordered=sorted({normalize_phi(p) for p in phi})
    if len(ordered)<2:return None,[]
    diff=np.diff(ordered+[ordered[0]+360.])
    common=float(np.median(diff))
    return common,diff.tolist()

def infer_series(rows):
    """Derive count, grid and missing angles from metadata. Never fill missing frames."""
    if not rows:raise ValueError("No frames")
    phi=[normalize_phi(x['phi']) for x in rows]
    step,spacings=_median_step(phi)
    if step is None:raise ValueError("Insufficient azimuth samples")
    unique=sorted(set(phi))
    nominal=360./step
    regular=abs(nominal-round(nominal))<.001 and all(abs((p/step)-round(p/step))<1e-3 for p in unique)
    expected=int(round(nominal)) if regular else None
    present=set(round(p/step)%expected for p in unique) if regular else set()
    missing=[round(i*step,5) for i in range(expected) if i not in present] if regular else []
    dup=sorted({p for p in phi if phi.count(p)>1})
    continuous=all(x.get('mode','').upper()=="CONTINUOUS" for x in rows)
    stationary=all(x.get('mode','').upper()!="CONTINUOUS" and abs(float(x.get('sweep_deg',0)))<1e-6 for x in rows)
    mode='continuous' if continuous else 'stationary' if stationary else 'mixed/unknown'
    sweep=np.median([float(x.get('sweep_deg',0)) for x in rows])
    # For step-and-shoot, we only know the sampled grid, not continuous azimuthal coverage.
    full_turn=bool(regular and not missing)
    complete=bool(full_turn and mode!='mixed/unknown')
    return dict(n_frames=len(rows),unique_frames=len(unique),phi_step_deg=step,
      expected_unique_frames=expected,missing_phi_deg=missing,
      duplicate_canonical_angles=dup,mode=mode,observed_sweep_deg=float(sweep),
      full_turn=full_turn,complete_grid=complete,phi_start_deg=min(unique),
      phi_max_deg=max(unique),sampled_positions_deg=unique)

RASX_FRAME_NAME=re.compile(
    r'^(?P<sample>.+)_image_(?P<seq>\d+)_(?P<chi>-?\d+(?:\.\d+)?)_(?P<frame>\d+)(?:_(?P<point>\d+))?$',
    flags=re.I)

def select_point_groups(rows):
    """Separate XY points even when SmartLab stores all RASX frames in one folder.

    The trailing _01, _02, ... after the phi-frame index is an optional
    acquisition/position suffix. It is *not* interpreted as a measured angle.
    Actual X/Y coordinates are independently retained in the group key.
    """
    groups=collections.defaultdict(list)
    for x in rows:
        p=pathlib.Path(x.get('path','unknown'))
        xy=(round(float(x['x']),3) if x.get('x') is not None else None,
            round(float(x['y']),3) if x.get('y') is not None else None)
        match=RASX_FRAME_NAME.match(p.stem)
        sample=match.group('sample') if match else p.stem
        position_suffix=match.group('point') if match else None
        key=(str(p.parent),sample,xy,position_suffix)
        groups[key].append(x)
    return dict(groups)

def discover(root, progress=None):
    """Read lightweight RASX headers; optionally report each discovery/validation step.

    The callback receives immutable snapshots. It is intended for a background
    worker and must not access Tk widgets. Frame image binaries are not loaded.
    """
    root=pathlib.Path(root)
    if not root.is_dir():raise FileNotFoundError(root)
    state=dict(stage='enumerate',found=0,total=None,checked=0,accepted=0,
               rejected=0,current='',last_ok='',last_error='')
    if progress:progress(state.copy())
    files=[]
    for p in root.rglob('*.rasx'):
        files.append(p)
        state.update(found=len(files),current=p.name)
        if progress:progress(state.copy())
    state.update(stage='sort',total=len(files),current='')
    if progress:progress(state.copy())
    files.sort()
    state.update(stage='headers')
    if progress:progress(state.copy())
    frames=[]
    failures=[]
    for i,p in enumerate(files,1):
        state.update(current=p.name)
        if progress:progress(state.copy())
        try:
            row=parse_rasx_header(p)
            row['path']=str(p)
            frames.append(row)
            state['last_ok']=p.name
        except Exception as e:
            failures.append(dict(file=str(p),error=str(e)))
            state['last_error']=p.name
        state.update(checked=i,accepted=len(frames),rejected=len(failures))
        if progress:progress(state.copy())
    state.update(stage='group',current='')
    if progress:progress(state.copy())
    groups=select_point_groups(frames)
    if progress:progress(dict(state,stage='group_done',groups_total=len(groups)))
    return groups,failures

def validate_group(rows):
    bychi=collections.defaultdict(list)
    for row in rows:bychi[round(row['chi'],4)].append(row)
    if len(bychi)<1:raise ValueError("No chi series")
    info={}
    for chi,items in sorted(bychi.items()):
        v=infer_series(items)
        if not v['complete_grid']:raise ValueError(f"Chi {chi}: incomplete/mixed phi series: {v['missing_phi_deg'][:20]}")
        if v['duplicate_canonical_angles']:
            # Endpoint 360 may duplicate 0: cannot silently double count.
            unique={}
            for row in items:unique.setdefault(row['phi'],row)
            items=list(unique.values())
            bychi[chi]=items
        if v['mode']=='continuous' and abs(v['observed_sweep_deg']-v['phi_step_deg'])>0.05:
            raise ValueError(f"Chi {chi}: sweep {v['observed_sweep_deg']} != step {v['phi_step_deg']}, possible angular gaps/overlaps")
        if v['mode']=='mixed/unknown':raise ValueError("Mixed/unknown scan mode; needs manual review")
        info[str(chi)]=v
    allrows=[x for grp in bychi.values() for x in grp]
    for key,tol in [('two_theta',.02),('omega',.02),('detector_distance',.2)]:
        vals=[x.get(key) for x in allrows]
        if any(v is None for v in vals) or max(vals)-min(vals)>tol:
            raise ValueError(f"Inconsistent {key}; cannot calibrate map")
    if abs(allrows[0]['two_theta']-32.6)>.5 or abs(allrows[0]['omega']-16.3)>.5:
        raise ValueError("Unvalidated reflection geometry; only calibrated LMO (011) near 32.6/16.3 supported")
    return bychi,info

def geometry(ch,ref):
    """Pixel-to-pole transform from verified HyPix H geometry, calibrated layout."""
    U,V=np.meshgrid(np.arange(775),np.arange(385))
    origin=ref['detector_origin']
    dist=ref['detector_distance']
    if abs(origin[0]-382.5)>1 or abs(origin[1]-197.5)>1 or abs(origin[2]-.1)>.003 or abs(origin[3]-.1)>.003 or abs(dist-150.608)>2:
        raise ValueError("Detector calibration differs from validated setup; cannot project without calibration")
    tt=np.deg2rad(ref['two_theta']); omega=np.deg2rad(ref['omega'])
    du=(U-origin[0])*origin[2];dv=(V-origin[1])*origin[3]
    r=np.stack((-du,-dv*np.cos(tt)-dist*np.sin(tt),dv*np.sin(tt)-dist*np.cos(tt)),axis=-1)
    s=r/np.linalg.norm(r,axis=-1,keepdims=True)
    twotheta=np.rad2deg(np.arccos(np.clip(-s[...,2],-1,1)))
    q=s+np.array([0.,0.,1.]);q/=np.linalg.norm(q,axis=-1,keepdims=True)
    q0=np.stack((q[...,0],np.cos(omega)*q[...,1]-np.sin(omega)*q[...,2],
                 np.sin(omega)*q[...,1]+np.cos(omega)*q[...,2]),axis=-1)
    al=np.deg2rad(90-ch)
    qs=np.stack((q0[...,0]*np.cos(al)+q0[...,1]*np.sin(al),-q0[...,0]*np.sin(al)+q0[...,1]*np.cos(al),q0[...,2]),axis=-1)
    qs=np.where((qs[...,0]>0)[...,None],-qs,qs)
    tilt=np.rad2deg(np.arccos(np.clip(-qs[...,0],-1,1)))
    azi=np.rad2deg(np.arctan2(qs[...,2],-qs[...,1]))%360
    signal=(twotheta>=31.6)&(twotheta<33.3)&(tilt>=20)&(tilt<65)
    bg1=(twotheta>=30.8)&(twotheta<31.3)&(tilt>=20)&(tilt<65)
    bg2=(twotheta>=33.7)&(twotheta<34.2)&(tilt>=20)&(tilt<65)
    return tilt,azi,signal,bg1,bg2

def _image(p):
    with zipfile.ZipFile(p) as z:raw=z.read('Data0/Image0.bin')
    if len(raw)!=775*385*4:raise ValueError(f"Unexpected image size: {p}")
    return np.frombuffer(raw,dtype='<u4').reshape(385,775)

def reconstruct(rows,output):
    bychi,info=validate_group(rows)
    sums=np.zeros((90,180),float); coverage=np.zeros_like(sums)
    # Stationary exposures are point samples. Continuous exposures cover true angular intervals.
    # We use an explicit quadrature over each frame's recorded angle range.
    for ch,series in sorted(bychi.items()):
        ref=series[0]
        tilt,az,signal,bg1,bg2=geometry(ch,ref)
        idx=np.flatnonzero(signal)
        if not len(idx):raise ValueError(f"No LMO(011) pixel coverage at chi={ch}")
        # Apply background correction using counts-per-pixel in neighbouring radial bands.
        masks=[signal,bg1,bg2]
        norm=[int(mask.sum()) for mask in masks]
        if not all(norm):raise ValueError("No background/signal pixel ROI")
        selected=np.flatnonzero(signal)
        t=np.floor(tilt.flat[selected]).astype(int)
        a=az.flat[selected]
        s_factors=np.array([np.mean(tilt[m]) for m in masks])
        # Estimate background linear in two-theta using actual local detector angular means.
        # 2theta average for each annulus is derived from the calibrated detector, not the tilt.
        theta=ref['two_theta']
        for j,row in enumerate(series):
            image=_image(row['path'])
            # Preserve signal shape; subtract a scalar background per pixel for this frame.
            bg=(float(np.mean(image[bg1]))+float(np.mean(image[bg2])))/2
            pixels=image.flat[selected].astype(float)-bg
            phi=float(row['phi']); sweep=float(row.get('sweep_deg',0))
            if row['mode']=='CONTINUOUS' and sweep>0:
                lo=float(row['scan_start'])
                hi=float(row['scan_stop'])
                phis=np.linspace(lo,hi,5)
            else:phis=[phi]
            for angle in phis:
                b=np.floor(((a+angle+1)%360)/2).astype(int)
                linear=t*180+b
                sums+=np.bincount(linear,weights=pixels/len(phis),minlength=90*180).reshape(90,180)
                coverage+=np.bincount(linear,weights=np.ones(len(linear))/len(phis),minlength=90*180).reshape(90,180)
    data=np.divide(sums,coverage,out=np.full_like(sums,np.nan),where=coverage>0)
    np.savez_compressed(output/'map.npz',intensity=data,coverage=coverage,
                        tilt=np.arange(90)+.5,azimuth=np.arange(180)*2+1)
    return data,coverage,info

_rich_cache=None

def rich_reports():
    global _rich_cache
    if _rich_cache is None:
        script=pathlib.Path(__file__).with_name("rich_reports.py")
        if not script.is_file():raise FileNotFoundError(f"Rich reports module missing: {script}")
        _rich_cache=runpy.run_path(str(script),run_name="_rich_reports_")
    return _rich_cache

def report(out,rows,data,cov,info,frame_count=None,original_acquisition=None):
    """Independent 2D Gaussian fits; all reports use exactly these parameters."""
    frame_count=len(rows) if frame_count is None else int(frame_count)
    acquisition=(dict(original_acquisition) if original_acquisition else dict(
        n_frames=frame_count,per_chi=info,
        xy=dict(x=rows[0].get("x"),y=rows[0].get("y")),
        method="Acquisition-aware pixel projection; measured Phi/Chi; five-point quadrature for continuous sweep",
        limitations="Validated LMO(011) HyPix horizontal only. No instrumental or tape-axis calibration."))
    acquisition["n_frames"]=frame_count
    with (out/"acquisition.json").open("w",encoding="utf-8") as f:
        json.dump(acquisition,f,indent=2,ensure_ascii=False)
    peaks=[]
    tilt,az=np.meshgrid(np.arange(90)+.5,np.arange(180)*2+1,indexing="ij")
    for seed in (82.,177.,272.,357.):
        beta=((az-seed+180.)%360.)-180.
        mask=np.isfinite(data)&(tilt>=30)&(tilt<=60)&(np.abs(beta)<25)&(cov>0)
        if mask.sum()<12:
            peaks.append(dict(seed=seed,status="INSUFFICIENT_COVERAGE"));continue
        x=beta[mask];y=tilt[mask];z=data[mask]
        med=float(np.nanmedian(z))
        peakid=int(np.nanargmax(z))
        def residual(p):
            amp,b0,t0,sb,st,bg=p
            model=amp*np.exp(-.5*((x-b0)/sb)**2-.5*((y-t0)/st)**2)+bg
            return model-z
        initial=[max(float(z.max()-med),1.),float(x[peakid]),float(y[peakid]),5.,3.,med]
        try:
            q=least_squares(residual,initial,
                bounds=([0.,-20.,25.,.5,.5,-np.inf],[np.inf,20.,65.,25.,20.,np.inf]),
                max_nfev=500)
            amp,b0,t0,sb,st,bg=[float(v) for v in q.x]
            err=residual(q.x)
            sse=float(np.sum(err**2))
            total=float(np.sum((z-np.mean(z))**2))
            r2=float(1.-sse/total) if total>1e-14 else None
            factor=2.*math.sqrt(2.*math.log(2.))
            peaks.append(dict(
                seed=float(seed),status="OK" if q.success else "FIT_WARNING",
                beta_deg=float((seed+b0)%360.),tilt_deg=t0,
                phi_fwhm_observed_deg=sb*factor,
                tilt_fwhm_observed_deg=st*factor,
                sigma_beta_deg=sb,sigma_tilt_deg=st,
                amplitude=amp,background=bg,n_bins=int(mask.sum()),
                fit_r2=r2,fit_rmse=math.sqrt(sse/len(z))))
        except Exception as e:
            peaks.append(dict(seed=float(seed),status="FIT_FAILED",error=str(e)))
    with (out/"peaks.json").open("w",encoding="utf-8") as f:
        json.dump(peaks,f,indent=2,ensure_ascii=False)
    with (out/"peaks.csv").open("w",newline="",encoding="utf-8-sig") as f:
        headers=sorted(set().union(*(p.keys() for p in peaks)))
        w=csv.DictWriter(f,fieldnames=headers);w.writeheader();w.writerows(peaks)
    rich_reports()["write_point_report"](out,rows,data,cov,info,peaks,frame_count)
    return peaks

def write_batch_summary(out,records):
    """Numeric IDs, four polar coordinates and widths, cone tilt and visual maps."""
    return rich_reports()["write_batch_report"](out,records)

def run(source,output,inspect=False,group_key=None):
    out=pathlib.Path(output)
    if out.exists():raise FileExistsError(f"Output must not exist: {out}")
    groups,failures=discover(source)
    if not groups:raise ValueError("No supported RASX frames found")
    if group_key is not None:
        groups={k:v for k,v in groups.items() if str(k)==group_key}
        if not groups:raise ValueError(f"Selected series not found: {group_key}")
    # Process from read-only inputs into brand-new directory only.
    out.mkdir(parents=True)
    records=[]
    for i,(key,rows) in enumerate(groups.items(),1):
        item=dict(point=i,source_group=str(key),frames=len(rows),
                  position_suffix=(int(key[3]) if key[3] is not None else None),x_mm=key[2][0],y_mm=key[2][1])
        try:
            series,info=validate_group(rows)
            item['acquisition']=info
            if not inspect:
                dest=out/f'point_{i:03d}'
                dest.mkdir()
                data,cov,_=reconstruct(rows,dest)
                peaks=report(dest,rows,data,cov,info)
                item['peaks']=peaks
            item['status']='OK'
        except Exception as e:item.update(status='FAILED',reason=str(e))
        records.append(item)
    with (out/'inspection.json').open('w') as f:json.dump(dict(groups=records,rejected_files=failures),f,indent=2)
    with (out/'inspection.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['point','position_suffix','x_mm','y_mm','source_group','frames','status','reason'],extrasaction='ignore')
        w.writeheader();w.writerows(records)
    if not inspect:write_batch_summary(out,records)
    for item in records:print(item['point'],item['status'],item.get('reason',''),item['frames'])
    return 0 if all(q['status']=='OK' for q in records) else 2

def rebuild_reports_only(source,output):
    """Replot existing map.npz and acquisition.json WITHOUT opening raw RASX.

    Regenerates all Gaussian fits from the measured reconstructed 2D maps.
    Output must not exist and MUST be separate from the existing result folder.
    """
    source=pathlib.Path(source)
    out=pathlib.Path(output)
    if out.exists():raise FileExistsError(f"Output must not exist: {out}")
    meta=source/"inspection.json"
    if not meta.is_file():raise FileNotFoundError(f"Existing batch inspection.json not found: {meta}")
    info=json.loads(meta.read_text(encoding="utf-8"))
    original=info.get("groups",[])
    if not original:raise ValueError("No recorded groups in inspection.json")
    out.mkdir(parents=True)
    rows=[]
    n_groups=len(original)
    for ordinal,old in enumerate(original,1):
        item=dict(old)
        item["point"]=int(old.get("point",ordinal))
        point=int(item["point"])
        suffix=item.get("position_suffix")
        item["position_suffix"]=int(suffix) if suffix not in (None,"") else None
        item["status"]=old.get("status","FAILED")
        if item["status"]!="OK":
            rows.append(item)
            print(f"REPORT {ordinal}/{n_groups}: point {point} unavailable",flush=True)
            continue
        src=source/f"point_{point:03d}"
        dest=out/f"point_{point:03d}"
        try:
            if not (src/"map.npz").is_file():
                raise FileNotFoundError(f"Saved map missing: {src/'map.npz'}")
            acq=json.loads((src/"acquisition.json").read_text(encoding="utf-8"))
            with np.load(src/"map.npz",allow_pickle=False) as z:
                data=z["intensity"].copy()
                cov=z["coverage"].copy()
            if data.shape!=(90,180) or cov.shape!=data.shape:
                raise ValueError("Unsupported saved map shape; raw data left untouched")
            dest.mkdir()
            shutil.copy2(src/"map.npz",dest/"map.npz")
            n_frames=int(acq.get("n_frames",item.get("frames",0)))
            x=item.get("x_mm",acq.get("xy",{}).get("x"))
            y=item.get("y_mm",acq.get("xy",{}).get("y"))
            dummy=[dict(x=x,y=y)]
            p=report(dest,dummy,data,cov,acq.get("per_chi",{}),
                     frame_count=n_frames,original_acquisition=acq)
            item["peaks"]=p
            item["status"]="OK"
        except Exception as e:
            item.update(status="FAILED",reason=f"Reports-only rebuild: {e}")
        rows.append(item)
        print(f"REPORT {ordinal}/{n_groups}: point {point} {item['status']}",flush=True)
    with (out/"inspection.json").open("w",encoding="utf-8") as f:
        json.dump(dict(groups=rows,rejected_files=info.get("rejected_files",[]),
                       reports_only_source=str(source)),f,indent=2,ensure_ascii=False)
    with (out/"inspection.csv").open("w",newline="",encoding="utf-8-sig") as f:
        cols=["point","position_suffix","x_mm","y_mm","frames","status","reason"]
        w=csv.DictWriter(f,fieldnames=cols,extrasaction="ignore")
        w.writeheader();w.writerows(rows)
    write_batch_summary(out,rows)
    bad=sum(r.get("status")!="OK" for r in rows)
    print(f"REPORTS_READY: {out}, OK={len(rows)-bad}, FAILED={bad}",flush=True)
    return 0 if bad==0 else 2

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--inspect',action='store_true')
    p.add_argument('--group-key',default=None)
    args=p.parse_args()
    if args.reports_only:raise SystemExit(rebuild_reports_only(args.source,args.output))
    raise SystemExit(run(args.source,args.output,args.inspect,args.group_key))

if __name__=='__main__':main()
