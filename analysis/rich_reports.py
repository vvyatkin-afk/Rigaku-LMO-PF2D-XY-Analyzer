"""Rich, reproducible PDF and CSV reports from already reconstructed LMO(011) maps.

No hardware access. 2D peak widths are observed, NOT instrument deconvolved.
Pole-axis tilt is a conditional four-pole cone fit, not a calibrated tape-axis tilt.
"""
from __future__ import annotations
import csv
import datetime as dt
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse

COLORS=("#2286bc","#8b55ac","#e98a35","#009b88")
DARK="#18334d"
MUTED="#587085"
STEREO="Stereographic radius r = tan(t / 2), dimensionless"
NOTE="Observed peak widths include instrumental broadening and angular sampling; beta is the instrument azimuth reference."

def _number(value, ndigits=3):
    try:
        n=float(value)
        if math.isfinite(n):return round(n,ndigits)
    except (TypeError,ValueError,OverflowError):
        pass
    return None

def numeric_position(value):
    """Produce genuine Excel-compatible integer cells rather than '01' strings."""
    if value is None or str(value).strip()=="":return None
    try:return int(str(value),10)
    except (TypeError,ValueError):return None

def stereo_radius(t_deg):
    try:return math.tan(math.radians(float(t_deg)/2))
    except (ValueError,TypeError):return None

def estimate_structure_tilt(peaks):
    """Fit a common cone axis through four measured scattering-vector directions.

    This only yields the apparent tilt relative to the measurement normal
    *if* the four LMO(011) poles represent an equivalent cone, with reliable
    center fits. Tape-longitudinal orientation is NOT calibrated.
    """
    vals=[]
    for p in peaks[:4]:
        if p.get("status")!="OK" or p.get("beta_deg") is None or p.get("tilt_deg") is None:
            return dict(status="INSUFFICIENT_PEAKS",n_poles=len(vals))
        b=math.radians(float(p["beta_deg"]))
        t=math.radians(float(p["tilt_deg"]))
        vals.append([math.sin(t)*math.cos(b),math.sin(t)*math.sin(b),math.cos(t)])
    if len(vals)!=4:return dict(status="INSUFFICIENT_PEAKS",n_poles=len(vals))
    q=np.asarray(vals,float)
    if not np.isfinite(q).all():return dict(status="NONFINITE_CENTERS",n_poles=4)
    centered=q-q.mean(axis=0)
    _,svd_s,vt=np.linalg.svd(centered,full_matrices=False)
    n=vt[-1,:]
    if n[2]<0:n=-n
    # A collapsed set of nearly identical poles has an indeterminate axis.
    if svd_s[1]<1e-5:return dict(status="DEGENERATE_POLES",n_poles=4)
    t_model=np.rad2deg(np.arccos(np.clip(q@n,-1,1)))
    axis_tilt=math.degrees(math.acos(float(np.clip(n[2],-1,1))))
    az=(math.degrees(math.atan2(float(n[1]),float(n[0])))+360.)%360.
    return dict(status="OK",n_poles=4,
                axis_tilt_deg=float(axis_tilt),
                axis_azimuth_deg=float(az),
                cone_halfangle_deg=float(np.mean(t_model)),
                cone_rms_deg=float(np.std(t_model)),
                method="four equivalent poles: least-squares plane/cone normal",
                azimuth_reference="instrument beta; tape longitudinal direction not calibrated")

def _head(fig, title,subtitle=None):
    fig.patch.set_facecolor("white")
    fig.text(.065,.965,title,fontsize=18,fontweight="bold",color=DARK,va="top")
    if subtitle:fig.text(.065,.925,subtitle,fontsize=9,color=MUTED,va="top")

def _foot(fig,text=NOTE):
    fig.text(.065,.025,text,fontsize=8,color=MUTED,va="bottom")

def peak_profiles(data,cov,p):
    """Cuts through the actual fitted 2D Gaussian; NO independent 1D refit.

    The modeled value is averaged over exactly the same sampled map cells
    as the observed cut. Missing map bins never receive fabricated values.
    """
    if p.get("status") not in ("OK","FIT_WARNING") or p.get("sigma_beta_deg") is None:
        return {}
    amp=float(p["amplitude"])
    bg=float(p["background"])
    b0=float(p["beta_deg"]);t0=float(p["tilt_deg"])
    sb=float(p["sigma_beta_deg"]);st=float(p["sigma_tilt_deg"])
    az=np.arange(180,dtype=float)*2.+1.
    tilt=np.arange(90,dtype=float)+.5
    beta=np.asarray((az-b0+180.)%360.-180.)
    yy,xx=np.meshgrid(tilt,beta,indexing="ij")
    pred=bg+amp*np.exp(-.5*((xx/sb)**2+(yy-t0)**2/st**2))
    ok=np.isfinite(data)&(cov>0)
    results={}
    # Cut +/-1.5 deg about fitted t center; averaged across sampled t bins.
    tsel=np.where(np.abs(tilt-t0)<=1.6)[0]
    if not len(tsel):tsel=np.array([int(np.argmin(abs(tilt-t0)))])
    beta_cols=np.where(abs(beta)<=24)[0]
    for axis in ("beta","tilt"):
        output=[]
        if axis=="beta":
            for j in beta_cols:
                rows_idx=tsel[ok[tsel,j]]
                if not len(rows_idx):continue
                output.append((float(beta[j]),float(np.mean(data[rows_idx,j])),
                               float(np.mean(pred[rows_idx,j])),int(len(rows_idx))))
        else:
            bsel=np.where(abs(beta)<=3.1)[0]
            if not len(bsel):bsel=np.array([int(np.argmin(abs(beta)))])
            for i,t in enumerate(tilt):
                if t<30 or t>60:continue
                cols=bsel[ok[i,bsel]]
                if not len(cols):continue
                output.append((float(t),float(np.mean(data[i,cols])),
                               float(np.mean(pred[i,cols])),int(len(cols))))
        output.sort(key=lambda x:x[0])
        results[axis]=output
    return results

def _pole_page(pdf,rows,data,cov,peaks,n_frames):
    fig=plt.figure(figsize=(8.5,10.5))
    _head(fig,"LMO (011)  |  measured 2D pole figure",
          f"X = {rows[0].get('x')} mm     Y = {rows[0].get('y')} mm     {n_frames} frames     4 independent peak fits")
    ax=fig.add_axes([.11,.17,.70,.67],projection="polar")
    phi=np.deg2rad(np.arange(180)*2.+1.)
    r=np.tan(np.deg2rad(np.arange(90)+.5)/2)
    a,rr=np.meshgrid(phi,r)
    im=ax.pcolormesh(a,rr,np.ma.masked_invalid(data),shading="nearest",cmap="magma",rasterized=True)
    ax.set_ylim(0,math.tan(math.radians(70)/2))
    ax.set_thetagrids(np.arange(0,360,45))
    ax.tick_params(labelsize=9)
    cb=fig.colorbar(im,ax=ax,shrink=.68,pad=.09)
    cb.set_label("Background-subtracted detector counts / contributing pixel",fontsize=8)
    for i,p in enumerate(peaks[:4]):
        if p.get("status") not in ("OK","FIT_WARNING"):continue
        theta=math.radians(float(p["beta_deg"]))
        radius=stereo_radius(p["tilt_deg"])
        ax.plot(theta,radius,marker="o",ms=8,mfc="none",mec="white",mew=1.5)
        ax.annotate(f"P{i+1}",xy=(theta,radius),xytext=(10,8),
                    textcoords="offset points",color="white",fontsize=9,fontweight="bold")
    fig.text(.12,.12,"White circles: centers of the independent 2D Gaussian fits.",fontsize=10,color=DARK)
    fig.text(.12,.092,STEREO,fontsize=9,color=MUTED)
    _foot(fig)
    pdf.savefig(fig,dpi=180);plt.close(fig)

def _parameters_page(pdf,rows,peaks,info,axis,n_frames):
    fig=plt.figure(figsize=(11.7,8.3))
    _head(fig,"Four-pole measurement and structural tilt",
          f"Stage X={rows[0].get('x')} mm  |  Y={rows[0].get('y')} mm  |  {n_frames} HyPix images")
    axs=fig.add_axes([.065,.47,.87,.35]);axs.axis("off")
    headers=["Pole","Azimuth beta\n(deg)","Radius t\n(deg)","Stereo r\n(unitless)","FWHM beta\n(deg)","FWHM t\n(deg)","2D fit R2","Status"]
    cells=[]
    for i,p in enumerate(peaks[:4],1):
        def f(k,fmt=".2f"):
            v=p.get(k)
            return "--" if v is None else format(float(v),fmt)
        cells.append([f"P{i}",f("beta_deg"),f("tilt_deg"),
                      "--" if p.get("tilt_deg") is None else f"{stereo_radius(p['tilt_deg']):.4f}",
                      f("phi_fwhm_observed_deg"),f("tilt_fwhm_observed_deg"),
                      f("fit_r2",".3f"),p.get("status","--")])
    tab=axs.table(cellText=cells,colLabels=headers,cellLoc="center",loc="center",
                   colWidths=[.10,.125,.12,.12,.13,.12,.11,.17])
    tab.auto_set_font_size(False);tab.set_fontsize(9);tab.scale(1,1.8)
    for (r,c),cell in tab.get_celld().items():
        cell.set_edgecolor("#dfe8ef")
        if r==0:
            cell.set_facecolor(DARK);cell.get_text().set_color("white");cell.get_text().set_weight("bold")
        elif r%2==1:cell.set_facecolor("#eff5f9")
        if r>0 and c==0:cell.get_text().set_color(COLORS[r-1]);cell.get_text().set_weight("bold")
    if axis.get("status")=="OK":
        fig.text(.077,.412,
            f"Conditional cone-axis tilt: {axis['axis_tilt_deg']:.3f} deg    |    "
            f"Direction beta: {axis['axis_azimuth_deg']:.1f} deg",
            fontsize=12,color=DARK,fontweight="bold")
        fig.text(.077,.386,
            f"Cone half-angle: {axis['cone_halfangle_deg']:.3f} deg    |    "
            f"Fit RMS: {axis['cone_rms_deg']:.3f} deg",
            fontsize=11,color=DARK,fontweight="bold")
    else:
        fig.text(.077,.405,"Conditional cone-axis tilt: unavailable ("+
                 axis.get("status","unknown")+")",fontsize=12,color=DARK,fontweight="bold")
    fig.text(.077,.350,"Tilt is estimated from the 3D directions of all four pole centers, not from one radius or FWHM.",
             fontsize=10,color=MUTED)
    fig.text(.077,.327,"Valid only if the four reflections belong to an equivalent cone; absolute tape tilt needs axis calibration.",
             fontsize=9.3,color=MUTED)
    fig.text(.077,.292,"Actual acquisition per chi (from RASX):",fontsize=11,fontweight="bold",color=DARK)
    for j,(chi,v) in enumerate(sorted(info.items(),key=lambda kv:float(kv[0]))):
        fig.text(.09,.255-j*.032,f"Chi {chi} deg: {v.get('n_frames')} frames; {v.get('mode')}; "
                 f"phi step {v.get('phi_step_deg')} deg; sweep {v.get('observed_sweep_deg')} deg",
                 fontsize=9.8,color=DARK)
    fig.text(.077,.105,"FWHM values are observed 2D Gaussian fit widths, not resolution-corrected mosaic spread.",
             fontsize=9.4,color=DARK)
    _foot(fig)
    pdf.savefig(fig);plt.close(fig)

def _fit_curves_page(pdf,data,cov,peaks,out):
    fig,axs=plt.subplots(4,2,figsize=(11.7,11),constrained_layout=False)
    fig.subplots_adjust(left=.08,right=.96,top=.88,bottom=.09,wspace=.27,hspace=.52)
    _head(fig,"Measured peak cuts and the fitted model",
          "Points = measured map intensities; continuous lines = cuts through the actual 2D Gaussian fit.")
    cut_rows=[]
    for i,p in enumerate(peaks[:4]):
        prof=peak_profiles(data,cov,p)
        for j,axisname in enumerate(("beta","tilt")):
            ax=axs[i,j]
            vals=prof.get(axisname,[])
            if not vals:
                ax.text(.5,.5,"No calibrated fit / coverage",ha="center",va="center",transform=ax.transAxes)
                ax.set_axis_off();continue
            x=np.array([v[0] for v in vals]);ym=np.array([v[1] for v in vals]);yf=np.array([v[2] for v in vals])
            ax.plot(x,yf,color=COLORS[i],lw=2.5,label="same 2D fit")
            ax.plot(x,ym,color=DARK,ls="none",marker="o",ms=3.1,alpha=.75,label="measured")
            if axisname=="beta":
                half=p["phi_fwhm_observed_deg"]/2
                center=0.
                xlabel="Azimuth offset from fitted beta (deg)"
                ttl=f"P{i+1}: beta {p['beta_deg']:.2f} deg  |  FWHM {2*half:.2f} deg"
            else:
                half=p["tilt_fwhm_observed_deg"]/2
                center=float(p["tilt_deg"])
                xlabel="Radial position t (deg)"
                ttl=f"P{i+1}: t {center:.2f} deg  |  FWHM {2*half:.2f} deg"
            ax.axvline(center,color=COLORS[i],lw=.9,alpha=.6,ls=":")
            ax.axvline(center-half,color=COLORS[i],lw=.9,alpha=.55,ls="--")
            ax.axvline(center+half,color=COLORS[i],lw=.9,alpha=.55,ls="--")
            ax.set(title=ttl,xlabel=xlabel,ylabel="Counts / pixel")
            ax.title.set_fontsize(10)
            ax.xaxis.label.set_fontsize(9);ax.yaxis.label.set_fontsize(9)
            ax.tick_params(labelsize=8)
            ax.grid(alpha=.15)
            if i==0:ax.legend(fontsize=8,loc="upper right")
            for x0,y0,z0,n in vals:
                cut_rows.append(dict(peak=i+1,axis=axisname,coordinate_deg=x0,
                                     measured=y0,model=z0,n_map_bins=n))
    fig.text(.08,.047,"Model is evaluated at the original fitted parameters and averaged across the same sampled map bins.",
             fontsize=8.7,color=MUTED)
    pdf.savefig(fig);plt.close(fig)
    with (out/"peak_profile_cuts.csv").open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=["peak","axis","coordinate_deg","measured","model","n_map_bins"])
        w.writeheader();w.writerows(cut_rows)

def _peak_heatmaps_page(pdf,data,cov,peaks):
    fig,axs=plt.subplots(2,2,figsize=(11.7,9))
    fig.subplots_adjust(left=.08,right=.91,top=.88,bottom=.1,hspace=.30,wspace=.32)
    _head(fig,"2D measured intensity and fitted FWHM contours",
          "Dashed ellipse = full width at half maximum of the SAME 2D Gaussian fit, at each peak center.")
    az=np.arange(180)*2.+1.; tt=np.arange(90)+.5
    for i,p in enumerate(peaks[:4]):
        ax=axs[i//2,i%2]
        if p.get("status") not in ("OK","FIT_WARNING"):
            ax.text(.5,.5,"Peak unavailable",ha="center",va="center",transform=ax.transAxes);continue
        b=np.asarray((az-float(p["beta_deg"])+180.)%360.-180.)
        inds=np.where(abs(b)<=22)[0]
        inds=inds[np.argsort(b[inds])]
        tid=np.where((tt>=30)&(tt<=60))[0]
        sub=data[np.ix_(tid,inds)]
        local=np.ma.masked_invalid(sub)
        im=ax.pcolormesh(b[inds],tt[tid],local,cmap="magma",shading="nearest",rasterized=True)
        e=Ellipse((0.,float(p["tilt_deg"])),width=p["phi_fwhm_observed_deg"],
                  height=p["tilt_fwhm_observed_deg"],fill=False,edgecolor="cyan",lw=1.8,ls="--")
        ax.add_patch(e)
        ax.plot(0.,float(p["tilt_deg"]),marker="+",ms=12,mew=2,color="white")
        ax.set(xlim=(-23,23),ylim=(30,60),xlabel="Beta offset (deg)",ylabel="Radial t (deg)",
               title=f"P{i+1}  |  fit R2 {p.get('fit_r2',float('nan')):.3f}")
        ax.grid(alpha=.12)
        fig.colorbar(im,ax=ax,fraction=.045,pad=.03)
    _foot(fig,"Cyan ellipses are model half-maximum contours; uncovered pixel bins remain masked.")
    pdf.savefig(fig);plt.close(fig)

def write_point_report(out,rows,data,cov,info,peaks,n_frames):
    """Create a four-page visual report plus a replot-friendly profile CSV."""
    out=Path(out)
    axis=estimate_structure_tilt(peaks)
    with (out/"structure_tilt.json").open("w",encoding="utf-8") as f:
        json.dump(axis,f,indent=2,ensure_ascii=False)
    with PdfPages(out/"LMO011_adaptive_pole_figure.pdf") as pdf:
        _pole_page(pdf,rows,data,cov,peaks,n_frames)
        _parameters_page(pdf,rows,peaks,info,axis,n_frames)
        _fit_curves_page(pdf,data,cov,peaks,out)
        _peak_heatmaps_page(pdf,data,cov,peaks)
    return axis

def summarize_records(records):
    """Flatten four independent peak parameter sets to numeric CSV columns."""
    rows=[]
    for item in records:
        position=numeric_position(item.get("position_suffix"))
        r=dict(point=int(item.get("point",len(rows)+1)),position_suffix=position,
               x_mm=_number(item.get("x_mm")),y_mm=_number(item.get("y_mm")),
               frames=int(item.get("frames",0)),status=item.get("status","FAILED"),
               reason=item.get("reason",""))
        peaks=item.get("peaks") or []
        for i in range(1,5):
            p=peaks[i-1] if len(peaks)>=i else {}
            r.update({f"P{i}_status":p.get("status","MISSING"),
                      f"P{i}_beta_deg":_number(p.get("beta_deg")),
                      f"P{i}_tilt_deg":_number(p.get("tilt_deg")),
                      f"P{i}_stereo_radius":_number(stereo_radius(p.get("tilt_deg")),6),
                      f"P{i}_phi_fwhm_deg":_number(p.get("phi_fwhm_observed_deg")),
                      f"P{i}_tilt_fwhm_deg":_number(p.get("tilt_fwhm_observed_deg")),
                      f"P{i}_fit_r2":_number(p.get("fit_r2"),4)})
        axis=estimate_structure_tilt(peaks)
        r.update(axis_status=axis.get("status","INSUFFICIENT_PEAKS"),
                 axis_tilt_deg=_number(axis.get("axis_tilt_deg")),
                 axis_azimuth_deg=_number(axis.get("axis_azimuth_deg")),
                 cone_halfangle_deg=_number(axis.get("cone_halfangle_deg")),
                 cone_rms_deg=_number(axis.get("cone_rms_deg")))
        phi=[r.get(f"P{i}_phi_fwhm_deg") for i in range(1,5)]
        rad=[r.get(f"P{i}_tilt_fwhm_deg") for i in range(1,5)]
        r["mean_phi_fwhm_deg"]=_number(np.mean([v for v in phi if v is not None])) if any(v is not None for v in phi) else None
        r["mean_tilt_fwhm_deg"]=_number(np.mean([v for v in rad if v is not None])) if any(v is not None for v in rad) else None
        rows.append(r)
    return rows

def _matrix_map(pdf,rows,key_template,title,unit,cmap,global_scale=True):
    fig,axs=plt.subplots(2,2,figsize=(11.7,8.5))
    fig.subplots_adjust(left=.08,right=.95,top=.865,bottom=.12,hspace=.33,wspace=.27)
    _head(fig,title,f"Spatial distribution across measured X/Y; labels are numeric acquisition positions. Unit: {unit}.")
    allvalues=[r.get(key_template.format(i)) for r in rows for i in range(1,5)]
    allvalues=[float(x) for x in allvalues if x is not None and math.isfinite(x)]
    mn=min(allvalues) if allvalues else 0.
    mx=max(allvalues) if allvalues else 1.
    if mn==mx:mx=mn+1.
    for i,ax in enumerate(axs.flat,1):
        key=key_template.format(i)
        good=[r for r in rows if r.get("x_mm") is not None and r.get("y_mm") is not None and r.get(key) is not None]
        if not good:
            ax.text(.5,.5,"No valid fitted peaks",ha="center",va="center",transform=ax.transAxes)
        else:
            xs=[r["x_mm"] for r in good];ys=[r["y_mm"] for r in good]
            vals=[r[key] for r in good]
            kw=dict(c=vals,cmap=cmap,s=130,edgecolors="white",linewidth=.9)
            if global_scale:kw.update(vmin=mn,vmax=mx)
            sc=ax.scatter(xs,ys,**kw)
            for row,x,y in zip(good,xs,ys):
                ax.annotate(str(row["position_suffix"] if row["position_suffix"] is not None else row["point"]),
                            (x,y),xytext=(3,3),textcoords="offset points",
                            fontsize=7.8,color=DARK)
            fig.colorbar(sc,ax=ax,shrink=.8,pad=.02)
        ax.set(xlabel="X stage (mm)",ylabel="Y stage (mm)",title=f"P{i}")
        ax.grid(alpha=.16)
        ax.set_aspect("equal",adjustable="datalim")
    _foot(fig)
    pdf.savefig(fig);plt.close(fig)

def _trend_page(pdf,rows):
    good=[r for r in rows if r["status"]=="OK"]
    fig,axs=plt.subplots(2,1,figsize=(11.7,8.5),sharex=True)
    fig.subplots_adjust(left=.09,right=.95,top=.87,bottom=.12,hspace=.25)
    _head(fig,"All four pole widths at each measurement point",
          "Four independent observed widths are shown. The curves connect acquisition indices, not a spatial interpolation.")
    xx=[r["position_suffix"] if r["position_suffix"] is not None else r["point"] for r in good]
    for i in range(1,5):
        for ax,key in zip(axs,(f"P{i}_phi_fwhm_deg",f"P{i}_tilt_fwhm_deg")):
            vals=[np.nan if r.get(key) is None else r[key] for r in good]
            ax.plot(xx,vals,marker="o",markersize=3,linewidth=1.5,
                    color=COLORS[i-1],label=f"P{i}")
    axs[0].set_ylabel("Azimuth FWHM beta (deg)")
    axs[1].set_ylabel("Radial FWHM t (deg)")
    axs[1].set_xlabel("Measured XY position number")
    for ax in axs:
        ax.grid(alpha=.2)
        ax.legend(ncol=4,fontsize=9)
    _foot(fig)
    pdf.savefig(fig);plt.close(fig)

def _tilt_page(pdf,rows):
    fig,axs=plt.subplots(1,2,figsize=(11.7,8.5))
    fig.subplots_adjust(left=.08,right=.95,top=.87,bottom=.12,wspace=.32)
    _head(fig,"Apparent four-pole cone-axis tilt",
          "Conditional geometric estimate. Axis direction uses instrument beta, not a calibrated tape direction.")
    good=[r for r in rows if r.get("axis_tilt_deg") is not None and
          r.get("x_mm") is not None and r.get("y_mm") is not None]
    if good:
        xs=[r["x_mm"] for r in good];ys=[r["y_mm"] for r in good];zs=[r["axis_tilt_deg"] for r in good]
        sc=axs[0].scatter(xs,ys,c=zs,s=150,cmap="plasma",edgecolors="white",linewidth=1)
        fig.colorbar(sc,ax=axs[0],shrink=.64,label="Axis tilt (deg)")
        for r,x,y in zip(good,xs,ys):
            axs[0].annotate(str(r["position_suffix"] or r["point"]),(x,y),
                            xytext=(4,4),textcoords="offset points",fontsize=8)
        uu=np.sin(np.deg2rad([r["axis_tilt_deg"] for r in good]))*np.cos(np.deg2rad([r["axis_azimuth_deg"] for r in good]))
        vv=np.sin(np.deg2rad([r["axis_tilt_deg"] for r in good]))*np.sin(np.deg2rad([r["axis_azimuth_deg"] for r in good]))
        axs[1].quiver(xs,ys,uu,vv,angles="xy",scale_units="xy",scale=.025,color="#147f7d",width=.006)
    else:
        for ax in axs:ax.text(.5,.5,"Not enough valid four-pole fits",ha="center",va="center",transform=ax.transAxes)
    for ax in axs:
        ax.set(xlabel="X stage (mm)",ylabel="Y stage (mm)")
        ax.grid(alpha=.17)
        ax.set_aspect("equal",adjustable="datalim")
    axs[0].set_title("Axis tilt magnitude")
    axs[1].set_title("Relative tilt directions (un-calibrated beta axis)")
    _foot(fig,"Interpretation requires an independent instrument-to-tape azimuth reference and valid equivalent-pole indexing.")
    pdf.savefig(fig);plt.close(fig)

def _detail_table_pages(pdf,rows,title,columns,widths=None):
    # Portrait-facing, readable landscape tables; repeat headers on every page.
    for start in range(0,len(rows),17):
        part=rows[start:start+17]
        fig=plt.figure(figsize=(11.7,8.3))
        _head(fig,title,f"Measured positions {start+1}-{start+len(part)} of {len(rows)}")
        ax=fig.add_axes([.065,.13,.87,.72]);ax.axis("off")
        titles=[label for label,key in columns]
        data=[]
        for r in part:
            cells=[]
            for _,key in columns:
                val=r.get(key)
                if val is None:cells.append("--")
                elif isinstance(val,(float,np.floating)):cells.append(f"{val:.2f}")
                else:cells.append(str(val))
            data.append(cells)
        t=ax.table(cellText=data,colLabels=titles,cellLoc="center",loc="center",
                   colWidths=widths)
        t.auto_set_font_size(False);t.set_fontsize(9);t.scale(1,1.44)
        for (rr,c),cell in t.get_celld().items():
            cell.set_edgecolor("#d8e2ea")
            if rr==0:
                cell.set_facecolor(DARK);cell.get_text().set_color("white")
                cell.get_text().set_weight("bold")
            elif rr%2==0:cell.set_facecolor("#edf4f8")
            if rr>0 and c>=3 and c<=6:cell.get_text().set_color(COLORS[c-3])
        _foot(fig)
        pdf.savefig(fig);plt.close(fig)

def write_batch_report(out,records):
    """Write numeric CSV/JSON plus full PDF maps, fit parameters and tilt checks."""
    out=Path(out)
    rows=summarize_records(records)
    fields=list(rows[0]) if rows else ["point","position_suffix","x_mm","y_mm","status"]
    with (out/"XY_summary.csv").open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    with (out/"XY_summary.json").open("w",encoding="utf-8") as f:
        json.dump(rows,f,indent=2,ensure_ascii=False)
    with PdfPages(out/"LMO011_XY_Overall_Report.pdf") as pdf:
        fig=plt.figure(figsize=(11.7,8.3))
        _head(fig,"LMO (011) | Batch pole figures and XY texture",
              f"Generated {dt.datetime.now().isoformat(timespec='minutes')} | real measured positions | read-only source data")
        ok=[r for r in rows if r["status"]=="OK"]
        fitok=[r for r in rows if r.get("axis_tilt_deg") is not None]
        fig.text(.08,.79,f"{len(rows)}",fontsize=49,color="#1e7eac",weight="bold")
        fig.text(.27,.80,"measured XY series",fontsize=16,color=DARK)
        fig.text(.08,.67,f"{len(ok)} / {len(rows)}",fontsize=29,color="#148977",weight="bold")
        fig.text(.32,.68,"point maps reconstructed",fontsize=14,color=DARK)
        fig.text(.08,.59,f"{len(fitok)} / {len(rows)}",fontsize=29,color="#a15d3a",weight="bold")
        fig.text(.32,.60,"valid four-pole cone-axis estimates",fontsize=13,color=DARK)
        fig.text(.08,.47,"Each point has four independent fitted beta centers, radial t centers, and both observed widths.",fontsize=11,color=DARK)
        fig.text(.08,.42,"All four values, rather than only their mean, appear in this report and in XY_summary.csv.",fontsize=11,color=DARK)
        fig.text(.08,.36,"The report includes measured/model cut curves and 2D peak-fit diagnostics for each point.",fontsize=11,color=DARK)
        fig.text(.08,.27,"Coordinate convention: t is polar angle (deg); r = tan(t / 2) is the stereographic map radius.",fontsize=10,color=MUTED)
        fig.text(.08,.22,"Caution: axis tilt is a model-dependent geometry estimate, not a calibrated tape-axis direction.",fontsize=10,color=MUTED)
        _foot(fig)
        pdf.savefig(fig);plt.close(fig)
        if rows:
            _matrix_map(pdf,rows,"P{}_tilt_deg","Four radial peak centers","t (degrees)","viridis")
            _matrix_map(pdf,rows,"P{}_phi_fwhm_deg","Four observed azimuthal FWHM","beta width (degrees)","plasma")
            _matrix_map(pdf,rows,"P{}_tilt_fwhm_deg","Four observed radial FWHM","t width (degrees)","magma")
            _matrix_map(pdf,rows,"P{}_stereo_radius","Four stereographic pole radii","dimensionless r","cividis")
            _trend_page(pdf,rows)
            _tilt_page(pdf,rows)
            base=[("Point","point"),("Position","position_suffix"),("X mm","x_mm"),("Y mm","y_mm")]
            quantities=[
                ("Measured azimuth peak centers beta (degrees)","P{}_beta_deg"),
                ("Measured radial pole positions t (degrees)","P{}_tilt_deg"),
                ("Observed azimuthal peak widths FWHM beta (degrees)","P{}_phi_fwhm_deg"),
                ("Observed radial peak widths FWHM t (degrees)","P{}_tilt_fwhm_deg"),
                ("Stereographic pole radii r = tan(t / 2)","P{}_stereo_radius"),
                ("2D Gaussian quality: R2 for four independent peak fits","P{}_fit_r2"),
            ]
            for title,key in quantities:
                cols=base+[(f"P{i}",key.format(i)) for i in range(1,5)]
                _detail_table_pages(pdf,rows,title,cols)
            tilt_cols=base+[("Axis tilt","axis_tilt_deg"),("Axis beta","axis_azimuth_deg"),
                            ("Cone angle","cone_halfangle_deg"),("Cone RMS","cone_rms_deg")]
            _detail_table_pages(pdf,rows,"Conditional four-pole texture-axis tilt (degrees)",tilt_cols)
            failure=[r for r in rows if r["status"]!="OK" or r["axis_status"]!="OK"]
            if failure:
                fig=plt.figure(figsize=(11.7,8.3))
                _head(fig,"Unresolved positions and fit caveats","Rows below were not silently interpolated.")
                for i,r in enumerate(failure[:22]):
                    fig.text(.075,.85-.031*i,
                             f"Point {r['point']} | position {r['position_suffix']} | {r['status']} "
                             f"| axis {r['axis_status']} | {str(r.get('reason') or '')[:95]}",
                             fontsize=9.3,color=DARK)
                _foot(fig)
                pdf.savefig(fig);plt.close(fig)
    return rows
