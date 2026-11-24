"""Append-only lattice and detailed pole-report pages for the universal analyzer."""
from __future__ import annotations
import collections
import itertools
import json
import math
from pathlib import Path
import textwrap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from analysis import lattice_strain as ls
from analysis.rich_reports import estimate_structure_tilt

def peak_profiles(data,cov,p):
    """Generalized cuts at any tilt, consistent with the 2D fitted Gaussian."""
    nt,na=data.shape
    t=np.arange(nt)+.5
    beta=(np.arange(na)+.5)*360./na
    d=(beta-p["beta_deg"]+180)%360-180
    model=p["background"]+p["amplitude"]*np.exp(-.5*((t[:,None]-p["tilt_deg"])/p["sigma_tilt_deg"])**2-.5*(d[None,:]/p["sigma_beta_deg"])**2)
    ok=np.isfinite(data)&(cov>0)
    r={}
    for name,coord,lim,aperture in (
        ("beta",d,24.,np.abs(t-p["tilt_deg"])<=1.6),
        ("tilt",t-p["tilt_deg"],15.,np.abs(d)<=3.1)):
        rows=[]
        for i in np.flatnonzero(np.abs(coord)<=lim):
            if name=="beta":
                sel=aperture&ok[:,i]; y=data[:,i][sel]; m=model[:,i][sel]
            else:
                sel=aperture&ok[i,:]; y=data[i,:][sel]; m=model[i,:][sel]
            if len(y):rows.append((float(coord[i]),float(y.mean()),float(m.mean()),len(y)))
        r[name]=sorted(rows)
    return r

def new_page(title,subtitle):
    fig=plt.figure(figsize=(11.7,8.3),facecolor="white")
    fig.text(.055,.963,title,va="top",fontsize=18,fontweight="bold",color="#18334d")
    fig.text(.055,.914,textwrap.fill(str(subtitle),140),va="top",fontsize=8.5,color="#587085")
    return fig

def finish(pdf,fig):
    pdf.savefig(fig);plt.close(fig)

def fmt(v,n=5):
    return f"{v:.{n}f}" if isinstance(v,(int,float)) and np.isfinite(v) else "N/A"

def table(fig,headers,rows,rect=(.05,.20,.90,.61),widths=None,fontsize=8):
    ax=fig.add_axes(rect);ax.axis("off")
    if not rows:
        ax.text(.5,.5,"No reliable values",ha="center",va="center",transform=ax.transAxes);return
    tab=ax.table(cellText=rows,colLabels=headers,cellLoc="center",loc="center",
        bbox=[0,0,1,1],colWidths=widths)
    tab.auto_set_font_size(False);tab.set_fontsize(fontsize)
    for (i,j),c in tab.get_celld().items():
        c.set_edgecolor("#d9e2eb");c.set_linewidth(.5)
        c.set_facecolor("#264967" if i==0 else ("#eef3f8" if i%2==0 else "#ffffff"))
        if i==0:c.set_text_props(color="white",fontweight="bold")

def foot(fig,lines):
    fig.text(.055,.035,"\n".join(textwrap.fill(str(t),145) for t in lines),va="bottom",
             fontsize=8,color="#52677d",linespacing=1.3)

def cell_plot(ax,A,title):
    vertices={p:A@np.asarray(p,float) for p in itertools.product((0,1),repeat=3)}
    for p,v in vertices.items():
        for i in range(3):
            if p[i]==0:
                q=list(p);q[i]=1;w=vertices[tuple(q)]
                ax.plot(*np.array([v,w]).T,color="#3578a3",lw=1.5)
    origin=np.zeros(3)
    for i,c in enumerate(("#b7524b","#35865c","#7656a0")):
        v=A[:,i]
        ax.quiver(*origin,*v,color=c,arrow_length_ratio=.12)
        ax.text(*v,["a","b","c"][i],color=c)
    pts=np.array(list(vertices.values()))
    ax.set_box_aspect(np.maximum(np.ptp(pts,axis=0),.1))
    ax.set_title(title,fontsize=10)
    ax.set_xlabel("X (A)",fontsize=7);ax.set_ylabel("Y (A)",fontsize=7);ax.set_zlabel("Z (A)",fontsize=7)
    ax.tick_params(labelsize=6)

def cell_pages(pdf,cell,cfg,subtitle):
    fig=new_page("Indexed lattice and strain",subtitle)
    fig.text(.055,.858,("HYPOTHETICAL SINGLE-DOMAIN CELL" if cell.get("hypothetical") else cell.get("status","N/A").replace("_"," ")),fontsize=12,fontweight="bold",color="#254f7c")
    if cell.get("status")!="CONDITIONAL_INDEXED_CELL":
        fig.text(.055,.78,textwrap.fill("A full direct-lattice basis requires at least three independent, "
            "signed and correctly indexed reciprocal vectors from the same single domain. "
            "A reflection family or filename alone does not identify individual spots. "
            "Configure explicit indexing or a confirmed cube-texture orientation in the reference JSON.",115),
            va="top",fontsize=12,linespacing=1.5)
        fig.text(.055,.53,"Usable indexed vectors: "+str(cell.get("n_vectors",0)),fontsize=11)
        foot(fig,[ls.NOTICE]);finish(pdf,fig);return
    sg=cell.get("model_sigma",{})
    rows=[[lab,fmt(cell.get(k)),fmt(sg.get(k)),fmt(cell.get(s))]
          for lab,k,s in (("a (A)","a_A","strain_a_pct"),("b (A)","b_A","strain_b_pct"),
             ("c (A)","c_A","strain_c_pct"),("alpha (deg)","alpha_deg",""),
             ("beta (deg)","beta_deg",""),("gamma (deg)","gamma_deg",""),
             ("Volume (A^3)","volume_A3","volume_change_pct"),("c/a","c_over_a",""),("a/b","a_over_b",""))]
    table(fig,["Quantity","Measured","Model sigma","Change (%)"],rows,rect=(.055,.36,.49,.44),fontsize=8)
    ax=fig.add_axes([.64,.40,.30,.35])
    x=cell["principal_strain_pct"];ax.barh(["e1","e2","e3"],x,color=["#367bb1" if v<0 else "#c76c47" for v in x])
    ax.axvline(0,color="#405364",lw=.8);ax.grid(axis="x",alpha=.2);ax.set_xlabel("Principal stretch strain (%)")
    ax.set_title("Rotation removed",fontsize=10)
    fig.text(.055,.29,f"Rank {cell['rank']}; vectors {cell['n_vectors']}; condition {cell['condition']:.2f}; "
        f"relative vector RMS {100*cell['relative_rms']:.3f}%.",fontsize=10)
    vals=[("Specimen normal",cell.get("normal_strain_pct")),("Tape longitudinal",cell.get("longitudinal_strain_pct")),
          ("Tape transverse",cell.get("transverse_strain_pct"))]
    fig.text(.055,.23,"     ".join(f"{k}: {fmt(v,3)}%" for k,v in vals),fontsize=9)
    foot(fig,[cell.get("warning",""),cell["uncertainty_note"],ls.NOTICE]);finish(pdf,fig)
    fig=new_page("Cell geometry and tensor",subtitle)
    cell_plot(fig.add_axes([.06,.39,.37,.43],projection="3d"),ls.basis(cfg["cell_A_deg"]),"Reference cell (crystal frame)")
    cell_plot(fig.add_axes([.52,.39,.37,.43],projection="3d"),np.asarray(cell["direct_basis_A"]),"Measured cell (specimen frame)")
    E=np.asarray(cell["stretch_strain_specimen_pct"])
    table(fig,["Specimen axes","X","Y","Z"],[[lab,*[fmt(v,4) for v in e]] for lab,e in zip("XYZ",E)],
          rect=(.055,.10,.43,.22),fontsize=8)
    pv=cell["principal_directions_specimen"]
    table(fig,["Principal (%)","X","Y","Z"],[[fmt(e,4),*[fmt(v,4) for v in g]]
        for e,g in zip(cell["principal_strain_pct"],pv)],rect=(.52,.10,.43,.22),fontsize=8)
    foot(fig,["Tensor = left stretch minus identity, expressed in the instrument specimen frame. Principal vector signs are arbitrary.",
              "Axes X/Y are not tape axes unless tape azimuth is independently calibrated."]);finish(pdf,fig)

def point_pages(pdf,dest,item,data,cov,peaks,cfg):
    subtitle=str(item["sample"])+" | group "+str(item["point"])
    # Restore detailed fitted cuts and conditional cone information from v1.3.
    cut_rows=[]
    for p in peaks:
        if p.get("status")!="OK":continue
        cuts=peak_profiles(data,cov,p)
        if not cuts: continue
        for axis,rows in cuts.items():
            cut_rows.extend((p["rank"],axis,*row) for row in rows)
        fig=new_page("Fitted pole cuts",subtitle+f" | peak {p['rank']}")
        for ax,(name,rows) in zip(fig.subplots(1,2),cuts.items()):
            if not rows:continue
            x,y,model,n=np.asarray(rows,float).T
            ax.plot(x,y,"o",ms=3,label="Measured")
            ax.plot(x,model,lw=1.5,label="2D fit")
            ax.legend(fontsize=8)
            ax.set_title(name);ax.set_xlabel("Angle (deg)");ax.set_ylabel("Net rate / pixel");ax.grid(alpha=.2)
        fig.subplots_adjust(left=.09,right=.95,bottom=.15,top=.78,wspace=.30)
        foot(fig,["Cuts through the fitted 2D pole Gaussian, not independent width fits."]);finish(pdf,fig)
    import csv
    with (Path(dest)/"peak_profile_cuts.csv").open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.writer(f);w.writerow(["rank","axis","offset_deg","measured","model","n_bins"]);w.writerows(cut_rows)
    family=item.get("lattice_analysis",{}).get("family_hkl")
    cone=dict(status="NOT_APPLICABLE_OR_INSUFFICIENT_POLES")
    if family and sorted(map(abs,family)) in ([0,1,1],[1,1,1],[2,2,2]) and len(peaks)==4:
        cone=estimate_structure_tilt(peaks)
        if cone.get("status")=="OK":
            fig=new_page("Conditional equivalent-pole cone fit",subtitle)
            fig.text(.07,.75,"\n".join(f"{k.replace('_',' ')}: {fmt(cone.get(k),4)} deg"
                for k in ("axis_tilt_deg","axis_azimuth_deg","cone_halfangle_deg","cone_rms_deg")),fontsize=14,va="top",linespacing=2.0)
            foot(fig,["Only valid under the equivalent-pole cone model. Not a separately calibrated tape tilt.",
                     "A distorted cone may reflect strain, unequal domains, or geometry errors."]);finish(pdf,fig)
    item["structure_tilt"]=cone
    (Path(dest)/"structure_tilt.json").write_text(json.dumps(cone,indent=2),encoding="utf-8")
    r=item.get("lattice_analysis",{})
    pp=r.get("peaks",[])
    fig=new_page("Measured radial diffraction and d-spacing",subtitle)
    fig.text(.055,.854,"Analysis QC: "+r.get("status","N/A").replace("_"," "),fontsize=12,fontweight="bold",color="#254f7c")
    rows=[]
    for p in pp:
        rows.append(["All" if p["rank"]==0 else p["rank"],str(p.get("indexed_hkl","N/A")),fmt(p.get("corrected_two_theta_deg"),4),
            fmt(p.get("d_A")),fmt(p.get("d_fit_sigma_A")),fmt(p.get("reference_d_A")),
            fmt(p.get("normal_strain_pct"),3),fmt(p.get("radial_fwhm_deg"),4),
            p.get("radial_status","N/A").replace("_"," ")])
    table(fig,["Pole","Signed hkl","2theta (deg)","d (A)","Fit sigma d","d ref (A)","Change (%)","FWHM (deg)","Radial QC"],
        rows,rect=(.04,.24,.92,.55),widths=[.055,.12,.11,.105,.11,.105,.10,.105,.19],fontsize=7.4)
    fig.text(.055,.17,f"Material/reference: {cfg['material']} | wavelength {cfg['wavelength_A']:.6f} A | "
        f"2theta correction {cfg['two_theta_offset_deg']:+.4f} deg",fontsize=8.5)
    fig.text(.055,.13,"Filename family: "+str(family)+" | All = aggregate, not an individual scattering vector.",fontsize=8.5)
    foot(fig,[ls.NOTICE,"Fit sigma is statistical model uncertainty, not absolute instrument accuracy."]);finish(pdf,fig)
    if pp:
        fig=new_page("Reciprocal vectors and radial fit details",subtitle)
        rows=[]
        for p in pp:
            g=p.get("g_inv_A",[None]*3)
            rows.append(["All" if p["rank"]==0 else p["rank"],fmt(p.get("observed_two_theta_deg"),4),
                fmt(p.get("reference_two_theta_deg"),4),fmt(p.get("two_theta_fit_sigma_deg"),5),
                fmt(p.get("d_angle_sigma_A"),5),fmt(p.get("radial_fit_r2"),4),
                *[fmt(v,5) for v in g],textwrap.fill(p.get("indexing_status","N/A").replace("_"," "),17)])
        table(fig,["Pole","Raw 2theta","Ref 2theta","Fit sigma deg","Angle sigma d","Radial R2","gx (1/A)","gy (1/A)","gz (1/A)","Indexing QC"],
            rows,rect=(.035,.22,.93,.57),widths=[.05,.10,.10,.105,.105,.09,.085,.085,.085,.195],fontsize=7.5)
        foot(fig,["g = 1/d times the hemisphere-oriented pole unit vector in instrument specimen coordinates.",
            "Raw and reference 2theta in degrees. Angle sigma d includes supplied angular systematic sigma only; N/A if unspecified."])
        finish(pdf,fig)
    f=Path(dest)/"Radial_Profiles.npz"
    if f.is_file():
        with np.load(f,allow_pickle=False) as z:
            x=z["two_theta_deg"];y=z["intensity"];models=z["model"];cc=z["coverage"];ranks=z["rank"]
        for start in range(0,len(ranks),4):
            fig=new_page("Radial intensity profiles",subtitle)
            axs=fig.subplots(2,2).flat
            for k,ax in enumerate(axs):
                i=start+k
                if i>=len(ranks):ax.axis("off");continue
                p=next(q for q in pp if q["rank"]==ranks[i])
                ax.plot(x,y[i],"o",ms=2.2,color="#27567d",label="Measured")
                ax.plot(x,models[i],lw=1.4,color="#d27b39",label="Model")
                ax.set_title(f"{'Aggregate' if ranks[i]==0 else 'Pole '+str(ranks[i])} | {p.get('radial_status','').replace('_',' ')}",fontsize=9)
                ax.set_xlabel("Raw pixel 2theta (deg)",fontsize=8);ax.set_ylabel("Net rate / pixel",fontsize=8)
                ax.tick_params(labelsize=7);ax.grid(alpha=.2);ax.legend(fontsize=7)
            fig.subplots_adjust(left=.09,right=.96,top=.82,bottom=.12,wspace=.30,hspace=.48)
            foot(fig,["Uncovered angle bins remain missing. Configured single/doublet Gaussian + linear background; clipped or poor fits excluded."]);finish(pdf,fig)
    cell_pages(pdf,r.get("lattice",{"status":"UNAVAILABLE"}),cfg,subtitle)

def export_results(out,records,cfg,export_csv):
    """Keep original series separate; only merge explicitly assigned physical points."""
    rows=[];cells=[];tensors=[];groups=collections.defaultdict(list)
    for r in records:
        a=r.get("lattice_analysis",{})
        key=dict(point=r["point"],sample=r["sample"],reflection_hint=r.get("reflection_hint"),
                 x_mm=r.get("x_mm"),y_mm=r.get("y_mm"),position_suffix=r.get("position_suffix"))
        for p in a.get("peaks",[]):rows.append(dict(key,**p))
        cell=a.get("lattice",{})
        cells.append(dict(key,lattice_status=cell.get("status",a.get("status","UNAVAILABLE")),
                          **{k:v for k,v in cell.items() if k!="status"}))
        tensors.append(dict(key,cell=cell))
        if a.get("physical_point_id") is not None:
            groups[str(a["physical_point_id"])].append(r)
    combined=[]
    for name,rs in groups.items():
        if len(rs)<2:continue
        # A physical ID is an explicit user assertion of same specimen/domain/point.
        xy={(r.get("x_mm"),r.get("y_mm")) for r in rs}
        if len(xy)>1:
            cell=dict(status="PHYSICAL_POINT_XY_MISMATCH")
        else:
            cell=ls.solve_from_peaks([p for r in rs for p in r["lattice_analysis"]["peaks"]],cfg)
        entry=dict(physical_point_id=name,groups=[r["point"] for r in rs],cell=cell)
        combined.append(entry)
        tensors.append(entry)
    export_csv(out/"Reflection_dSpacing.csv",rows)
    export_csv(out/"Batch_Lattice_Summary.csv",cells)
    export_csv(out/"Lattice_Parameters.csv",cells)
    export_csv(out/"Lattice_Strain.csv",[{k:v for k,v in c.items() if k in ("point","sample","x_mm","y_mm","lattice_status")
        or "strain" in k or "volume_change" in k} for c in cells])
    (out/"Strain_Tensor.json").write_text(json.dumps(tensors,indent=2),encoding="utf-8")
    (out/"Combined_Physical_Points.json").write_text(json.dumps(combined,indent=2),encoding="utf-8")
    return combined

def batch_pages(pdf,out,records,cfg,combined):
    fig=new_page("Batch lattice analysis | quality inventory","Independent acquisition groups; no automatic mixing across samples/reflections.")
    counts=collections.Counter(r.get("lattice_analysis",{}).get("status","UNAVAILABLE") for r in records)
    fig.text(.06,.83,"\n".join(f"{k.replace('_',' ')}: {v}" for k,v in counts.items()),fontsize=12,va="top",linespacing=1.6)
    fig.text(.06,.53,"Full cells: "+str(sum(r.get("lattice_analysis",{}).get("lattice",{}).get("status")=="CONDITIONAL_INDEXED_CELL" for r in records)),fontsize=12)
    fig.text(.06,.46,"Explicit cross-reflection physical point groups: "+str(len(combined)),fontsize=12)
    foot(fig,[ls.NOTICE,"All results and unavailable statuses are exported in the CSV and JSON tables."]);finish(pdf,fig)
    for start in range(0,len(records),18):
        fig=new_page("Batch lattice parameters",f"Groups {start+1}-{min(start+18,len(records))}")
        rows=[]
        for r in records[start:start+18]:
            a=r.get("lattice_analysis",{});c=a.get("lattice",{})
            pp=[p for p in a.get("peaks",[]) if "d_A" in p]
            rows.append([r["point"],str(r.get("reflection_hint",""))[:12],len(pp),
                fmt(np.mean([p["d_A"] for p in pp]) if pp else None),
                *[fmt(c.get(k),4) for k in ("a_A","b_A","c_A","strain_a_pct","strain_b_pct","strain_c_pct")],
                textwrap.fill(("HYPOTHETICAL CELL" if c.get("hypothetical") else c.get("status",a.get("status","N/A"))).replace("_"," "),22)])
        table(fig,["Group","Family","Valid d","Mean d (A)","a (A)","b (A)","c (A)","ea (%)","eb (%)","ec (%)","Cell QC"],
            rows,rect=(.04,.12,.92,.70),widths=[.055,.075,.06,.09,.075,.075,.075,.075,.075,.075,.235],fontsize=7.2)
        foot(fig,["Mean d is descriptive within one family only; it is not a universal lattice constant."]);finish(pdf,fig)
    all_poles=[(r,p) for r in records for p in r.get("peaks",[])]
    for start in range(0,len(all_poles),21):
        fig=new_page("All fitted pole metrics",f"Poles {start+1}-{min(start+21,len(all_poles))}")
        rows=[]
        for r,p in all_poles[start:start+21]:
            radial=next((q for q in r.get("lattice_analysis",{}).get("peaks",[]) if q["rank"]==p["rank"]),{})
            rows.append([r["point"],p["rank"],fmt(p.get("beta_deg"),3),fmt(p.get("tilt_deg"),3),
                fmt(p.get("phi_fwhm_observed_deg"),3),fmt(p.get("tilt_fwhm_observed_deg"),3),
                fmt(p.get("stereo_radius"),4),fmt(p.get("fit_r2"),3),fmt(radial.get("d_A"),5),
                textwrap.fill(p.get("coverage_qc","N/A").replace("_"," "),15)])
        table(fig,["Group","Pole","Beta deg","Tilt deg","FWHM beta","FWHM tilt","Stereo r","Fit R2","d (A)","Coverage QC"],
            rows,rect=(.04,.12,.92,.70),widths=[.065,.06,.09,.09,.105,.105,.10,.075,.10,.21],fontsize=7.4)
        foot(fig,["Observed pole widths are not instrument-deconvolved. Intensity rank is not crystallographic identity across points."])
        finish(pdf,fig)
    by=collections.defaultdict(list)
    for r in records:
        if r.get("status") in ("OK","PARTIAL"):
            by[(r.get("source_folder"),r["sample"],r.get("reflection_hint"),r.get("two_theta_deg"))].append(r)
    for key,rs in by.items():
        if len(rs)<2: continue
        metrics=collections.OrderedDict()
        for k,lab in (("a_A","a (A)"),("b_A","b (A)"),("c_A","c (A)"),("strain_a_pct","a change (%)"),
            ("strain_b_pct","b change (%)"),("strain_c_pct","c change (%)"),("volume_A3","Volume (A^3)"),
            ("volume_change_pct","Volume change (%)"),("c_over_a","c/a"),("a_over_b","a/b"),
            ("longitudinal_strain_pct","Tape longitudinal (%)"),("transverse_strain_pct","Tape transverse (%)"),
            ("normal_strain_pct","Specimen normal stretch (%)")):
            vals=[r.get("lattice_analysis",{}).get("lattice",{}).get(k,np.nan) for r in rs]
            if np.any(np.isfinite(vals)):metrics[lab]=np.asarray(vals,float)
        for k,lab in (("axis_tilt_deg","Conditional cone axis tilt (deg)"),("cone_rms_deg","Cone RMS (deg)")):
            vals=[r.get("structure_tilt",{}).get(k,np.nan) for r in rs]
            if np.any(np.isfinite(vals)):metrics[lab]=np.asarray(vals,float)
        for i in range(3):
            vals=[r.get("lattice_analysis",{}).get("lattice",{}).get("principal_strain_pct",[np.nan]*3)[i] for r in rs]
            if np.any(np.isfinite(vals)):metrics[f"Principal {i+1} strain (%)"]=np.asarray(vals,float)
        # Every signed reflection (or pole rank when unassigned) gets its own
        # series. Rank varies with intensity: do not claim identity across XY.
        labels=collections.defaultdict(dict)
        for i,r in enumerate(rs):
            for p in r.get("peaks",[]):
                if p.get("status")!="OK":continue
                for k,unit in (("beta_deg","beta center (deg)"),("tilt_deg","tilt center (deg)"),
                    ("phi_fwhm_observed_deg","observed beta FWHM (deg)"),("tilt_fwhm_observed_deg","observed tilt FWHM (deg)"),
                    ("stereo_radius","stereographic radius"),("fit_r2","fit R2")):
                    if k in p:labels[f"rank {p['rank']} (not tracked) | {unit}"][i]=p[k]
            for p in r.get("lattice_analysis",{}).get("peaks",[]):
                h=p.get("indexed_hkl")
                label=("hkl "+str(h)) if h and p.get("indexing_status")=="USER_ASSIGNED" else ("aggregate (not a pole)" if p["rank"]==0 else "rank "+str(p["rank"])+" (not tracked)")
                for k,unit in (("d_A","d (A)"),("normal_strain_pct","d change (%)"),("radial_fwhm_deg","Radial FWHM (deg)"),("d_fit_sigma_A","Fit sigma d (A)")):
                    if p.get("radial_status")=="OK" and k in p:
                        labels[label+" | "+unit][i]=p[k]
        for label,values in labels.items():
            v=np.full(len(rs),np.nan)
            for i,val in values.items():v[i]=val
            metrics[label]=v
        if not metrics:continue
        # Conditional directional strain vs texture correlations; never combine unlike series.
        if len(rs)>=3:
            fig=new_page("Within-series strain and texture correlations",str(key[1])+" | "+str(key[2]))
            pairs=[]
            for label,vals in metrics.items():
                if "(%)" not in label or "Fit sigma" in label:continue
                for metric,arr in metrics.items():
                    if not any(k in metric for k in ("cone axis tilt","observed beta FWHM","observed tilt FWHM")):continue
                    good=np.isfinite(vals)&np.isfinite(arr)
                    if good.sum()>=3 and np.std(vals[good])>1e-10 and np.std(arr[good])>1e-10:
                        pairs.append((label,metric,int(good.sum()),float(np.corrcoef(vals[good],arr[good])[0,1])))
            if pairs:
                for start in range(0,len(pairs),18):
                    f=fig if start==0 else new_page("Within-series strain and texture correlations",str(key[1]))
                    table(f,["Strain metric","Texture metric","n","Pearson r"],
                        [[textwrap.fill(a,38),textwrap.fill(b,38),n,fmt(c,3)] for a,b,n,c in pairs[start:start+18]],
                        rect=(.05,.12,.90,.69),widths=[.42,.42,.07,.09],fontsize=7.5)
                    foot(f,["Descriptive correlation only. Hypothetical cells and untracked intensity ranks must be interpreted with caution."])
                    finish(pdf,f)
            else:plt.close(fig)
        xy=np.array([[r.get("x_mm") if r.get("x_mm") is not None else np.nan,r.get("y_mm") if r.get("y_mm") is not None else np.nan] for r in rs],float)
        for start in range(0,len(metrics),3):
            fig=new_page("Lattice trends and measured XY locations",str(key[1])+" | "+str(key[2]))
            axs=fig.subplots(3,2)
            for row,(lab,values) in enumerate(list(metrics.items())[start:start+3]):
                ax,bx=axs[row]
                ax.plot([r["point"] for r in rs],values,"o-",ms=3,lw=1)
                ax.set_title(lab,fontsize=8.5);ax.set_xlabel("Acquisition group",fontsize=7);ax.grid(alpha=.2)
                ok=np.isfinite(values)&np.isfinite(xy).all(axis=1)
                if ok.any():
                    signed="(%)" in lab
                    lim=max(float(np.nanmax(np.abs(values))),1e-9)
                    kw=dict(cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim)) if signed else dict(cmap="viridis")
                    im=bx.scatter(xy[ok,0],xy[ok,1],c=values[ok],s=35,edgecolors="#405364",linewidths=.3,**kw)
                    fig.colorbar(im,ax=bx,pad=.04).ax.tick_params(labelsize=6)
                    bx.set_xlabel("Stage X (mm)",fontsize=7);bx.set_ylabel("Stage Y (mm)",fontsize=7)
                else:bx.text(.5,.5,"No XY coordinate",ha="center",transform=bx.transAxes)
                bx.set_title("Measured points only",fontsize=8.5)
                for axis in (ax,bx):axis.tick_params(labelsize=7)
            for row in range(min(3,len(metrics)-start),3):
                for ax in axs[row]:ax.axis("off")
            fig.subplots_adjust(left=.09,right=.92,top=.84,bottom=.10,wspace=.35,hspace=.65)
            foot(fig,["Signed maps use a zero-centered scale within each metric. No interpolation or cross-reflection merging."]);finish(pdf,fig)
        metric_items=list(metrics.items())
        for start in range(0,len(metric_items),9):
            fig=new_page("Lattice distributions and descriptive statistics",str(key[1])+" | "+str(key[2]))
            axs=fig.subplots(3,3).flat
            for j,ax in enumerate(axs):
                if start+j>=len(metric_items):ax.axis("off");continue
                lab,v=metric_items[start+j];v=v[np.isfinite(v)]
                ax.hist(v,bins=min(12,max(1,len(v))),color="#4e83a8")
                ax.set_title(textwrap.fill(lab,35),fontsize=8)
                ax.set_xlabel(f"n={len(v)} mean={np.mean(v):.5g} SD={np.std(v):.3g}\nmin={np.min(v):.5g} max={np.max(v):.5g}",fontsize=7)
                ax.tick_params(labelsize=6)
            fig.subplots_adjust(left=.08,right=.96,top=.84,bottom=.12,wspace=.36,hspace=.80)
            finish(pdf,fig)
    for r in combined:
        subtitle="Explicit physical point "+str(r["physical_point_id"])+" | groups "+str(r["groups"])
        cell_pages(pdf,r["cell"],cfg,subtitle)
        folder=out/"Combined_Physical_Points"
        folder.mkdir(exist_ok=True)
        name="Physical_Point_"+"_".join(str(x) for x in r["groups"])+".pdf"
        with PdfPages(folder/name) as individual:
            cell_pages(individual,r["cell"],cfg,subtitle)
        for group in records:
            if group["point"] in r["groups"] and (out/group["folder"]).is_dir():
                (out/group["folder"]/"Combined_Physical_Point.json").write_text(
                    json.dumps(dict(r,report="../Combined_Physical_Points/"+name),indent=2),encoding="utf-8")
