"""Lattice/strain analysis from measured radial diffraction profiles.

Reciprocal basis convention: g=(s_out-s_in)/lambda, |g|=1/d, no 2*pi.
No phase identification, instrument calibration, or stress inference is implied.
"""
from __future__ import annotations
import itertools
import json
import math
import re
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

DEFAULT = dict(material="LMO (nominal pseudocubic reference)", wavelength_A=1.5406,
    wavelength_confirmed=False, cell_A_deg=[3.880,3.880,3.880,90.,90.,90.],
    reference_confirmed=False, calibration_confirmed=False,
    two_theta_offset_deg=0., two_theta_systematic_sigma_deg=None,
    tape_azimuth_deg=None, cube_texture_azimuth_deg=None,
    assume_single_domain=False, angular_sigma_deg=None,
    indexing=[], detector_center_offset_px=[0.,0.], detector_distance_offset_mm=0.,
    secondary_wavelength_A=None, secondary_intensity_ratio=.5,
    allow_hypothetical_four_pole_cell=True, use_header_wavelength=True)
RADIAL_STEP = .025
RADIAL_HALF = 2.4
NOTICE = ("Conditional diffraction estimates; phase, wavelength, reference and instrument "
          "calibration require confirmation. Observed radial widths are not resolution corrected.")

def hkl_parse(value):
    if isinstance(value,(list,tuple)):
        if len(value)!=3 or any(int(x)!=x for x in value): raise ValueError("hkl must contain three integers")
        h=tuple(int(x) for x in value)
    else:
        s=str(value).strip().strip("(){}[]").strip()
        if re.fullmatch(r"\d{3}",s): h=tuple(map(int,s))
        elif re.fullmatch(r"[+-]?\d+[\s,]+[+-]?\d+[\s,]+[+-]?\d+",s):
            h=tuple(map(int,re.split(r"[\s,]+",s)))
        else: return None
    return h if any(h) else None

def hint_from_name(s):
    # Explicit LMO/PF/hkl or a leading hkl followed by _2D. Do not use sample IDs.
    for pattern in (r"(?i)(?:^|[^A-Za-z0-9])(?:LMO|PF|hkl)[_ -]*\(?(\d{3}|[+-]?\d+[, ]+[+-]?\d+[, ]+[+-]?\d+)\)?(?=$|[^0-9])",
                    r"(?i)(?:^|[^A-Za-z0-9])(\d{3})_2D(?:_|$)"):
        m=re.search(pattern,s)
        if m and hkl_parse(m[1]): return m[1]
    return "not specified"

def basis(cell):
    a,b,c,al,be,ga=map(float,cell)
    if not all(np.isfinite([a,b,c,al,be,ga])) or min(a,b,c)<=0 or not all(0<x<180 for x in (al,be,ga)):
        raise ValueError("Invalid reference cell")
    ca,cb,cg=np.cos(np.deg2rad([al,be,ga])); sg=np.sin(np.deg2rad(ga))
    v=(ca-cb*cg)/sg; z=1-cb*cb-v*v
    if z<=1e-10: raise ValueError("Reference cell is not positive definite")
    return np.array([[a,b*cg,c*cb],[0,b*sg,c*v],[0,0,c*np.sqrt(z)]])

def settings(path=None):
    cfg=json.loads(json.dumps(DEFAULT))
    if path:
        incoming=json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if not isinstance(incoming,dict): raise ValueError("Reference JSON must be an object")
        unknown=set(incoming)-set(cfg)
        if unknown: raise ValueError("Unknown reference settings: "+", ".join(sorted(unknown)))
        cfg.update(incoming)
    basis(cfg["cell_A_deg"])
    for k in ("wavelength_A","two_theta_offset_deg","detector_distance_offset_mm"):
        cfg[k]=float(cfg[k])
        if not np.isfinite(cfg[k]): raise ValueError("Nonfinite "+k)
    if cfg["wavelength_A"]<=0: raise ValueError("Positive wavelength required")
    for k in ("tape_azimuth_deg","cube_texture_azimuth_deg","two_theta_systematic_sigma_deg","angular_sigma_deg"):
        if cfg[k] is not None:
            cfg[k]=float(cfg[k])
            if not np.isfinite(cfg[k]) or ("sigma" in k and cfg[k]<=0): raise ValueError("Invalid "+k)
    if cfg["secondary_wavelength_A"] is not None:
        cfg["secondary_wavelength_A"]=float(cfg["secondary_wavelength_A"])
        if not np.isfinite(cfg["secondary_wavelength_A"]) or cfg["secondary_wavelength_A"]<=cfg["wavelength_A"]:
            raise ValueError("Secondary wavelength must be greater than the primary wavelength")
    cfg["secondary_intensity_ratio"]=float(cfg["secondary_intensity_ratio"])
    if not np.isfinite(cfg["secondary_intensity_ratio"]) or not 0<cfg["secondary_intensity_ratio"]<=1:
        raise ValueError("Secondary intensity ratio must be in (0,1]")
    for k in ("wavelength_confirmed","reference_confirmed","calibration_confirmed","assume_single_domain","allow_hypothetical_four_pole_cell","use_header_wavelength"):
        if type(cfg[k]) is not bool: raise ValueError(k+" must be true or false")
    shift=np.asarray(cfg["detector_center_offset_px"],float)
    if shift.shape!=(2,) or not np.isfinite(shift).all(): raise ValueError("Detector center offset must have two finite values")
    cfg["detector_center_offset_px"]=shift.tolist()
    if not isinstance(cfg["indexing"],list): raise ValueError("indexing must be a list")
    for rule in cfg["indexing"]:
        if not isinstance(rule,dict) or not isinstance(rule.get("sample"),str) or not rule["sample"]:
            raise ValueError("Each indexing rule requires an exact sample string")
        if "family_hkl" in rule and hkl_parse(rule["family_hkl"]) is None: raise ValueError("Invalid family_hkl")
        for rank,h in rule.get("hkl_by_rank",{}).items():
            if not str(rank).isdigit() or int(rank)<1 or hkl_parse(h) is None:
                raise ValueError("Invalid signed hkl_by_rank entry")
    return cfg

def effective_settings(cfg,header):
    """Explicit confirmed/disabled manual wavelength wins; otherwise use recorded line.
    Generic Ka metadata lists both wavelengths but does not establish actual optics.
    """
    out=dict(cfg)
    wave=(header.get("wave_type") or "").lower().replace("_","").replace("-","").replace(" ","")
    lines=header.get("source_wavelengths_A") or {}
    if cfg.get("use_header_wavelength",True) and not cfg.get("wavelength_confirmed"):
        key="Kbeta" if wave in ("kb","kbeta","kb1") else ("Kalpha2" if wave in ("ka2","kalpha2") else "Kalpha1" if wave in ("ka","kalpha","ka1","kalpha1") else None)
        if key and lines.get(key) is not None and lines[key]>0:
            out["wavelength_A"]=lines[key]
    return out

def write_settings(dest,cfg):
    (Path(dest)/"Lattice_Settings.json").write_text(json.dumps(cfg,indent=2),encoding="utf-8")

def reference_d(hkl,cfg):
    if hkl is None: return None
    g=np.linalg.solve(basis(cfg["cell_A_deg"]).T,np.asarray(hkl,float))
    return float(1/np.linalg.norm(g))

def radial_axis(center):
    n=round(2*RADIAL_HALF/RADIAL_STEP)
    return center-RADIAL_HALF+(np.arange(n)+.5)*RADIAL_STEP

class RadialCube:
    """Rate sum / pixel coverage, accumulated during the existing raw image pass."""
    def __init__(self,center,nt=90,na=180):
        self.axis=radial_axis(center)
        self.shape=(nt,na,len(self.axis))
        self.sums=np.zeros(self.shape,np.float32)
        self.weights=np.zeros_like(self.sums)
        self.lower=center-RADIAL_HALF
    def bins(self,tth):
        return np.floor((tth-self.lower)/RADIAL_STEP).astype(np.intp)
    def add(self,pole_bins,radial_bins,rate,weight):
        valid=(radial_bins>=0)&(radial_bins<self.shape[2])&(pole_bins>=0)&(pole_bins<self.shape[0]*self.shape[1])
        idx=pole_bins[valid]*self.shape[2]+radial_bins[valid]
        np.add.at(self.sums.ravel(),idx,rate[valid]*weight)
        np.add.at(self.weights.ravel(),idx,weight)
    def save(self,path):
        np.savez_compressed(path, two_theta_deg=self.axis, sums=self.sums, coverage=self.weights,
            tilt_deg=np.arange(self.shape[0])+.5,
            beta_deg=(np.arange(self.shape[1])+.5)*360/self.shape[1])

def direction(tilt,beta):
    t,b=np.deg2rad([tilt,beta])
    return np.array([np.sin(t)*np.cos(b),np.sin(t)*np.sin(b),np.cos(t)])

def fit_radial(x,y,coverage,cfg=None):
    """Observed peak center, free linear background; reject edge/unresolved/poor fits."""
    cfg=cfg or DEFAULT
    secondary=cfg.get("secondary_wavelength_A")
    correction=cfg.get("two_theta_offset_deg",0.)
    ratio=cfg.get("secondary_intensity_ratio",.5)
    def second_center(mu):
        arg=secondary/cfg["wavelength_A"]*math.sin(math.radians((mu+correction)/2))
        return math.degrees(2*math.asin(min(arg,1.)))-correction
    mask=np.isfinite(y)&(coverage>0)
    result=dict(radial_status="INSUFFICIENT_RADIAL_COVERAGE",n_radial_bins=int(mask.sum()))
    if mask.sum()<12: return result,None
    xx=x[mask]; yy=y[mask]; x0=float(np.mean(xx))
    step=float(np.median(np.diff(x)))
    bg=float(np.percentile(yy,20)); amp=float(np.max(yy)-bg)
    if not np.isfinite(amp) or amp<=0: return result,None
    seed=float(xx[np.argmax(yy)])
    def prediction(p,x):
        a,mu,s,b,m=p
        profile=b+m*(x-x0)+a*np.exp(-.5*((x-mu)/s)**2)
        if secondary is not None:
            profile+=a*ratio*np.exp(-.5*((x-second_center(mu))/s)**2)
        return profile
    # Account for unequal counts of sampled detector pixels without claiming Poisson errors.
    w=np.sqrt(coverage[mask]/max(float(np.max(coverage[mask])),1.))
    try:
        fit=least_squares(lambda p:(prediction(p,xx)-yy)*w,
             [amp,seed,max(.12,step),bg,0.],
             bounds=([0.,xx[0],step/3,-np.inf,-np.inf],
                     [np.inf,xx[-1],.65,np.inf,np.inf]),max_nfev=500)
        a,mu,s,b,m=map(float,fit.x); width=2.354820045*s
        resid=prediction(fit.x,xx)-yy
        denom=np.sum((yy-yy.mean())**2)
        r2=float(1-np.sum(resid**2)/denom) if denom>0 else 0.
        sigma=None
        if np.linalg.matrix_rank(fit.jac)==5 and len(xx)>5:
            cov=np.linalg.pinv(fit.jac.T@fit.jac)*(2*fit.cost/(len(xx)-5))
            sigma=float(np.sqrt(max(0.,cov[1,1])))
        status="OK"
        if not fit.success or r2<.80 or sigma is None: status="POOR_RADIAL_FIT"
        elif width<2*step: status="UNRESOLVED_RADIAL_PEAK"
        elif mu-1.3*width<=xx[0] or (second_center(mu) if secondary is not None else mu)+1.3*width>=xx[-1]: status="CLIPPED_RADIAL_PEAK"
        else:
            within=(x>=mu-width)&(x<=(second_center(mu) if secondary is not None else mu)+width)
            if not np.all(mask[within]): status="GAPPED_RADIAL_PEAK"
        result.update(radial_status=status, observed_two_theta_deg=mu,
            radial_fwhm_deg=width, radial_fit_r2=r2, two_theta_fit_sigma_deg=sigma,
            radial_coverage_min_deg=float(xx[0]),radial_coverage_max_deg=float(xx[-1]),
            radial_model="fixed-ratio equal-width wavelength doublet" if secondary is not None else "single Gaussian",
            secondary_two_theta_deg=second_center(mu) if secondary is not None else None)
        return result,prediction(fit.x,x)
    except (ValueError,np.linalg.LinAlgError,FloatingPointError) as exc:
        return dict(result,radial_status="RADIAL_FIT_FAILED",radial_error=str(exc)),None

def matching_rule(item,cfg):
    hits=[r for r in cfg["indexing"] if r.get("sample")==item["sample"] and
          ("position_suffix" not in r or str(r["position_suffix"])==str(item.get("position_suffix"))) and
          ("source_folder" not in r or r["source_folder"]==item.get("source_folder"))]
    if len(hits)>1: raise ValueError("Multiple indexing rules match "+item["sample"])
    return hits[0] if hits else {}

def analyze_profiles(dest,item,peaks,cfg):
    path=Path(dest)/"Radial_Cube.npz"
    family_hint=hkl_parse(item.get("reflection_hint",""))
    result=dict(status="RAW_RADIAL_DATA_REQUIRED",peaks=[],notice=NOTICE,
                family_hkl=list(family_hint) if family_hint else None)
    if not path.is_file(): return result
    with np.load(path,allow_pickle=False) as z:
        x=z["two_theta_deg"]; sums=z["sums"]; cov=z["coverage"]
        tilt=z["tilt_deg"]; beta=z["beta_deg"]
    rule=matching_rule(item,cfg)
    family=hkl_parse(rule.get("family_hkl",item.get("reflection_hint","")))
    result["family_hkl"]=list(family) if family else None
    result["indexing_basis"]="explicit rule" if rule else "filename hint only"
    result["physical_point_id"]=rule.get("physical_point_id")
    profiles=[]; fits=[]; coverages=[]; ranks=[]
    for p in peaks:
        if p.get("status")!="OK" or "sigma_tilt_deg" not in p: continue
        # Central angular aperture, consistent with the old measured pole map.
        ta=np.abs(tilt-p["tilt_deg"]); ba=np.abs((beta-p["beta_deg"]+180)%360-180)
        aperture=(ta[:,None]<=max(1.5,p["sigma_tilt_deg"]))&(ba[None,:]<=max(3.,p["sigma_beta_deg"]))
        cc=cov[aperture].sum(axis=0,dtype=float)
        ss=sums[aperture].sum(axis=0,dtype=float)
        y=np.divide(ss,cc,out=np.full_like(ss,np.nan),where=cc>0)
        r,model=fit_radial(x,y,cc,cfg)
        r.update(rank=int(p["rank"]),beta_deg=float(p["beta_deg"]),tilt_deg=float(p["tilt_deg"]),
                 pole_coverage_qc=p.get("coverage_qc","UNKNOWN"))
        # Historical azimuth map bins have centers offset +half-bin from their
        # actual rounding centers. Preserve old pole reports; correct only g.
        r["vector_beta_deg"]=(r["beta_deg"]-180/len(beta))%360
        if r["radial_status"]=="OK":
            tt=r["observed_two_theta_deg"]+cfg["two_theta_offset_deg"]
            if 0<tt<180:
                d=cfg["wavelength_A"]/(2*math.sin(math.radians(tt/2)))
                r.update(corrected_two_theta_deg=tt,bragg_theta_deg=tt/2,d_A=d,wavelength_A=cfg["wavelength_A"],
                         g_inv_A=(direction(r["tilt_deg"],r["vector_beta_deg"])/d).tolist())
                sig=r.get("two_theta_fit_sigma_deg")
                if sig is not None:
                    r["d_fit_sigma_A"]=abs(d*math.radians(sig)/2/math.tan(math.radians(tt/2)))
                sy=cfg["two_theta_systematic_sigma_deg"]
                if sy is not None:
                    r["d_angle_sigma_A"]=abs(d*math.radians(math.hypot(sig or 0,sy))/2/math.tan(math.radians(tt/2)))
            else: r["radial_status"]="INVALID_CORRECTED_ANGLE"
        ranks.append(p["rank"]); profiles.append(y); coverages.append(cc)
        fits.append(model if model is not None else np.full_like(x,np.nan))
        result["peaks"].append(r)
    if not any("d_A" in p for p in result["peaks"]):
        cc=cov.sum(axis=(0,1),dtype=float);ss=sums.sum(axis=(0,1),dtype=float)
        y=np.divide(ss,cc,out=np.full_like(ss,np.nan),where=cc>0)
        r,model=fit_radial(x,y,cc,cfg)
        r.update(rank=0,pole_coverage_qc="AGGREGATE_NOT_A_POLE",indexing_status="AGGREGATE_NOT_A_POLE")
        if r["radial_status"]=="OK":
            tt=r["observed_two_theta_deg"]+cfg["two_theta_offset_deg"]
            if 0<tt<180:
                d=cfg["wavelength_A"]/(2*math.sin(math.radians(tt/2)))
                r.update(corrected_two_theta_deg=tt,bragg_theta_deg=tt/2,d_A=d,wavelength_A=cfg["wavelength_A"])
                sig=r.get("two_theta_fit_sigma_deg")
                if sig is not None:r["d_fit_sigma_A"]=abs(d*math.radians(sig)/2/math.tan(math.radians(tt/2)))
            else:r["radial_status"]="INVALID_CORRECTED_ANGLE"
        result["peaks"].append(r);ranks.append(0);profiles.append(y);coverages.append(cc)
        fits.append(model if model is not None else np.full_like(x,np.nan))
    shape=(0,len(x))
    np.savez_compressed(Path(dest)/"Radial_Profiles.npz",two_theta_deg=x,
        rank=np.asarray(ranks,int),intensity=np.array(profiles) if profiles else np.empty(shape),
        model=np.array(fits) if fits else np.empty(shape),coverage=np.array(coverages) if coverages else np.empty(shape))
    # For generic families, only a cubic reference makes all permutations equivalent.
    A0=basis(cfg["cell_A_deg"])
    cubic=np.allclose(A0,np.eye(3)*A0[0,0],rtol=1e-8,atol=1e-8)
    mapping=rule.get("hkl_by_rank",{})
    explicit=bool(mapping)
    auto=auto_index(result["peaks"],family,cfg) if cubic and not explicit else {}
    hypothetical=False
    solve_cfg=cfg
    if cubic and not explicit and not auto and cfg["allow_hypothetical_four_pole_cell"] and cfg["cube_texture_azimuth_deg"] is None:
        az=infer_four_pole_orientation(result["peaks"],family)
        if az is not None:
            model_cfg=dict(cfg,assume_single_domain=True,cube_texture_azimuth_deg=az,_allow_hypothetical_solution=True)
            candidate=auto_index(result["peaks"],family,model_cfg)
            if len(candidate)==4:
                auto=candidate;hypothetical=True;solve_cfg=model_cfg
                result["indexing_basis"]="hypothetical cubic single-domain texture; a/b labels symmetry-ambiguous"
    for r in result["peaks"]:
        rank=str(r["rank"])
        h=hkl_parse(mapping[rank]) if rank in mapping else (None if explicit else auto.get(r["rank"]))
        if h and r["rank"]>0:
            r["indexed_hkl"]=list(h)
            r["indexing_status"]="USER_ASSIGNED" if explicit else ("HYPOTHETICAL_CUBE_TEXTURE" if hypothetical else "CONDITIONAL_CUBE_TEXTURE")
        elif r["rank"]!=0: r["indexing_status"]="UNASSIGNED"
        d0=reference_d(h or (family if cubic else None),cfg)
        if d0 and "d_A" in r:
            r.update(reference_d_A=d0,normal_strain_pct=100*(r["d_A"]/d0-1),
                     reference_two_theta_deg=(math.degrees(2*math.asin(cfg["wavelength_A"]/(2*d0)))
                        if cfg["wavelength_A"]<2*d0 else None))
            if cubic and family:
                r["apparent_cubic_a_A"]=float(r["d_A"]*np.linalg.norm(h or family))
    result["status"]="CONDITIONAL" if any("d_A" in r for r in result["peaks"]) else "NO_RELIABLE_RADIAL_FIT"
    if all(cfg[k] for k in ("wavelength_confirmed","reference_confirmed","calibration_confirmed")) and result["status"]=="CONDITIONAL":
        result["status"]="USER_CALIBRATION_CONFIRMED"
    result["lattice"]=solve_from_peaks(result["peaks"],solve_cfg)
    if hypothetical:
        result["lattice"]["hypothetical"]=True
        result["lattice"]["warning"]="Hypothetical single-domain cube texture; confirm phase/domain, indexing and orientation. a/b labels symmetry-ambiguous."
    result["lattice"]["indexing_basis"]=result["indexing_basis"]
    return result

def infer_four_pole_orientation(peaks,family):
    if family is None:return None
    kind=sorted(map(abs,family))
    if kind not in ([0,1,1],[1,1,1],[2,2,2]):return None
    good=[p for p in peaks if "g_inv_A" in p and p.get("pole_coverage_qc")=="OK"]
    if len(good)!=4:return None
    expected=45. if kind==[0,1,1] else math.degrees(math.acos(1/math.sqrt(3)))
    if any(abs(p["tilt_deg"]-expected)>15 for p in good):return None
    beta=np.array([p["vector_beta_deg"] for p in good])
    fourth=np.mean(np.exp(4j*np.deg2rad(beta)))
    if abs(fourth)<.75:return None
    return float((np.rad2deg(np.angle(fourth))/4-(0 if kind==[0,1,1] else 45))%90)

def auto_index(peaks,family,cfg):
    """Conservative cube-texture indexing near user-supplied crystal a azimuth.
    Explicit single-domain assumption; no phase identification. Unique nearest pole.
    """
    az=cfg["cube_texture_azimuth_deg"]
    if family is None or az is None or not cfg["assume_single_domain"]: return {}
    candidates=set()
    for p in set(itertools.permutations(map(abs,family))):
        for signs in itertools.product((-1,1),repeat=3):
            h=tuple(v*s for v,s in zip(p,signs))
            if h[2]>0: candidates.add(h)
            # In-plane Friedel sign is ambiguous in the hemisphere projection.
    if not candidates: return {}
    a=math.radians(az)
    R=np.array([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1.]])
    hs=sorted(candidates); G=np.array([R@np.array(h)/np.linalg.norm(h) for h in hs])
    assigned={}
    for r in peaks:
        if "g_inv_A" not in r or r.get("pole_coverage_qc")!="OK": continue
        u=np.asarray(r["g_inv_A"]);u/=np.linalg.norm(u)
        angles=np.rad2deg(np.arccos(np.clip(G@u,-1,1)))
        order=np.argsort(angles)
        if angles[order[0]]>15 or (len(order)>1 and angles[order[1]]-angles[order[0]]<8): continue
        assigned[r["rank"]]=hs[order[0]]
    if len(set(assigned.values()))!=len(assigned): return {}
    return assigned

def solve_cell(H,G,cfg):
    """H rows signed indexed reflections, G rows measured reciprocal vectors.
    G=H B^T, A=B^-T, F=A A0^-1. Left stretch in specimen coordinates.
    """
    H=np.asarray(H,float);G=np.asarray(G,float)
    if H.ndim!=2 or H.shape!=G.shape or H.shape[1]!=3 or len(H)<3:
        return dict(status="INSUFFICIENT_INDEXED_VECTORS",n_vectors=len(H))
    rank=int(np.linalg.matrix_rank(H));cond=float(np.linalg.cond(H))
    if rank<3 or not np.isfinite(cond) or cond>30:
        return dict(status="RANK_OR_CONDITION_FAILURE",rank=rank,condition=cond if np.isfinite(cond) else None,n_vectors=len(H))
    BT=np.linalg.lstsq(H,G,rcond=None)[0]; B=BT.T
    residual=np.linalg.norm(H@BT-G,axis=1)
    rms=float(np.sqrt(np.mean(residual**2))/np.mean(np.linalg.norm(G,axis=1)))
    if np.linalg.det(B)<=0 or np.linalg.cond(B)>30:
        return dict(status="INVALID_RECIPROCAL_BASIS",rank=rank,n_vectors=len(H))
    if rms>.015: return dict(status="INCONSISTENT_INDEXED_VECTORS",rank=rank,n_vectors=len(H),relative_rms=rms)
    A=np.linalg.inv(B).T; A0=basis(cfg["cell_A_deg"]); F=A@np.linalg.inv(A0)
    U,s,Vt=np.linalg.svd(F); R=U@Vt
    if np.linalg.det(R)<=0: return dict(status="IMPROPER_ROTATION")
    stretch=(U*s)@U.T; E=stretch-np.eye(3)
    lens=np.linalg.norm(A,axis=0); l0=np.linalg.norm(A0,axis=0)
    angles=[math.degrees(math.acos(float(np.clip(A[:,j]@A[:,k]/(lens[j]*lens[k]),-1,1)))) for j,k in ((1,2),(0,2),(0,1))]
    eigen,vectors=np.linalg.eigh(E)
    out=dict(status="CONDITIONAL_INDEXED_CELL",n_vectors=len(H),rank=rank,condition=cond,
        relative_rms=rms,redundant_vectors=len(H)-3, a_A=float(lens[0]),b_A=float(lens[1]),c_A=float(lens[2]),
        alpha_deg=angles[0],beta_deg=angles[1],gamma_deg=angles[2],
        volume_A3=float(np.linalg.det(A)),volume_change_pct=float(100*(np.linalg.det(F)-1)),
        strain_a_pct=float(100*(lens[0]/l0[0]-1)),strain_b_pct=float(100*(lens[1]/l0[1]-1)),
        strain_c_pct=float(100*(lens[2]/l0[2]-1)),c_over_a=float(lens[2]/lens[0]),a_over_b=float(lens[0]/lens[1]),
        direct_basis_A=A.tolist(),reciprocal_basis_inv_A=B.tolist(),rotation=R.tolist(),
        deformation_gradient=F.tolist(),stretch_strain_specimen_pct=(100*E).tolist(),
        principal_strain_pct=(100*eigen).tolist(),principal_directions_specimen=vectors.T.tolist(),
        normal_strain_pct=float(100*E[2,2]),
        uncertainty_note="Not a calibrated uncertainty budget; angular texture and model systematics may dominate.")
    if cfg["tape_azimuth_deg"] is not None:
        a=math.radians(cfg["tape_azimuth_deg"]); v=np.array([math.cos(a),math.sin(a),0]); w=np.array([-v[1],v[0],0])
        out.update(longitudinal_strain_pct=float(100*v@E@v),transverse_strain_pct=float(100*w@E@w))
    if len(H)==3: out["warning"]="Exactly determined; no independent residual check."
    return out

def solve_from_peaks(peaks,cfg):
    valid=[p for p in peaks if p.get("indexed_hkl") and "g_inv_A" in p and p.get("pole_coverage_qc")=="OK"
           and (p.get("indexing_status")!="HYPOTHETICAL_CUBE_TEXTURE" or cfg.get("_allow_hypothetical_solution",False))]
    if not cfg["assume_single_domain"]:
        return dict(status="SINGLE_DOMAIN_NOT_CONFIRMED",n_vectors=len(valid))
    if len(valid)<3: return dict(status="INSUFFICIENT_INDEXED_VECTORS",n_vectors=len(valid))
    H=np.array([p["indexed_hkl"] for p in valid]); G=np.array([p["g_inv_A"] for p in valid])
    out=solve_cell(H,G,cfg)
    # Propagate only specified angular and 2theta uncertainties. No invented precision.
    if out["status"]=="CONDITIONAL_INDEXED_CELL" and cfg["angular_sigma_deg"] is not None and cfg["two_theta_systematic_sigma_deg"] is not None:
        rng=np.random.default_rng(153); trials=[]
        for _ in range(150):
            offset=rng.normal(0,cfg["two_theta_systematic_sigma_deg"])
            gg=[]
            for p in valid:
                tt=p["corrected_two_theta_deg"]+offset+rng.normal(0,p.get("two_theta_fit_sigma_deg") or 0)
                t=p["tilt_deg"]+rng.normal(0,cfg["angular_sigma_deg"])
                b=p["vector_beta_deg"]+rng.normal(0,cfg["angular_sigma_deg"])
                gg.append(direction(t,b)*2*math.sin(math.radians(tt/2))/p.get("wavelength_A",cfg["wavelength_A"]))
            z=solve_cell(H,gg,cfg)
            if z["status"]=="CONDITIONAL_INDEXED_CELL": trials.append([z[k] for k in ("a_A","b_A","c_A","strain_a_pct","strain_b_pct","strain_c_pct")])
        if len(trials)>=100:
            out["model_sigma"]={k:float(v) for k,v in zip(("a_A","b_A","c_A","strain_a_pct","strain_b_pct","strain_c_pct"),np.std(trials,axis=0,ddof=1))}
            out["uncertainty_note"]="Monte Carlo of supplied angle sigmas only; excludes phase, reference, wavelength and domain uncertainty."
        out["monte_carlo_accepted"]=len(trials)
    return out
