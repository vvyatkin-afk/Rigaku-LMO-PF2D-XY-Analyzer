"""Reconstruct measured LMO(011) pole intensities and phi profiles, read-only originals.
Geometry follows Rigaku d*TREK Cspatial/Cdetector/Cgoniometer source; no symmetry imposed.
"""
import pathlib,numpy as np,json,csv,sys
from scipy.optimize import least_squares,brentq
from scipy.special import erf
P=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(__file__).parent
D=np.load(P/"raw_frames.npz")
A=D["counts"].reshape(2,180,385,775)
PHI=D["phi"].reshape(2,180)
assert np.array_equal(PHI,np.tile(np.arange(0,360,2),(2,1)))
assert np.array_equal(D["chi"],np.repeat([38.,42.],180))
V,U=np.mgrid[:385,:775]
du=(U-382.5)*.1;dv=(V-197.5)*.1; th=np.deg2rad(32.6);om=np.deg2rad(16.3)
r=np.stack((-du,-dv*np.cos(th)-150.608*np.sin(th),dv*np.sin(th)-150.608*np.cos(th)),axis=-1)
sf=r/np.linalg.norm(r,axis=-1,keepdims=True)
TT=np.rad2deg(np.arccos(-sf[...,2]))
q=sf+[0,0,1]; q/=np.linalg.norm(q,axis=-1,keepdims=True)
q0=np.stack((q[...,0],np.cos(om)*q[...,1]-np.sin(om)*q[...,2],np.sin(om)*q[...,1]+np.cos(om)*q[...,2]),axis=-1)
GEOM=[]
for chi in [38,42]:
 al=np.deg2rad(90-chi)
 qs=np.stack((q0[...,0]*np.cos(al)+q0[...,1]*np.sin(al),-q0[...,0]*np.sin(al)+q0[...,1]*np.cos(al),q0[...,2]),axis=-1)
 qs=np.where((qs[...,0]>0)[...,None],-qs,qs)
 GEOM.append((np.rad2deg(np.arccos(-qs[...,0])),np.rad2deg(np.arctan2(qs[...,2],-qs[...,1]))%360))
BANDS={"wide":[(31.6,33.3),(30.8,31.3),(33.7,34.2)],"narrow":[(31.9,33.0),(30.8,31.3),(33.7,34.2)]}
PROFILES={}; TABLE=[]; MAP={}
for label,bands in BANDS.items():
 for lo,hi in [(30,60),(35,55),(38,53),(40,50),(20,65)]:
  means=[];variances=[]
  for g,chi in enumerate([38,42]):
   tilt,az=GEOM[g]; sel=(tilt>=lo)&(tilt<hi)
   masks=[sel&(TT>=a)&(TT<b) for a,b in bands]
   nn=np.array([m.sum() for m in masks])
   cen=np.array([TT[m].mean() for m in masks])
   f=(cen[0]-cen[1])/(cen[2]-cen[1]);coeff=np.array([1.,-(1-f),-f])
   cnt=np.array([A[g][:,m].sum(axis=1) for m in masks],dtype=float)
   net=(cnt/nn[:,None]*coeff[:,None]).sum(axis=0)
   var=(cnt/nn[:,None]**2*coeff[:,None]**2).sum(axis=0)
   key=f"{label}_t{lo}_{hi}_chi{chi}"
   PROFILES[key]=(net,var)
   means.append(net);variances.append(var)
   TABLE.append(dict(key=key,pixels=nn.tolist(),mean_two_theta=cen.tolist(),background_weight_right=float(f)))
  PROFILES[f"{label}_t{lo}_{hi}_combined"]=(np.mean(means,axis=0),np.sum(variances,axis=0)/4)
 if label!="wide": continue
 for g,chi in enumerate([38,42]):
  tilt,az=GEOM[g]; out=[]
  for ba,bb in bands:
   m=(TT>=ba)&(TT<bb)&(tilt<90)
   ti=np.floor(tilt[m]).astype(int)
   base=az[m];vals=A[g][:,m]
   sums=np.zeros((90,180));num=np.zeros_like(sums);ttsum=np.zeros_like(sums)
   for j,phi in enumerate(PHI[g]):
    bi=np.floor(((base+phi+1)%360)/2).astype(int); ix=ti*180+bi
    sums+=np.bincount(ix,weights=vals[j],minlength=90*180).reshape(90,180)
    num+=np.bincount(ix,minlength=90*180).reshape(90,180)
    ttsum+=np.bincount(ix,weights=TT[m],minlength=90*180).reshape(90,180)
   out.append((sums,num,ttsum))
  sums=np.array([x[0] for x in out]);num=np.array([x[1] for x in out]);ttsum=np.array([x[2] for x in out])
  valid=(num>0).all(axis=0)
  with np.errstate(divide="ignore",invalid="ignore"):
   cen=ttsum/num;fr=(cen[0]-cen[1])/(cen[2]-cen[1])
   net=sums[0]/num[0]-(1-fr)*sums[1]/num[1]-fr*sums[2]/num[2]
   var=sums[0]/num[0]**2+(1-fr)**2*sums[1]/num[1]**2+fr**2*sums[2]/num[2]**2
  net[~valid]=np.nan;var[~valid]=np.nan
  MAP[f"chi{chi}_net"]=net;MAP[f"chi{chi}_variance"]=var;MAP[f"chi{chi}_coverage"]=num[0]
  MAP[f"chi{chi}_signal_counts"]=sums[0];MAP[f"chi{chi}_background_left_counts"]=sums[1];MAP[f"chi{chi}_background_right_counts"]=sums[2]
# Equal exposure and near-identical ray sampling in overlap; use pixel-coverage weights.
nets=np.array([MAP[f"chi{c}_net"] for c in [38,42]])
vars=np.array([MAP[f"chi{c}_variance"] for c in [38,42]])
cov=np.array([MAP[f"chi{c}_coverage"] for c in [38,42]])
cov[~np.isfinite(nets)]=0
den=cov.sum(axis=0);weights=np.divide(cov,den,out=np.zeros_like(cov),where=den>0)
MAP["combined_net"]=np.nansum(nets*weights,axis=0);MAP["combined_variance"]=np.nansum(vars*weights**2,axis=0)
MAP["combined_net"][den==0]=np.nan;MAP["combined_variance"][den==0]=np.nan
MAP["combined_coverage"]=den
MAP["tilt_centers"]=np.arange(90)+.5;MAP["azimuth_centers"]=np.arange(0,360,2)
np.savez_compressed(P/"pole_intensities.npz",**MAP)
# Circular local peak fit: Gaussian averaged through the actual continuous 2 degree phi sweep.
def dist(x): return (x+180)%360-180
def peak(x,mu,width):
 s=width/(2*np.sqrt(2*np.log(2))); h=1.
 return s*np.sqrt(np.pi/2)/(2*h)*(erf((x-mu+h)/(np.sqrt(2)*s))-erf((x-mu-h)/(np.sqrt(2)*s)))
def observed_width(w):
 half=float(peak(0,0,w))/2
 return 2*brentq(lambda x:float(peak(x,0,w))-half,0,3*w+2)
def fitone(y,var,center,window=22):
 xx=dist(np.arange(0,360,2)-center);m=np.abs(xx)<=window
 x=xx[m]; yy=y[m];sig=np.sqrt(var[m]); baseline=np.quantile(yy,.1)
 p0=[float(max(yy)-baseline),0.,10.,float(baseline),0.]
 def model(p,xx): return p[0]*peak(xx,p[1],p[2])+p[3]+p[4]*xx
 r=least_squares(lambda p:(model(p,x)-yy)/sig,p0,bounds=([0,-10,2,-.3,-.02],[5,10,35,.3,.02]),x_scale="jac")
 n=len(x);dof=n-5;red=float(np.dot(r.fun,r.fun)/dof)
 covar=np.linalg.pinv(r.jac.T@r.jac)
 se=np.sqrt(np.diag(covar)); se_scaled=se*np.sqrt(max(1,red))
 wobs=observed_width(r.x[2]);der=(observed_width(r.x[2]+.001)-observed_width(r.x[2]-.001))/.002
 return dict(center_deg=float((center+r.x[1])%360),center_se_deg=float(se[1]),center_se_scaled_deg=float(se_scaled[1]),fwhm_observed_deg=float(wobs),fwhm_se_deg=float(se[2]*der),fwhm_se_scaled_deg=float(se_scaled[2]*der),gaussian_component_fwhm_deg=float(r.x[2]),amplitude=float(r.x[0]),background=float(r.x[3]),background_slope_per_degree=float(r.x[4]),reduced_chi_square=red,n_points=n,fit_half_window_deg=window,params=r.x.tolist(),seed_center_deg=center)
FITS={}
for key,(yy,var) in PROFILES.items():
 FITS[key]=[fitone(yy,var,c) for c in [82,177,272,357]]
# Window sensitivity for primary combined profile.
for w in [18,26]:
 yy,var=PROFILES["wide_t30_60_combined"]
 FITS[f"wide_t30_60_combined_window{w}"]=[fitone(yy,var,c,w) for c in [82,177,272,357]]
np.savez_compressed(P/"phi_profiles.npz",phi=np.arange(0,360,2),**{k:v[0] for k,v in PROFILES.items()},**{k+"_variance":v[1] for k,v in PROFILES.items()})
OUT=dict(method="Independent circular peak fits, Gaussian convolved with measured 2 degree continuous scan window, linear local background; Poisson propagated count uncertainties; observed FWHM reported, no instrumental deconvolution",signal_bands=BANDS,profile_settings=TABLE,fits=FITS,source_count=360,frames_per_chi=180,exposure_seconds=5,raw_shape=list(A.shape),intensity_units="counts per pixel per 5-second frame, background subtracted",warnings=["Detector correction flags in RASX are retained; count errors are conditional counting estimates and do not include calibration or line-shape systematics","Absolute pole-azimuth zero follows the measured motor phi reference, not a crystallographic RD assignment","No imposed 90 degree spacing, equal widths, fixed 45 degree tilt, or symmetry replication","Unmeasured pole-angle regions remain masked"])
(P/"analysis_results.json").write_text(json.dumps(OUT,indent=2),encoding="utf-8")
compact={k:np.asarray(v,np.float32) for k,v in MAP.items() if k in ["combined_net","combined_variance","combined_coverage","chi38_net","chi42_net","tilt_centers","azimuth_centers"]}
compact.update({k:np.asarray(v[0],np.float32) for k,v in PROFILES.items()});compact.update({k+"_variance":np.asarray(v[1],np.float32) for k,v in PROFILES.items()});compact["phi"]=np.arange(0,360,2)
np.savez_compressed(P/"deliverable_numeric.npz",**compact)
print("RESULTS_START"+json.dumps({"primary":FITS["wide_t30_60_combined"],"chi38":FITS["wide_t30_60_chi38"],"chi42":FITS["wide_t30_60_chi42"],"sensitivity":[{"profile":k,"fwhm":[round(x["fwhm_observed_deg"],3) for x in f]} for k,f in FITS.items() if k.endswith("combined") or "_window" in k],"numeric_bytes":(P/"deliverable_numeric.npz").stat().st_size})+"RESULTS_END")

