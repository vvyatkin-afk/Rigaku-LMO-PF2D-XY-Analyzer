"""Render scientific pole figure and phi fits from measured, unsymmetrized data."""
from pathlib import Path
import sys,json,csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.special import erf
P=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent
D=np.load(P/"deliverable_numeric.npz")
if (P/"analysis_results.json").exists():
 R=json.loads((P/"analysis_results.json").read_text())["fits"]["wide_t30_60_combined"]
else:R=json.loads((P/"primary_results.json").read_text())["primary"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.titlesize":13,"axes.labelsize":11,"pdf.fonttype":42,"ps.fonttype":42,"savefig.facecolor":"white"})
colours=["#347ca9","#9271b5","#e38a3f","#40a48d"]
def wrap(x):return (np.asarray(x)+180)%360-180
def kernel(x,mu,w):
 s=w/(2*np.sqrt(2*np.log(2)))
 return s*np.sqrt(np.pi/2)/2*(erf((x-mu+1)/(np.sqrt(2)*s))-erf((x-mu-1)/(np.sqrt(2)*s)))
def curve(x,r):
 a,m,w,b,k=r["params"];return a*kernel(x,m,w)+b+k*x
# Figure 1: stereographic pole figure, no interpolation or fabricated coverage.
net=D["combined_net"]
cmap=LinearSegmentedColormap.from_list("lmo_intensity",["#fffefd","#e8f1f9","#9dc7e4","#4b94c3","#415d9e","#8b3d8a","#d75253","#f5ab45"])
cmap.set_bad("#ebedef")
fig=plt.figure(figsize=(9.2,9.6))
ax=fig.add_axes([.09,.19,.76,.68],projection="polar")
ax.set_facecolor("#ebedef");ax.set_theta_zero_location("E");ax.set_theta_direction(1)
te=np.deg2rad(np.arange(-1,361,2));re=np.tan(np.deg2rad(np.arange(91))/2)
mesh=ax.pcolormesh(te,re,np.ma.masked_invalid(net),cmap=cmap,vmin=0,vmax=1.65,shading="flat",rasterized=True)
ax.set_ylim(0,1);ax.set_xticks(np.deg2rad(np.arange(0,360,30)));ax.set_xticklabels([f"{x}°" for x in range(0,360,30)])
ticks=np.array([15,30,45,60,75,90]);ax.set_yticks(np.tan(np.deg2rad(ticks)/2));ax.set_yticklabels([f"{x}°" for x in ticks],color="#4b5563",fontsize=9)
ax.set_rlabel_position(35);ax.grid(color="#7d8b98",alpha=.37,lw=.65)
ax.spines["polar"].set_color("#8e9aa5")
# Labels sit outside the measured spots, with real data maxima as an initial anchor.
for i,r in enumerate(R):
 delta=wrap(D["azimuth_centers"]-r["center_deg"]);m=(abs(delta)<12)[None,:]&(D["tilt_centers"][:,None]>30)&(D["tilt_centers"][:,None]<60)
 z=np.where(m,net,np.nan);ii,jj=np.unravel_index(np.nanargmax(z),z.shape)
 rr=np.tan(np.deg2rad(D["tilt_centers"][ii])/2)
 ax.text(np.deg2rad(r["center_deg"]),rr+.09,f"P{i+1}",ha="center",va="center",weight="bold",color="#18324b",fontsize=12,bbox={"facecolor":"white","alpha":.85,"edgecolor":"none","pad":1.8})
cax=fig.add_axes([.89,.30,.025,.50]);cb=fig.colorbar(mesh,cax=cax);cb.set_label("Background-subtracted counts / pixel / 5 s",labelpad=10);cb.ax.tick_params(labelsize=9)
fig.text(.08,.964,"LMO (011) | measured pole figure",fontsize=20,weight="bold",color="#17324d")
fig.text(.08,.935,"Stereographic projection  •  360 HyPix frames  •  no imposed crystal symmetry",fontsize=10.5,color="#536474")
fig.text(.08,.115,"Radial angle t: tilt from the specimen normal (0° at center, 90° at rim).",fontsize=10)
fig.text(.08,.090,"Azimuth β: measured φ reference; 0° at right, increasing counterclockwise.",fontsize=10)
fig.text(.08,.065,"Gray regions were not measured with valid signal and background coverage.",fontsize=10,color="#596774")
fig.text(.08,.040,"2θ signal 31.6–33.3°; background 30.8–31.3° and 33.7–34.2°. Display grid 1° × 2°.",fontsize=9,color="#596774")
fig.savefig(P/"LMO_011_pole_figure.png",dpi=320);fig.savefig(P/"LMO_011_pole_figure.pdf");plt.close(fig)
# Figure 2: measured phi profile and independent local fits.
phi=D["phi"];key="wide_t30_60_combined";y=D[key];sigma=np.sqrt(D[key+"_variance"])
fig=plt.figure(figsize=(12.5,9.0));gs=fig.add_gridspec(3,2,height_ratios=[1.25,1,1],hspace=.55,wspace=.25,left=.08,right=.96,top=.88,bottom=.14)
ax=fig.add_subplot(gs[0,:])
ax.errorbar(phi,y,yerr=sigma,fmt="o",ms=3.1,color="#344a5f",ecolor="#b3bfca",elinewidth=.8,capsize=0,label="Measured intensity (two χ settings)")
for i,r in enumerate(R):
 xx=np.linspace(-r["fit_half_window_deg"],r["fit_half_window_deg"],350);px=xx+r["seed_center_deg"];py=curve(xx,r)
 for shift in [-360,0,360]:
  a=px+shift;sel=(a>=0)&(a<=360);ax.plot(a[sel],py[sel],lw=2.1,color=colours[i])
ax.set_xlim(0,360);ax.set_xticks(np.arange(0,361,30));ax.set_ylabel("Counts / pixel / 5 s");ax.set_xlabel("Sample rotation φ (degrees)");ax.grid(alpha=.18);ax.spines[["top","right"]].set_visible(False)
for i,r in enumerate(R):
 ax=fig.add_subplot(gs[1+i//2,i%2]);x=wrap(phi-r["seed_center_deg"]);sel=abs(x)<=22;order=np.argsort(x[sel]);xx=x[sel][order];yy=y[sel][order];ss=sigma[sel][order]
 ax.errorbar(xx,yy,yerr=ss,fmt="o",ms=3.8,color=colours[i],ecolor="#bdc8d0",capsize=2)
 xf=np.linspace(-22,22,450);ax.plot(xf,curve(xf,r),color=colours[i],lw=2)
 a,m,w,b,k=r["params"];half=b+k*m+a*kernel(m,m,w)/2;fwhm=r["fwhm_observed_deg"]
 ax.plot([m-fwhm/2,m+fwhm/2],[half,half],color="#384c60",lw=1.6)
 ax.plot([m-fwhm/2,m-fwhm/2],[half*.96,half*1.04],color="#384c60");ax.plot([m+fwhm/2,m+fwhm/2],[half*.96,half*1.04],color="#384c60")
 ax.set_title(f"P{i+1}   φ = {r['center_deg']:.2f}°   |   FWHM = {fwhm:.2f} ± {r['fwhm_se_scaled_deg']:.2f}°",fontsize=11,loc="left",color="#17324d")
 ax.set_xlim(-22,22);ax.set_xlabel(f"φ − {r['seed_center_deg']:g}° (wrapped degrees)");ax.set_ylabel("Counts / pixel / 5 s");ax.grid(alpha=.18);ax.spines[["top","right"]].set_visible(False)
fig.text(.08,.953,"LMO (011) | observed azimuthal widths",fontsize=20,weight="bold",color="#17324d")
fig.text(.08,.916,"Four independent circular fits; 2° continuous φ sweep included in the profile model.",fontsize=11,color="#536474")
fig.text(.08,.041,"Integration: t = 30–60°. Gaussian plus linear local background. Uncertainties: residual-scaled 1σ fit estimates.",fontsize=9.3,color="#596774")
fig.text(.08,.020,"These are observed Δφ widths, with no instrumental-resolution subtraction. P4 uses measurements across the 360°/0° boundary.",fontsize=9.3,color="#596774")
fig.savefig(P/"LMO_011_phi_fits.png",dpi=300);fig.savefig(P/"LMO_011_phi_fits.pdf");plt.close(fig)
# Reusable numeric tables, retaining negative background fluctuations.
with (P/"LMO_011_phi_profile.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.writer(f);w.writerow(["phi_deg","chi38_net_counts_per_pixel_5s","chi38_standard_error","chi42_net_counts_per_pixel_5s","chi42_standard_error","combined_net_counts_per_pixel_5s","combined_standard_error"])
 for j,p in enumerate(phi):
  vals=[]
  for s in ["chi38","chi42","combined"]:
   k="wide_t30_60_"+s;vals.extend([float(D[k][j]),float(np.sqrt(D[k+"_variance"][j]))])
  w.writerow([float(p)]+vals)
with (P/"LMO_011_phi_peaks.csv").open("w",newline="",encoding="utf-8") as f:
 fields=["peak","center_deg","center_se_scaled_deg","fwhm_observed_deg","fwhm_se_scaled_deg","gaussian_component_fwhm_deg","reduced_chi_square","fit_half_window_deg"]
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for i,r in enumerate(R):w.writerow({"peak":i+1,**{k:r[k] for k in fields[1:]}})
with (P/"LMO_011_pole_grid.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.writer(f);w.writerow(["tilt_from_normal_deg","azimuth_deg","net_counts_per_pixel_5s","standard_error","pixel_contributions","stereographic_x","stereographic_y"])
 for i,t in enumerate(D["tilt_centers"]):
  rr=np.tan(np.deg2rad(t)/2)
  for j,b in enumerate(D["azimuth_centers"]):
   ok=np.isfinite(net[i,j]);w.writerow([float(t),float(b),float(net[i,j]) if ok else "",float(np.sqrt(D["combined_variance"][i,j])) if ok else "",int(D["combined_coverage"][i,j]),float(rr*np.cos(np.deg2rad(b))),float(rr*np.sin(np.deg2rad(b)))])
print("Rendered pole figure, phi fits, and CSV tables in",P)

