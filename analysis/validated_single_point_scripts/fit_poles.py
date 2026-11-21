"""Fit unsymmetrized measured pole features and test a common tilted orientation.

Input is the independently reconstructed, unsmoothed 1 degree tilt by 2 degree
sample-azimuth map. The effective Gaussian widths retain measurement and map
bin broadening; no instrument or phi-exposure deconvolution is performed.
"""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

BASE=Path.cwd()
OUT=BASE/'geometry'
DATA=BASE/'analysis/deliverable_numeric.npz'
FWHM=2*np.sqrt(2*np.log(2))


def wrap(a):
    return (a+180.)%360.-180.


def pole_vectors(t,b):
    t,b=np.deg2rad(t),np.deg2rad(b)
    return np.stack((np.sin(t)*np.cos(b),np.sin(t)*np.sin(b),np.cos(t)),axis=-1)


def rigid_fit(observed, ideal):
    u,s,vt=np.linalg.svd(ideal.T@observed)
    fix=np.diag([1.,1.,np.linalg.det(vt.T@u.T)])
    R=vt.T@fix@u.T
    pred=ideal@R.T
    angles=np.rad2deg(np.arccos(np.clip(np.sum(pred*observed,axis=1),-1.,1.)))
    axis=R@np.array([0.,0.,1.])
    tilt=np.rad2deg(np.arccos(np.clip(axis[2],-1.,1.)))
    azimuth=np.rad2deg(np.arctan2(axis[1],axis[0]))%360.
    assert np.allclose(R.T@R,np.eye(3),atol=1e-10)
    assert np.isclose(np.linalg.det(R),1.,atol=1e-10)
    return R,pred,angles,axis,tilt,azimuth


def fit_peak(d,reference_beta,azhalf=25.,tilt_limits=(30.,60.)):
    t,b=np.meshgrid(d['tilt_centers'],d['azimuth_centers'],indexing='ij')
    localb=wrap(b-reference_beta)
    mask=(np.abs(localb)<=azhalf)&(t>=tilt_limits[0])&(t<=tilt_limits[1])&(d['combined_coverage']>0)&np.isfinite(d['combined_net'])&(d['combined_variance']>0)
    x=t[mask].astype(float); y=localb[mask].astype(float)
    z=d['combined_net'][mask].astype(float); error=np.sqrt(d['combined_variance'][mask].astype(float))
    imax=np.argmax(z)
    # Parameters: amplitude, tilt center, local-beta center, two marginal sigmas,
    # correlation, background offset, background slopes per degree in tilt/beta.
    p0=[float(z.max()),float(x[imax]),float(y[imax]),2.,4.5,0.,0.,0.,0.]
    lo=[0.,35.,-10.,.3,.5,-.97,-.15,-.015,-.015]
    hi=[10.,55.,10.,9.,16.,.97,.15,.015,.015]
    def gaussian(p,x,y):
        a,tc,bc,st,sb,rho,_,_,_=p
        xx=(x-tc)/st; yy=(y-bc)/sb
        return a*np.exp(-(xx*xx-2*rho*xx*yy+yy*yy)/(2*(1-rho*rho)))
    def model(p,x,y):
        return gaussian(p,x,y)+p[6]+p[7]*(x-45.)+p[8]*y
    def residual(p):
        return (model(p,x,y)-z)/error
    fit=least_squares(residual,p0,bounds=(lo,hi),max_nfev=3000,xtol=1e-11,ftol=1e-11,gtol=1e-11)
    p=fit.x
    dof=len(z)-len(p)
    reduced=float(np.dot(fit.fun,fit.fun)/dof)
    cov=np.linalg.pinv(fit.jac.T@fit.jac)*max(1.,reduced)
    st,sb,rho=p[3:6]
    fac=np.sqrt(1-rho*rho)
    # Delta-method errors retain covariance between sigma and correlation.
    def se_of_gradient(g):
        return float(np.sqrt(max(0.,g@cov@g)))
    gt=np.zeros(9);gt[3]=FWHM
    gb=np.zeros(9);gb[4]=FWHM
    gtc=np.zeros(9);gtc[3]=FWHM*fac;gtc[5]=-FWHM*st*rho/fac
    gbc=np.zeros(9);gbc[4]=FWHM*fac;gbc[5]=-FWHM*sb*rho/fac
    result=dict(reference_beta_deg=reference_beta,tilt_center_deg=float(p[1]),azimuth_center_deg=float((reference_beta+p[2])%360.),tilt_center_se_deg=float(np.sqrt(cov[1,1])),azimuth_center_se_deg=float(np.sqrt(cov[2,2])),tilt_beta_center_covariance_deg2=cov[np.ix_([1,2],[1,2])].tolist(),tilt_marginal_fwhm_deg=float(FWHM*st),tilt_marginal_fwhm_se_deg=se_of_gradient(gt),azimuth_marginal_fwhm_deg=float(FWHM*sb),azimuth_marginal_fwhm_se_deg=se_of_gradient(gb),tilt_center_cut_fwhm_deg=float(FWHM*st*fac),tilt_center_cut_fwhm_se_deg=se_of_gradient(gtc),azimuth_center_cut_fwhm_deg=float(FWHM*sb*fac),azimuth_center_cut_fwhm_se_deg=se_of_gradient(gbc),tilt_azimuth_correlation=float(rho),correlation_se=float(np.sqrt(cov[5,5])),amplitude=float(p[0]),background_at_t45_beta_reference=float(p[6]),background_tilt_slope=float(p[7]),background_beta_slope=float(p[8]),reduced_chi_squared=reduced,n_bins=int(len(z)),degrees_of_freedom=dof,success=bool(fit.success),fit_window_tilt_deg=list(tilt_limits),fit_window_beta_offset_deg=[-azhalf,azhalf],parameters=p.tolist(),parameter_covariance=cov.tolist())
    return result,dict(t=t,b=b,localb=localb,mask=mask,model=model,gaussian=gaussian,p=p,x=x,y=y,z=z,error=error)


def main(input_npz=DATA,output_dir=OUT):
    global DATA,OUT
    DATA=Path(input_npz).expanduser().resolve()
    OUT=Path(output_dir).expanduser().resolve()
    OUT.mkdir(parents=True,exist_ok=True)
    d=np.load(DATA)
    references=[82.5,176.5,271.2,356.7]
    fits=[];internals=[]
    for ref in references:
        fit,intern=fit_peak(d,ref)
        # A smaller window checks local-baseline/model sensitivity, without changing
        # the primary fit or claiming a second independent statistical test.
        sensitive,_=fit_peak(d,ref,azhalf=20.,tilt_limits=(33.,57.))
        fit['smaller_window_changes']={k:float(sensitive[k]-fit[k]) for k in ['tilt_center_deg','azimuth_center_deg','tilt_marginal_fwhm_deg','azimuth_marginal_fwhm_deg','tilt_center_cut_fwhm_deg','azimuth_center_cut_fwhm_deg']}
        fits.append(fit);internals.append(intern)
    notes={'input':str(DATA),'model':'Elliptical Gaussian plus background plane, free centers, marginal sigmas and correlation; weighted least squares using propagated map variances.','sampling':'Observed map: 1 degree polar-tilt bins and 2 degree sample-azimuth bins; original continuous phi exposures span 2 degrees. No deconvolution.','width_definitions':'Marginal widths integrate the fitted Gaussian over the other coordinate. Center-cut widths are straight radial or azimuthal cuts through the fitted Gaussian center and equal marginal widths times sqrt(1-rho^2).','uncertainties':'Conditional one-standard-error fit uncertainties, scaled by max(1,reduced chi squared). They exclude geometry-calibration errors and unmodeled covariance from detector corrections or binning.','orientation_assumptions':'Diagnostic assumes pseudocubic (001) texture and four ideal {011} normals at45 degrees, separated90 degrees. It does not establish phase identity or separate mounting tilt from texture-axis tilt.'}
    (OUT/'pole_2d_fits.json').write_text(json.dumps({'notes':notes,'peaks':fits},indent=2))
    columns=['azimuth_center_deg','tilt_center_deg','azimuth_center_se_deg','tilt_center_se_deg','azimuth_marginal_fwhm_deg','azimuth_marginal_fwhm_se_deg','tilt_marginal_fwhm_deg','tilt_marginal_fwhm_se_deg','azimuth_center_cut_fwhm_deg','azimuth_center_cut_fwhm_se_deg','tilt_center_cut_fwhm_deg','tilt_center_cut_fwhm_se_deg','tilt_azimuth_correlation','reduced_chi_squared']
    with (OUT/'pole_2d_fits.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows({k:p[k] for k in columns} for p in fits)
    fig,axs=plt.subplots(2,4,figsize=(18,8))
    fig.subplots_adjust(left=.06,right=.88,bottom=.09,top=.94,hspace=.34,wspace=.30)
    for col,(fit,intern) in enumerate(zip(fits,internals)):
        t=intern['t'];b=intern['localb'];mask=intern['mask'];p=intern['p'];ref=fit['reference_beta_deg']
        display=np.where(mask,d['combined_net'],np.nan)
        rows=np.where(mask.any(axis=1))[0];cols=np.where(mask.any(axis=0))[0]
        cols=cols[np.argsort(b[0,cols])]
        bx=b[0,cols]+ref;tx=t[rows,0]
        patch=display[np.ix_(rows,cols)]
        image_top=axs[0,col].pcolormesh(bx,tx,patch,shading='nearest',cmap='magma',vmin=0,vmax=1.75)
        X,Y=np.meshgrid(bx-ref,tx)
        peak=intern['gaussian'](p,Y,X)
        axs[0,col].contour(bx,tx,peak,levels=[p[0]/2],colors='cyan',linewidths=1.3)
        axs[0,col].plot(ref+p[2],p[1],'+',color='white',ms=10)
        axs[0,col].set_title(f'Pole near {ref:.0f} deg: data + fit')
        model=intern['model'](p,Y,X)
        errs=np.sqrt(d['combined_variance'][np.ix_(rows,cols)])
        image_bottom=axs[1,col].pcolormesh(bx,tx,(patch-model)/errs,shading='nearest',cmap='RdBu_r',vmin=-4,vmax=4)
        axs[1,col].set_title(f'Residual / sigma; chi2/dof={fit["reduced_chi_squared"]:.2f}')
        for row in [0,1]: axs[row,col].set(xlabel='Sample azimuth (deg, unwrapped)',ylabel='Polar tilt from normal (deg)')
    colorbar_top=fig.add_axes([.91,.57,.014,.30])
    colorbar_bottom=fig.add_axes([.91,.12,.014,.30])
    fig.colorbar(image_top,cax=colorbar_top,label='Net counts / pixel / 5 s')
    fig.colorbar(image_bottom,cax=colorbar_bottom,label='Residual / conditional sigma')
    fig.savefig(OUT/'pole_2d_fit_diagnostics.png',dpi=150)
    # Equal-weight proper rigid rotation, in the right-handed specimen plot basis.
    observed=pole_vectors([p['tilt_center_deg'] for p in fits],[p['azimuth_center_deg'] for p in fits])
    ideal=pole_vectors(np.repeat(45.,4),[90.,180.,270.,360.])
    R,pred,residual,axis,axis_tilt,axis_azimuth=rigid_fit(observed,ideal)
    pred_t=np.rad2deg(np.arccos(np.clip(pred[:,2],-1,1)))
    pred_b=np.rad2deg(np.arctan2(pred[:,1],pred[:,0]))%360.
    orient={'assumption':notes['orientation_assumptions'],'fit':'Equal-weight proper rigid rotation fitted to four unit pole vectors; no data or fitted peak center was modified.','axis_tilt_from_specimen_normal_deg':float(axis_tilt),'axis_azimuth_in_instrument_referenced_plot_deg':float(axis_azimuth),'axis_vector_in_plot_basis':axis.tolist(),'rotation_matrix_ideal_to_measured':R.tolist(),'rms_angular_residual_deg':float(np.sqrt(np.mean(residual**2))),'max_angular_residual_deg':float(residual.max()),'pole_residuals':[dict(observed_tilt_deg=float(fit['tilt_center_deg']),observed_azimuth_deg=float(fit['azimuth_center_deg']),predicted_tilt_deg=float(pt),predicted_azimuth_deg=float(pb),tilt_residual_deg=float(fit['tilt_center_deg']-pt),azimuth_residual_deg=float(wrap(fit['azimuth_center_deg']-pb)),spherical_angular_residual_deg=float(rr)) for fit,pt,pb,rr in zip(fits,pred_t,pred_b,residual)],'interpretation':'This diagnoses the direction of the approximate pseudocubic [001] texture axis relative to the mounted specimen coordinates. It cannot distinguish crystallographic texture tilt from mounting/alignment tilt. Azimuth zero comes from the instrument and is not an independently established tape/crystal direction.'}
    orient['agreement_qualification']='The common tilt explains the main pattern, but the exact 45 degree cubic model leaves a systematic outward radial residual of about 0.3 to 0.4 degree, larger than the conditional centroid fit errors. It is approximate geometric agreement, not an exact crystallographic or calibration validation.'
    orient['unresolved_radial_offset']='A departure from the assumed cubic metric or an absolute angular calibration offset could contribute. This dataset alone does not identify the cause; no correction has been applied.'
    (OUT/'orientation_check.json').write_text(json.dumps(orient,indent=2))
    diagram,ax=plt.subplots(figsize=(7.5,7.8),subplot_kw={'projection':'polar'},constrained_layout=True)
    tt=np.array([p['tilt_center_deg'] for p in fits]);bb=np.array([p['azimuth_center_deg'] for p in fits])
    ideal_circle=pole_vectors(np.repeat(45.,721),np.linspace(0.,360.,721))@R.T
    circle_t=np.rad2deg(np.arccos(np.clip(ideal_circle[:,2],-1,1)))
    circle_b=np.arctan2(ideal_circle[:,1],ideal_circle[:,0])
    ax.plot(circle_b,np.tan(np.deg2rad(circle_t)/2),color='#8795a6',lw=1.2,label='Tilted ideal 45 deg pole circle')
    ax.scatter(np.deg2rad(pred_b),np.tan(np.deg2rad(pred_t)/2),s=110,facecolors='none',edgecolors='#192f46',linewidths=1.5,label='Rigid cubic model centers')
    ax.scatter(np.deg2rad(bb),np.tan(np.deg2rad(tt)/2),s=26,c='#df6c2c',label='Measured fitted centers',zorder=4)
    ax.scatter(np.deg2rad(axis_azimuth),np.tan(np.deg2rad(axis_tilt)/2),s=65,c='#157b9f',marker='D',label='Approximate [001] axis')
    for i,(ti,bi) in enumerate(zip(tt,bb),1):
        ax.annotate(str(i),(np.deg2rad(bi),np.tan(np.deg2rad(ti)/2)),xytext=(7,7),textcoords='offset points',color='#a7440d',fontsize=12)
    ticks=[15.,30.,45.,60.]
    ax.set_rticks(np.tan(np.deg2rad(ticks)/2),labels=[f'{t:.0f} deg' for t in ticks])
    ax.set_ylim(0,np.tan(np.deg2rad(65.)/2));ax.set_rlabel_position(25)
    ax.grid(alpha=.35)
    ax.set_title(f'Orientation diagnostic\nAxis tilt {axis_tilt:.2f} deg; spherical RMS residual {np.sqrt(np.mean(residual**2)):.2f} deg',pad=20)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.08),frameon=False,fontsize=10)
    diagram.savefig(OUT/'orientation_check.png',dpi=170,bbox_inches='tight')
    (OUT/'orientation_diagnostic.txt').write_text(f'Independent geometry diagnostic, not a symmetry constraint.\n\nA proper rigid rotation of four ideal pseudocubic {{011}} directions at 45 degrees yields an approximate [001] axis tilt of {axis_tilt:.3f} degrees from the mounted specimen normal, toward instrument-referenced azimuth {axis_azimuth:.3f} degrees. The RMS spherical center residual is {np.sqrt(np.mean(residual**2)):.3f} degrees (maximum {residual.max():.3f} degrees).\n\nThe unequal radial centers and azimuthal spacings are therefore largely consistent with a tilted texture/sample axis. Mounting tilt and intrinsic texture tilt cannot be separated. The exact cubic 45 degree model leaves an outward radial discrepancy of about 0.3 to 0.4 degree, exceeding the conditional statistical centroid errors. A lattice-metric difference or absolute angular offset remains unresolved. No peak was moved, replicated, symmetrized, or corrected by this diagnostic.\n\nWidths in pole_2d_fits.csv are observed effective widths of the unsmoothed measured map. Marginal widths integrate over the other coordinate; center-cut widths are narrower by sqrt(1-rho^2) when tilt and azimuth are correlated. These radial pole widths are not independent omega rocking-curve widths.\n')
    print(json.dumps({'peaks':[{k:p[k] for k in columns} for p in fits],'orientation':orient},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_npz',nargs='?',default=str(DATA),help='Unsmoothed calibrated numerical map NPZ')
    parser.add_argument('output_dir',nargs='?',default=str(OUT),help='Directory for all generated fit tables, diagnostics and plots')
    args=parser.parse_args()
    main(args.input_npz,args.output_dir)
