"""Batch runner for measured LMO (011) pole figures at multiple XY positions.
Each point needs a unique RASX filename prefix or its own folder.
Requires the SAME acquisition as Ivan's verified example: chi=38,42 degrees,
phi=0,2,...,358, exactly 360 images / point. Does not invent missing data.
Original RASX files and Rigaku configuration are never modified.
"""
import argparse,csv,datetime,hashlib,json,pathlib,re,shutil,subprocess,sys,zipfile
from collections import defaultdict
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader,PdfWriter

P=pathlib.Path(__file__).resolve().parent
SCRIPT=P/'validated_single_point_scripts'
REG=re.compile(r'^(?P<stem>.+)_image_(?P<seq>\d+)_(?P<chi>-?\d+(?:\.\d+)?)_(?P<num>\d+)\.rasx$',re.I)
def safe(s):return re.sub(r'[^A-Za-z0-9_-]','_',s)[:58]+'_'+hashlib.sha256(s.encode('utf-8')).hexdigest()[:8]
def findgroups(root):
    groups=defaultdict(list)
    for f in root.rglob('*.rasx'):
        m=REG.match(f.name)
        if m: groups[str(f.parent.relative_to(root))+'__'+m['stem']].append((f,int(m['seq']),float(m['chi']),int(m['num'])))
    return dict(groups)
def axes(p):
    with zipfile.ZipFile(p) as z:r=ET.fromstring(z.read('Data0/MesurementConditions0.xml'))
    return {x.get('Name'):float(x.get('Position')) for x in r.findall('Axes/Axis') if x.get('Name') in ('X','Y') and x.get('Position') is not None}
def checked_group(frames):
    if len(frames)!=360:raise ValueError('Expected 360 frames (180 for each chi), found '+str(len(frames)))
    refs={}
    for p,seq,chi,n in frames:
        if (seq,chi) not in [(1,38.),(2,42.)] or not (1<=n<=180):
            raise ValueError('Unexpected chi/index/frame: '+p.name)
        if (seq,chi,n) in refs:raise ValueError('Duplicated index '+str((seq,chi,n)))
        refs[(seq,chi,n)]=p
    if len(refs)!=360:raise ValueError('Incomplete phi index grid')
    xy=[axes(refs[(seq,chi,n)]) for seq,chi,n in [(1,38.,1),(1,38.,180),(2,42.,1),(2,42.,180)]]
    coords=[(q.get('X',float('nan')),q.get('Y',float('nan'))) for q in xy]
    for q in coords[1:]:
        if any(np.isfinite(v) and np.isfinite(w) and abs(v-w)>0.05 for v,w in zip(coords[0],q)):
            raise ValueError('XY changes within a point group; separate the source datasets')
    return refs,coords[0]
def cover(pid,xy,phi,poles,dest):
    fig,ax=plt.subplots(figsize=(8.27,11.69));ax.axis('off')
    fig.text(.08,.94,'LMO (011) | XY point report',fontsize=21,weight='bold',color='#1e3e5d')
    fig.text(.08,.895,f'{pid}   |   X={xy[0]:.3f}, Y={xy[1]:.3f} mm',fontsize=13)
    fig.text(.08,.865,'Measured: 2 chi positions x 180 phi frames = 360 HyPix frames',fontsize=10)
    h=['Pole','beta (deg)','t (deg)','phi FWHM','beta cut FWHM','t cut FWHM']
    rows=[]
    for i,(a,b) in enumerate(zip(phi,poles),1):
        rows.append([f'P{i}',f"{b['azimuth_center_deg']:.2f}",f"{b['tilt_center_deg']:.2f}",
          f"{a['fwhm_observed_deg']:.2f}",f"{b['azimuth_center_cut_fwhm_deg']:.2f}",
          f"{b['tilt_center_cut_fwhm_deg']:.2f}"])
    tab=ax.table(cellText=rows,colLabels=h,bbox=[.06,.41,.88,.32],cellLoc='center')
    tab.auto_set_font_size(False);tab.set_fontsize(8.8)
    for (i,j),cell in tab.get_celld().items():
        cell.set_edgecolor('#d4e0eb')
        if i==0:cell.set_facecolor('#e7f0f8');cell.set_text_props(weight='bold')
    for i,line in enumerate([
     'FWHM values are observed, not instrument-resolution corrected.',
     'Azimuthal phi FWHM is integrated over t = 30-60 degrees.',
     'The 2D widths are the fitted center cuts (beta and tilt).',
     'Four peaks are fitted independently; no imposed crystal symmetry.',
     'For source geometry, each RASX is checked against the reference.',
     'If the metadata differs, the pipeline rejects the point.']):
        fig.text(.09,.34-.033*i,'- '+line,fontsize=9)
    fig.text(.08,.06,datetime.datetime.now().astimezone().isoformat(timespec='minutes'),fontsize=8,color='#657789')
    fig.savefig(dest/'point_cover.pdf');plt.close(fig)
    w=PdfWriter()
    for p in [dest/'point_cover.pdf',dest/'LMO_011_pole_figure.pdf',dest/'LMO_011_phi_fits.pdf']:
        for page in PdfReader(p).pages:w.add_page(page)
    with (dest/'LMO011_point_report.pdf').open('wb') as f:w.write(f)
def run_one(group,frames,root,work,dry):
    point=safe(group);refs,xy=checked_group(frames)
    if dry:return {'point':point,'x_mm':xy[0],'y_mm':xy[1],'frames':360,'status':'VALID'}
    staged=work/point;source=staged/'frames';scripts=staged/'scripts'
    source.mkdir(parents=True);scripts.mkdir(parents=True)
    for s in ['run_lmo_pf.py','reconstruct_poles.py','fit_poles.py','render_results.py']:
        if not (SCRIPT/s).is_file():raise FileNotFoundError('Prior validated script missing '+s)
        shutil.copy2(SCRIPT/s,scripts/s)
    with (scripts/'reference_manifest.csv').open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=['file','sha256']);wr.writeheader()
        for (seq,chi,n),p in sorted(refs.items()):
            name=f'011_2D_LMO_011_PF2D_TEST_image_{seq}_{chi:.3f}_{n:03d}.rasx'
            dest=source/name;shutil.copy2(p,dest)
            sha=hashlib.sha256(p.read_bytes()).hexdigest()
            if sha!=hashlib.sha256(dest.read_bytes()).hexdigest():raise ValueError('Copy hash mismatch')
            wr.writerow({'file':name,'sha256':sha})
    dest=root/'points'/point
    dest.parent.mkdir(parents=True,exist_ok=True)
    with (staged/'processing.log').open('w',encoding='utf-8') as f:
        cmd=([sys.executable,'--script','run_lmo_pf.py'] if getattr(sys,'frozen',False) else [sys.executable,str(scripts/'run_lmo_pf.py')])
        rc=subprocess.run([*cmd,'--source',str(source),'--output',str(dest)],stdout=f,stderr=subprocess.STDOUT)
    if rc.returncode:raise RuntimeError('Prior validated reconstruction failed; see '+str(staged/'processing.log'))
    p=json.loads((dest/'analysis_results.json').read_text(encoding='utf-8'))['fits']['wide_t30_60_combined']
    poles=json.loads((dest/'pole_2d_fits.json').read_text(encoding='utf-8'))['peaks']
    if len(p)!=4 or len(poles)!=4:raise ValueError('Cannot resolve four poles')
    if not all(q['fwhm_observed_deg']>0 and np.isfinite(q['fwhm_observed_deg']) for q in p):raise ValueError('Invalid FWHM')
    cover(point,xy,p,poles,dest)
    out={'point':point,'x_mm':xy[0],'y_mm':xy[1],'frames':360,'status':'OK',
         'mean_phi_fwhm_deg':float(np.mean([v['fwhm_observed_deg'] for v in p]))}
    for i,(a,b) in enumerate(zip(p,poles),1):
        out.update({f'P{i}_phi_deg':a['center_deg'],f'P{i}_phi_fwhm_deg':a['fwhm_observed_deg'],
         f'P{i}_beta_deg':b['azimuth_center_deg'],f'P{i}_tilt_deg':b['tilt_center_deg'],
         f'P{i}_beta_cut_fwhm_deg':b['azimuth_center_cut_fwhm_deg'],
         f'P{i}_tilt_cut_fwhm_deg':b['tilt_center_cut_fwhm_deg']})
    return out
def summary_pdf(out,rows):
    ok=[r for r in rows if r['status']=='OK']
    with PdfPages(out/'LMO011_XY_Overall_Report.pdf') as pdf:
        fig=plt.figure(figsize=(11.69,8.27));fig.patch.set_facecolor('white');plt.axis('off')
        fig.text(.07,.89,'LMO (011) | Multi-point pole figures',fontsize=21,weight='bold',color='#153d60')
        fig.text(.07,.80,f"Measured point groups: {len(ok)} complete; {len(rows)-len(ok)} rejected or incomplete",fontsize=12)
        for i,t in enumerate(['No extrapolation between XY points; missing points remain missing.',
          'Each completed point has a separate report with four independently fitted peaks.',
          'Observed phi widths retain angular broadening; no instrumental deconvolution.',
          'Required acquisition: chi 38/42 deg, phi step 2 deg, 360 frames per point.']):
            fig.text(.08,.69-.065*i,'- '+t,fontsize=11)
        fig.text(.08,.13,'Offline analysis of archived RASX data; not a native SmartLab Texture project.',fontsize=10,color='#ad493b')
        pdf.savefig(fig);plt.close(fig)
        if ok:
            fig,axs=plt.subplots(1,2,figsize=(11.69,8.27))
            for ax,key,ttl in zip(axs,['mean_phi_fwhm_deg','P1_tilt_deg'],
                                  ['Average phi FWHM (deg)','P1 tilt from sample normal (deg)']):
                xs=[q['x_mm'] for q in ok];ys=[q['y_mm'] for q in ok];zs=[q[key] for q in ok]
                im=ax.scatter(xs,ys,c=zs,cmap='viridis',s=200,edgecolors='white')
                for q in ok:ax.annotate(q['point'],(q['x_mm'],q['y_mm']),xytext=(5,5),textcoords='offset points',fontsize=7)
                ax.set(xlabel='X (mm)',ylabel='Y (mm)',title=ttl);ax.grid(alpha=.2)
                fig.colorbar(im,ax=ax,shrink=.75)
            fig.tight_layout();pdf.savefig(fig);plt.close(fig)
            for k in range(0,len(ok),12):
                part=ok[k:k+12];fig,ax=plt.subplots(figsize=(11.69,8.27));ax.axis('off')
                header=['Point','X','Y','Mean phi FWHM','P1','P2','P3','P4']
                data=[[q['point'],f"{q['x_mm']:.2f}",f"{q['y_mm']:.2f}",f"{q['mean_phi_fwhm_deg']:.2f}"]+
                 [f"{q[f'P{i}_phi_fwhm_deg']:.2f}" for i in range(1,5)] for q in part]
                t=ax.table(cellText=data,colLabels=header,loc='center',cellLoc='center')
                t.auto_set_font_size(False);t.set_fontsize(8.5);t.scale(1,1.5)
                for (i,j),c in t.get_celld().items():
                    c.set_edgecolor('#d4dfe9')
                    if i==0:c.set_facecolor('#e5eff7')
                ax.set_title('Observed phi widths by XY position',fontsize=15,pad=18)
                pdf.savefig(fig);plt.close(fig)
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=pathlib.Path,required=True)
    ap.add_argument('--output',type=pathlib.Path,required=True)
    ap.add_argument('--group-key',default=None,help='Process exactly one selected group')
    ap.add_argument('--dry-run',action='store_true')
    ap.add_argument('--keep-work',action='store_true')
    a=ap.parse_args()
    if not a.source.is_dir():ap.error('Source folder inaccessible')
    if a.output.exists():ap.error('Output folder already exists; select a NEW one')
    groups=findgroups(a.source)
    if a.group_key is not None:
        if a.group_key not in groups:ap.error('Selected group not found: '+a.group_key)
        groups={a.group_key:groups[a.group_key]}
    if not groups:ap.error('No <prefix>_image_<chi_index>_<chi>_<frame>.rasx files')
    a.output.mkdir(parents=True)
    work=a.output/'_staged';work.mkdir()
    rows=[]
    for name,files in groups.items():
        print('GROUP',name,'FRAMES',len(files),flush=True)
        try:row=run_one(name,files,a.output,work,a.dry_run)
        except Exception as e:row={'point':safe(name),'frames':len(files),'status':'FAILED','reason':str(e)};print('ERROR',str(e),file=sys.stderr)
        rows.append(row)
    if not a.dry_run:summary_pdf(a.output,rows)
    headers=sorted(set().union(*(r.keys() for r in rows)))
    with (a.output/'XY_summary.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=headers);w.writeheader();w.writerows(rows)
    (a.output/'XY_summary.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False),encoding='utf-8')
    if not a.keep_work:shutil.rmtree(work,ignore_errors=True)
    if rows and all(r.get('status')=='OK' for r in rows):
        print('RESULTS_READY',a.output,flush=True)
    print('DONE',a.output,'OK',sum(r['status']=='OK' for r in rows),'TOTAL',len(rows))
    if any(r['status']=='FAILED' for r in rows):sys.exit(2)
if __name__=='__main__':main()
