"""Scientific and end-to-end synthetic RASX checks. No production data."""
import json
import math
import tempfile
import unittest
import zipfile
from pathlib import Path
from xml.etree.ElementTree import Element,SubElement,tostring
import numpy as np
from pypdf import PdfReader
from scipy.spatial.transform import Rotation
from analysis import lattice_strain as l
from analysis import universal_pf2d as u

H=np.array([[1,0,1],[-1,0,1],[0,1,1],[0,-1,1]],float)

def reciprocal(H,A):
    return H@np.linalg.inv(A)

def simulated_frames(root, hkl="011", center=32.6, phi_step=5, xy=0., suffix="01", continuous=False):
    """Build a small 2D HyPix H geometry with four measured spots.
    Pixel geometry forward projection is the inverse-consistency test, not calibration validation.
    """
    root.mkdir(parents=True,exist_ok=True)
    dimensions=[320,300]; origin=[159.5,149.5,.15,.15]; distance=150.
    chi=45. if hkl=="011" else 54.7356
    row=dict(dimensions=dimensions,detector_origin=origin,detector_distance=distance,
             two_theta=center,omega=center/2,chi=chi)
    t,b,signal,bg1,bg2,tt=u.detector_map(row,True)
    centers=[0.,90.,180.,270.] if hkl=="011" else [45.,135.,225.,315.]
    shifts=[-.13,.06,-.13,.06]
    for i,phi in enumerate(np.arange(0,360,phi_step),1):
        inten=np.full_like(tt,12.)
        phis=np.linspace(phi-phi_step/2,phi+phi_step/2,7) if continuous else [phi]
        for phase in phis:
            for beta,delta in zip(centers,shifts):
                db=(b+phase-beta+180)%360-180
                inten+=1500./len(phis)*np.exp(-.5*((t-chi)/2.8)**2-.5*(db/5.8)**2-.5*((tt-center-delta)/.12)**2)
        image=np.rint(inten).astype("<u4")
        xml=Element("MeasurementConditions");axes=SubElement(xml,"Axes")
        for key,value in (("Phi",phi),("Chi",chi),("TwoTheta",center),("Omega",center/2),("X",0),("Y",xy)):
            SubElement(axes,"Axis",Name=key,Position=str(value))
        scan=SubElement(xml,"ScanInformation")
        for key,value in (("AxisName","Phi"),("Mode","CONTINUOUS" if continuous else "STEP"),
                          ("Start",phi-phi_step/2 if continuous else phi),("Stop",phi+phi_step/2 if continuous else phi)):
            SubElement(scan,key).text=str(value)
        imageinfo=SubElement(xml,"ImageInformation");hdr=SubElement(imageinfo,"DTHeader")
        values=dict(PXD_DETECTOR_DIMENSIONS="320 300",PXD_SPATIAL_DISTORTION_INFO="159.5 149.5 0.15 0.15",
            PXD_SPATIAL_DISTORTION_VECTORS="0 -1 1 0",PXD_GONIO_VALUES=f"270 {center} 0 0 0 150",
            PXD_DETECTOR_IDENTIFICATION="HyPix3000(H)",PXD_DETECTOR_IMAGE_DIRECTION="0",
            PXD_DETECTOR_ACQUISITION_EXPTIME="1000",BYTE_ORDER="little_endian",Data_type="unsigned long int")
        for k,v in values.items():
            pair=SubElement(hdr,"Pair")
            SubElement(pair,"string").text=k;SubElement(pair,"string").text=v
        p=root/f"LMO_{hkl}_image_1_{chi:.3f}_{i:03d}_{suffix}.rasx"
        with zipfile.ZipFile(p,"w",zipfile.ZIP_DEFLATED) as z:
            z.writestr("Data0/MesurementConditions0.xml",tostring(xml))
            z.writestr("Data0/Image0.bin",image.tobytes())

class LatticeMathTests(unittest.TestCase):
    def setUp(self):
        self.cfg=l.settings()
        self.A0=l.basis(self.cfg["cell_A_deg"])

    def test_bragg_and_generic_reference(self):
        for h,exp in (("011",3.88/math.sqrt(2)),("111",3.88/math.sqrt(3)),
                      ("222",3.88/math.sqrt(12)),("200",1.94),("204",3.88/math.sqrt(20))):
            self.assertAlmostEqual(l.reference_d(l.hkl_parse(h),self.cfg),exp)
        self.assertEqual(l.hkl_parse("-1, 0, 2"),(-1,0,2))
        self.assertEqual(l.hint_from_name("LMO_204_image"),"204")
        self.assertEqual(l.hint_from_name("222_2D_Test"),"222")
        self.assertEqual(l.hint_from_name("#3636 sample 123"),"not specified")
        cfg=dict(self.cfg,cell_A_deg=[3.7,3.9,4.1,88.,91.,92.])
        self.assertAlmostEqual(l.reference_d((2,0,0),cfg),1/np.linalg.norm(np.linalg.inv(l.basis(cfg["cell_A_deg"])).T@np.array([2,0,0])))

    def test_pure_rotation_no_strain(self):
        for euler in ([0,0,0],[4,-7,28],[15,20,89]):
            A=Rotation.from_euler("xyz",euler,degrees=True).as_matrix()@self.A0
            out=l.solve_cell(H,reciprocal(H,A),self.cfg)
            self.assertEqual(out["status"],"CONDITIONAL_INDEXED_CELL")
            np.testing.assert_allclose(out["principal_strain_pct"],[0,0,0],atol=1e-10)
            np.testing.assert_allclose(out["rotation"],A@np.linalg.inv(self.A0),atol=1e-10)

    def test_isotropic_sensitivity_and_anisotropy(self):
        # Exact reciprocal-vector inversion over requested strains. This is not
        # actual angular or instrument measurement precision.
        for e in (.0001,.0005,.001,.005,.01,-.0001,-.0005,-.001,-.005,-.01):
            A=(1+e)*self.A0
            out=l.solve_cell(H,reciprocal(H,A),self.cfg)
            self.assertAlmostEqual(out["strain_a_pct"],e*100,places=9)
            self.assertAlmostEqual(out["volume_change_pct"],100*((1+e)**3-1),places=8)
        R=Rotation.from_euler("xyz",[4,8,19],degrees=True).as_matrix()
        E=np.diag([-.004,.002,.006])
        A=R@(np.eye(3)+E)@self.A0
        out=l.solve_cell(H,reciprocal(H,A),self.cfg)
        np.testing.assert_allclose([out[k] for k in ("strain_a_pct","strain_b_pct","strain_c_pct")],[-.4,.2,.6],atol=1e-9)
        np.testing.assert_allclose(out["principal_strain_pct"],[-.4,.2,.6],atol=1e-9)
        np.testing.assert_allclose(out["stretch_strain_specimen_pct"],100*R@E@R.T,atol=1e-9)

    def test_shear_and_tape_axis(self):
        cfg=dict(self.cfg,tape_azimuth_deg=30.)
        F=np.array([[1,.01,0],[.01,1,0],[0,0,.995]])
        out=l.solve_cell(H,reciprocal(H,F@self.A0),cfg)
        np.testing.assert_allclose(out["stretch_strain_specimen_pct"],100*(F-np.eye(3)),atol=1e-9)
        a=np.deg2rad(30);v=np.array([np.cos(a),np.sin(a),0])
        self.assertAlmostEqual(out["longitudinal_strain_pct"],100*v@(F-np.eye(3))@v,places=8)
        self.assertNotAlmostEqual(out["gamma_deg"],90.,places=3)

    def test_rank_and_inconsistent_indexing(self):
        g=reciprocal(H,self.A0)
        self.assertEqual(l.solve_cell(H[:2],g[:2],self.cfg)["status"],"INSUFFICIENT_INDEXED_VECTORS")
        h=np.array([[1,1,1],[2,2,2],[3,3,3]])
        self.assertEqual(l.solve_cell(h,reciprocal(h,self.A0),self.cfg)["status"],"RANK_OR_CONDITION_FAILURE")
        bad=g.copy();bad[3]*=1.2
        self.assertEqual(l.solve_cell(H,bad,self.cfg)["status"],"INCONSISTENT_INDEXED_VECTORS")
        bad=g.copy();bad[:,0]*=-1
        self.assertEqual(l.solve_cell(H,bad,self.cfg)["status"],"INVALID_RECIPROCAL_BASIS")

    def test_radial_gaussian_and_calibration_sensitivity(self):
        for center in (32.6,38.9,86.9):
            x=l.radial_axis(center);sigma=.11
            y=8+.3*(x-center)+40*np.exp(-.5*((x-center-.137)/sigma)**2)
            r,model=l.fit_radial(x,y,np.ones_like(x)*40,self.cfg)
            self.assertEqual(r["radial_status"],"OK",r)
            self.assertAlmostEqual(r["observed_two_theta_deg"],center+.137,places=5)
            self.assertAlmostEqual(r["radial_fwhm_deg"],2.354820045*sigma,places=5)
            # A systematic angular offset must change d; pure peak center alone
            # cannot constrain it.
            t=r["observed_two_theta_deg"]
            d1=self.cfg["wavelength_A"]/(2*np.sin(np.deg2rad(t/2)))
            d2=self.cfg["wavelength_A"]/(2*np.sin(np.deg2rad((t+.01)/2)))
            self.assertLess(d2,d1)

    def test_doublet_and_bad_fits(self):
        cfg=dict(self.cfg,secondary_wavelength_A=1.54439,secondary_intensity_ratio=.5)
        x=l.radial_axis(86.9);mu=86.94
        mu2=np.rad2deg(2*np.arcsin(cfg["secondary_wavelength_A"]/cfg["wavelength_A"]*np.sin(np.deg2rad(mu/2))))
        y=10+100*np.exp(-.5*((x-mu)/.09)**2)+50*np.exp(-.5*((x-mu2)/.09)**2)
        r,_=l.fit_radial(x,y,np.ones_like(x),cfg)
        self.assertEqual(r["radial_status"],"OK",r)
        self.assertAlmostEqual(r["observed_two_theta_deg"],mu,places=5)
        y=10+100*np.exp(-.5*((x-x[-1])/.1)**2)
        r,_=l.fit_radial(x,y,np.ones_like(x),self.cfg)
        self.assertNotEqual(r["radial_status"],"OK")
        gapped=np.ones_like(x);gapped[(x>86.85)&(x<87.0)]=0
        r,_=l.fit_radial(x,10+100*np.exp(-.5*((x-86.9)/.1)**2),gapped,self.cfg)
        self.assertNotEqual(r["radial_status"],"OK")

    def test_old_maps_and_invalid_config(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            self.assertEqual(l.analyze_profiles(d,{"sample":"LMO_011"},[],self.cfg)["status"],"RAW_RADIAL_DATA_REQUIRED")
            p=d/"bad.json";p.write_text(json.dumps({"wavelength_A":0}))
            with self.assertRaises(ValueError):l.settings(p)

    def test_auto_index_4_poles_hypothesis(self):
        g=reciprocal(H,self.A0)
        ps=[]
        for i,v in enumerate(g,1):
            q=v/np.linalg.norm(v)
            t=np.rad2deg(np.arccos(q[2]));b=np.rad2deg(np.arctan2(q[1],q[0]))%360
            ps.append(dict(rank=i,g_inv_A=v.tolist(),tilt_deg=t,vector_beta_deg=b,pole_coverage_qc="OK"))
        az=l.infer_four_pole_orientation(ps,(0,1,1))
        self.assertIsNotNone(az)
        cfg=dict(self.cfg,cube_texture_azimuth_deg=az,assume_single_domain=True)
        idx=l.auto_index(ps,(0,1,1),cfg)
        self.assertEqual(len(idx),4)
        self.assertEqual(l.auto_index(ps,(0,1,1),self.cfg),{})
        ps[0]["pole_coverage_qc"]="CLIPPED_OR_GAPPED"
        self.assertIsNone(l.infer_four_pole_orientation(ps,(0,1,1)))

class LatticeRawIntegrationTests(unittest.TestCase):
    def test_end_to_end_one_pass_reports_cache(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);raw=root/"raw";out=root/"out"
            simulated_frames(raw,phi_step=5)
            original=list(raw.glob("*.rasx"))
            sizes=[p.stat().st_size for p in original]
            code=u.run(raw,out)
            self.assertEqual(code,0)
            meta=json.loads((out/"inspection.json").read_text())
            self.assertEqual(len(meta["groups"]),1)
            row=meta["groups"][0];point=out/row["folder"]
            self.assertEqual(row["status"],"OK")
            self.assertFalse((point/"Lattice_Analysis_Error.txt").exists())
            self.assertFalse((point/"Lattice_Report_Error.txt").exists())
            self.assertFalse((out/"Batch_Lattice_Report_Error.txt").exists())
            good=[p for p in row["lattice_analysis"]["peaks"] if p["radial_status"]=="OK"]
            self.assertGreaterEqual(len(good),4,row["lattice_analysis"])
            for p in good:
                self.assertAlmostEqual(p["corrected_two_theta_deg"],32.47 if abs((p["vector_beta_deg"]+45)%180-45)<10 else 32.66,delta=.07)
                self.assertTrue(2.6<p["d_A"]<2.9)
            for name in ("Radial_Cube.npz","Radial_Profiles.npz","Reflection_dSpacing.csv",
                         "Lattice_Parameters.csv","Lattice_Strain.csv","Strain_Tensor.json",
                         "peak_profile_cuts.csv","structure_tilt.json","Pole_Figure_Report.pdf"):
                self.assertTrue((point/name).is_file(),name)
            self.assertGreater(len(PdfReader(point/"Pole_Figure_Report.pdf").pages),6)
            self.assertTrue((out/"Batch_Lattice_Summary.csv").is_file())
            self.assertGreater(len(PdfReader(out/"Overall_Pole_Figure_Report.pdf").pages),5)
            new=root/"rebuild";cfg=dict(l.settings(),two_theta_offset_deg=.02)
            self.assertEqual(u.rebuild(out,new,cfg),0)
            meta2=json.loads((new/"inspection.json").read_text())
            good2=[p for p in meta2["groups"][0]["lattice_analysis"]["peaks"] if p["radial_status"]=="OK"]
            self.assertAlmostEqual(good2[0]["corrected_two_theta_deg"]-good[0]["corrected_two_theta_deg"],.02,places=6)
            self.assertEqual(sizes,[p.stat().st_size for p in original])
            with self.assertRaises(ValueError):u.rebuild(out,root/"bad",dict(cfg,detector_distance_offset_mm=1.))

if __name__=="__main__":unittest.main()
