"""Universal PF2D scientific and source-integrity regression tests (synthetic data only)."""
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
import numpy as np

from analysis.universal_pf2d import (
    FRAME_RE, detector_map, fit_peaks, group_frames,
    series_summary, _safe_out, parse_frame)

class UniversalPF2DTests(unittest.TestCase):
    def test_phi_steps_complete_and_partial(self):
        for step in (.5,2.,4.,5.,10.):
            n=round(360/step)
            rows=[dict(phi=(i*step+.25)%360,mode="CONTINUOUS",
                       exposure_ms=500.,sweep_deg=step) for i in range(n)]
            full=series_summary(rows)
            self.assertTrue(full["full_turn"],(step,full))
            del rows[5]
            incomplete=series_summary(rows)
            self.assertFalse(incomplete["full_turn"])
            self.assertTrue(incomplete["phi_gaps"])
    def test_one_chi_and_mixed_exposure(self):
        rows=[dict(phi=i*90.,mode="CONTINUOUS",sweep_deg=90.,exposure_ms=x)
              for i,x in enumerate((500.,2000.,500.,2000.))]
        s=series_summary(rows)
        self.assertEqual(s["exposure_ms_min"],500.)
        self.assertEqual(s["exposure_ms_max"],2000.)
        self.assertTrue(s["full_turn"])
    def test_hkl_and_xy_separate(self):
        rows=[]
        for name,tt,x in (("222_2D_Test_image_1_54.700_001_07.rasx",86.9,0),
                          ("111_2D_Test_image_1_54.700_001_07.rasx",38.4,0),
                          ("222_2D_Test_image_1_54.700_001_07.rasx",86.9,10)):
            rows.append(dict(path="/tmp/"+name,sample=FRAME_RE.match(Path(name).stem).group("sample"),
                 x=x,y=4.,acquisition_suffix="07",two_theta=tt,omega=tt/2,
                 detector="HyPix3000(H)",dimensions=[775,385]))
        self.assertEqual(len(group_frames(rows)),3)
    def test_geometry_general_2theta(self):
        base=dict(dimensions=[775,385],detector_origin=[382.5,197.5,.1,.1],
                  detector_distance=149.842,omega=44.3544,chi=54.7,two_theta=86.9)
        tilt,azi,signal,bg1,bg2=detector_map(base)
        self.assertEqual(signal.shape,(385,775))
        self.assertGreater(signal.sum(),120)
        self.assertGreater(bg1.sum(),50)
        self.assertGreater(bg2.sum(),50)
        self.assertTrue(np.isfinite(tilt).all())
        base.update(two_theta=32.6,omega=16.3,chi=38.0)
        self.assertGreater(detector_map(base)[2].sum(),120)
    def test_unmeasured_sectors_not_filled_and_dynamic_peaks(self):
        nt,na=90,180
        tt=(np.arange(nt)+.5)[:,None]
        aa=(np.arange(na)+.5)[None,:]*2.
        delta=((aa-124+180)%360)-180
        z=2*np.exp(-.5*((tt-54)/3)**2-.5*(delta/5)**2)
        cov=np.ones_like(z)
        cov[:,90:]=0
        data=np.where(cov>0,z,np.nan)
        peaks=fit_peaks(data,cov,n_max=8)
        self.assertGreaterEqual(len(peaks),1)
        self.assertAlmostEqual(peaks[0]["beta_deg"],124,delta=4)
        self.assertAlmostEqual(peaks[0]["tilt_deg"],54,delta=3)
        self.assertLess(len(peaks),8)
    def test_unresolved_edge_spikes_are_not_reported_as_peaks(self):
        d=np.full((90,180),-0.006)
        c=np.zeros_like(d)
        c[41:65,:]=1.
        for phi in (20,65,100,146,212,277,330):
            d[41,phi//2]=.45
        d[c==0]=np.nan
        peaks=fit_peaks(d,c)
        self.assertEqual(peaks,[],peaks)

    def test_four_well_resolved_poles(self):
        t=(np.arange(90)+.5)[:,None]
        a=(np.arange(180)+.5)[None,:]*2.
        d=np.zeros((90,180))
        for phi in (82.,177.,272.,357.):
            delta=((a-phi+180.)%360.)-180.
            d+=.3*np.exp(-.5*((t-47.)/2.5)**2-.5*(delta/4.5)**2)
        c=np.ones_like(d)
        peaks=fit_peaks(d,c,n_max=10)
        fitted=[p["beta_deg"] for p in peaks if p["status"]=="OK"]
        self.assertEqual(len(fitted),4,fitted)
        for wanted in (82.,177.,272.,357.):
            self.assertTrue(any(abs(((x-wanted+180)%360)-180)<3. for x in fitted),fitted)

    def test_refuse_output_inside_source(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)/"raw"
            root.mkdir()
            with self.assertRaises(ValueError):
                _safe_out(root,root/"new_results")
            self.assertFalse((root/"new_results").exists())
    def test_realistic_synthetic_rasx_header_without_real_data(self):
        from xml.etree.ElementTree import Element,SubElement,tostring
        root=Element("MeasurementConditions")
        axes=SubElement(root,"Axes")
        for k,v in (("Phi",0),("Chi",54.7),("TwoTheta",86.9),("Omega",44.3544),
                    ("X",0),("Y",4)):
            SubElement(axes,"Axis",Name=k,Position=str(v))
        scan=SubElement(root,"ScanInformation")
        for k,v in (("AxisName","Phi"),("Mode","CONTINUOUS"),("Start","-0.25"),
                    ("Stop","0.25"),("Step","0.5")):
            SubElement(scan,k).text=str(v)
        image=SubElement(root,"ImageInformation")
        h=SubElement(image,"DTHeader")
        values=dict(PXD_DETECTOR_DIMENSIONS="775 385",
            PXD_SPATIAL_DISTORTION_INFO="382.5 197.5 0.1 0.1",
            PXD_SPATIAL_DISTORTION_VECTORS="0 -1 1 0",
            PXD_GONIO_VALUES="270 86.9 0 0 0 149.842",
            PXD_DETECTOR_IDENTIFICATION="HyPix3000(H)",
            PXD_DETECTOR_IMAGE_DIRECTION="0",
            PXD_DETECTOR_ACQUISITION_EXPTIME="500",
            BYTE_ORDER="little_endian",Data_type="unsigned long int")
        for k,v in values.items():
            pair=SubElement(h,"Pair")
            SubElement(pair,"string").text=k
            SubElement(pair,"string").text=v
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/"222_2D_Test_image_1_54.700_001_07.rasx"
            with zipfile.ZipFile(p,"w",compression=zipfile.ZIP_DEFLATED) as z:
                z.writestr("Data0/MesurementConditions0.xml",tostring(root))
                z.writestr("Data0/Image0.bin",bytes(775*385*4))
            r=parse_frame(p)
            self.assertEqual(r["two_theta"],86.9)
            self.assertEqual(r["reflection_hint"],"222")
            self.assertEqual(r["acquisition_suffix"],"07")
            self.assertEqual(r["dimensions"],[775,385])

if __name__=="__main__":
    unittest.main()
