"""Regression tests for numeric XY exports, four-pole tilt and rich cached-map reports."""
import csv
import hashlib
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pypdf import PdfReader
from analysis.adaptive_pf2d import rebuild_reports_only, report
from analysis.rich_reports import estimate_structure_tilt, summarize_records, numeric_position


def synthetic_map():
    beta=np.arange(180)*2.+1.
    tilt=np.arange(90)+.5
    b,t=np.meshgrid(beta,tilt)
    data=np.full((90,180),.025)
    for seed,tc,sb,st,amp in [
        (82.,42.,4.5,2.5,1.3),(177.,44.,5.2,2.6,1.15),
        (272.,46.,4.2,2.3,1.2),(357.,43.,3.9,2.1,1.25)]:
        delta=(b-seed+180)%360-180
        data+=amp*np.exp(-.5*((delta/sb)**2+((t-tc)/st)**2))
    coverage=np.ones_like(data)
    return data,coverage


def fake_acquisition():
    return dict(n_frames=360,per_chi={
        "38.0":dict(n_frames=180,mode="continuous",phi_step_deg=2.,observed_sweep_deg=2.),
        "42.0":dict(n_frames=180,mode="continuous",phi_step_deg=2.,observed_sweep_deg=2.)},
        xy=dict(x=-30.5,y=4.),method="regression-only synthetic map")


class RichReportTests(unittest.TestCase):
    def test_numeric_excel_position(self):
        self.assertEqual(numeric_position("01"),1)
        self.assertEqual(numeric_position("39"),39)
        self.assertIsNone(numeric_position(None))

    def test_cone_normal_recovers_tilt_from_four_polar_positions(self):
        a=math.radians(5.)
        b=math.radians(35.)
        normal=np.array([math.sin(a)*math.cos(b),math.sin(a)*math.sin(b),math.cos(a)])
        aux=np.array([0.,0.,1.])
        e1=np.cross(normal,aux);e1/=np.linalg.norm(e1)
        e2=np.cross(normal,e1)
        poles=[]
        for theta in np.deg2rad([0,90,180,270]):
            q=normal*math.cos(math.radians(45.))+(e1*math.cos(theta)+e2*math.sin(theta))*math.sin(math.radians(45.))
            poles.append(dict(status="OK",beta_deg=math.degrees(math.atan2(q[1],q[0]))%360.,
                              tilt_deg=math.degrees(math.acos(q[2]))))
        result=estimate_structure_tilt(poles)
        self.assertEqual(result["status"],"OK")
        self.assertAlmostEqual(result["axis_tilt_deg"],5.,places=5)
        self.assertAlmostEqual(result["axis_azimuth_deg"],35.,places=5)
        self.assertLess(result["cone_rms_deg"],1e-5)
        poles[1]["status"]="FIT_FAILED"
        self.assertEqual(estimate_structure_tilt(poles)["status"],"INSUFFICIENT_PEAKS")

    def test_39_position_table_integers_and_four_radial_fields(self):
        records=[]
        for i in range(1,40):
            peaks=[dict(status="OK",beta_deg=90*k,tilt_deg=40+k,
                        phi_fwhm_observed_deg=6+k,tilt_fwhm_observed_deg=3+k,
                        fit_r2=.98) for k in range(4)]
            records.append(dict(point=i,position_suffix=f"{i:02d}",x_mm=i,
                                y_mm=0,frames=360,status="OK",peaks=peaks))
        summary=summarize_records(records)
        self.assertEqual(len(summary),39)
        self.assertEqual(summary[0]["position_suffix"],1)
        self.assertEqual(summary[-1]["position_suffix"],39)
        for i in range(1,5):
            self.assertIn(f"P{i}_tilt_deg",summary[0])
            self.assertIn(f"P{i}_tilt_fwhm_deg",summary[0])
            self.assertIn(f"P{i}_stereo_radius",summary[0])

    def test_reports_only_rebuild_uses_cached_maps_and_produces_rich_pdfs(self):
        data,cov=synthetic_map()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            original=root/"already_processed";original.mkdir()
            source=original/"point_001";source.mkdir()
            np.savez_compressed(source/"map.npz",intensity=data,coverage=cov)
            (source/"acquisition.json").write_text(json.dumps(fake_acquisition()),encoding="utf-8")
            metadata=dict(groups=[dict(point=1,position_suffix="01",x_mm=-30.5,y_mm=4.,
                                       status="OK",frames=360,source_group="synthetic")],
                          rejected_files=[])
            (original/"inspection.json").write_text(json.dumps(metadata),encoding="utf-8")
            before=hashlib.sha256((source/"map.npz").read_bytes()).hexdigest()
            target=root/"improved_reports"
            rc=rebuild_reports_only(original,target)
            self.assertEqual(rc,0)
            after=hashlib.sha256((source/"map.npz").read_bytes()).hexdigest()
            self.assertEqual(before,after,"Cached source map was modified")
            point=target/"point_001"
            self.assertTrue((point/"peak_profile_cuts.csv").is_file())
            self.assertTrue((point/"structure_tilt.json").is_file())
            peak=json.loads((point/"peaks.json").read_text())
            self.assertEqual(len(peak),4)
            self.assertTrue(all(p["status"]=="OK" for p in peak))
            self.assertTrue(all(p["fit_r2"]>.94 for p in peak))
            self.assertEqual(len(PdfReader(str(point/"LMO011_adaptive_pole_figure.pdf")).pages),4)
            self.assertGreaterEqual(len(PdfReader(str(target/"LMO011_XY_Overall_Report.pdf")).pages),8)
            with (target/"XY_summary.csv").open(encoding="utf-8-sig",newline="") as f:
                row=next(csv.DictReader(f))
            self.assertEqual(row["position_suffix"],"1")
            self.assertEqual(row["point"],"1")
            self.assertTrue(float(row["P4_tilt_deg"])>40)
            self.assertTrue(float(row["P1_tilt_fwhm_deg"])>0)
            with (point/"peak_profile_cuts.csv").open(encoding="utf-8-sig") as f:
                prof=list(csv.DictReader(f))
            self.assertTrue(any(r["axis"]=="tilt" for r in prof))
            self.assertTrue(any(r["axis"]=="beta" for r in prof))


if __name__=="__main__":
    unittest.main()
