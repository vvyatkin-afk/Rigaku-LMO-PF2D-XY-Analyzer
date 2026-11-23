"""The normal workflow must discover RASX metadata only once per execution."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import analysis.adaptive_pf2d as core


class CombinedRunTests(unittest.TestCase):
    @staticmethod
    def fake_groups():
        results={}
        for i in range(1,4):
            name=f"sample_image_1_38.000_001_{i:02d}.rasx"
            key=("/tmp/rasx","sample",(-30.0+i,0.0),f"{i:02d}")
            results[key]=[dict(path="/tmp/rasx/"+name,x=-30.0+i,y=0.0,chi=38.,phi=0.)]
        return results

    def test_one_discovery_for_all_positions_and_live_events(self):
        groups=self.fake_groups()
        info={"38.0": {"n_frames":1, "mode":"stationary",
                        "phi_step_deg":4., "observed_sweep_deg":0.}}
        callback=lambda progress: (progress(dict(stage="enumerate",found=3)),
                                    progress(dict(stage="headers",total=3,checked=3,
                                                 accepted=3,rejected=0)),
                                    progress(dict(stage="group_done",groups_total=3)))
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/"batch"
            output=io.StringIO()
            with patch.object(core,"discover",side_effect=lambda source,progress: (callback(progress),groups,[])[1:]) as discovery, \
                 patch.object(core,"validate_group",return_value=({38.:[{}]},info)) as validation, \
                 patch.object(core,"reconstruct",return_value=(np.ones((90,180)),np.ones((90,180)),info)) as reconstruction, \
                 patch.object(core,"report",return_value=[{"status":"OK"}]*4) as fits, \
                 patch.object(core,"write_batch_summary",return_value=[]), \
                 contextlib.redirect_stdout(output):
                result=core.run("/tmp/rasx",dest)
            self.assertEqual(result,0)
            self.assertEqual(discovery.call_count,1)
            self.assertEqual(validation.call_count,3)
            self.assertEqual(reconstruction.call_count,3)
            self.assertEqual(fits.call_count,3)
            events=[json.loads(line.split(":",1)[1]) for line in output.getvalue().splitlines()
                    if line.startswith("PROGRESS_JSON:")]
            self.assertEqual(len([e for e in events if e["stage"]=="point_done"]),3)
            self.assertIn("metadata_done",[e["stage"] for e in events])
            self.assertIn("file_scan",[e["stage"] for e in events])
            self.assertIn("complete",[e["stage"] for e in events])

    def test_reconstruction_reports_frame_count_and_filename(self):
        shape=(385,775)
        tilt=np.full(shape,42.,dtype=float)
        az=np.zeros(shape,dtype=float)
        signal=np.zeros(shape,dtype=bool);signal[0,0]=True
        bg1=np.zeros(shape,dtype=bool);bg1[0,1]=True
        bg2=np.zeros(shape,dtype=bool);bg2[0,2]=True
        image=np.zeros(shape,dtype=np.uint32)
        image[0,0]=10;image[0,1]=1;image[0,2]=1
        rows=[dict(path=f"/tmp/frame_{i:03d}.rasx",source_name=f"frame_{i:03d}.rasx",
                   phi=i*4.,chi=38.,two_theta=32.6,mode="STATIONARY",
                   sweep_deg=0.) for i in (1,2)]
        emitted=[]
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(core,"validate_group",return_value=({38.:rows},{})), \
                 patch.object(core,"geometry",return_value=(tilt,az,signal,bg1,bg2)), \
                 patch.object(core,"_image",return_value=image):
                data,coverage,_=core.reconstruct(rows,Path(tmp),progress=lambda n,total,name:
                    emitted.append((n,total,name)))
            self.assertTrue((Path(tmp)/"map.npz").exists())
        self.assertTrue(emitted)
        self.assertEqual(emitted[-1],(2,2,"frame_002.rasx"))

    def test_single_position_does_not_reconstruct_other_positions(self):
        groups=self.fake_groups()
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/"single"
            with patch.object(core,"discover",return_value=(groups,[])) as discovery, \
                 patch.object(core,"validate_group",return_value=({38.:[{}]},{"38.0":{}})), \
                 patch.object(core,"reconstruct",return_value=(np.ones((90,180)),np.ones((90,180)),{})) as reconstruct, \
                 patch.object(core,"report",return_value=[{"status":"OK"}]*4), \
                 patch.object(core,"write_batch_summary"):
                with contextlib.redirect_stdout(io.StringIO()):
                    rc=core.run("/tmp/rasx",dest,position=2)
            self.assertEqual(rc,0)
            self.assertEqual(discovery.call_count,1)
            self.assertEqual(reconstruct.call_count,1)
            self.assertTrue((dest/"point_002").is_dir())
            self.assertFalse((dest/"point_001").exists())


if __name__=="__main__":
    unittest.main()
