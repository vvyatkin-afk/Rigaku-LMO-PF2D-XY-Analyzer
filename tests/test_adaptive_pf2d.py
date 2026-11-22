import math
import unittest
from analysis.adaptive_pf2d import infer_series, parse_rasx_header, normalize_phi, select_point_groups

class AdaptiveAcquisitionTests(unittest.TestCase):
    def test_continuous_infer_step_and_frame_count(self):
        rows=[dict(phi=float(i*2),chi=38.,mode='CONTINUOUS',sweep_deg=2.,x=0.,y=0.) for i in range(180)]
        p=infer_series(rows)
        self.assertEqual(p['n_frames'],180)
        self.assertEqual(p['phi_step_deg'],2.)
        self.assertEqual(p['mode'],'continuous')
        self.assertTrue(p['full_turn'])
    def test_stationary_four_degree(self):
        rows=[dict(phi=float(i*4),chi=42.,mode='SINGLE',sweep_deg=0.,x=0.,y=0.) for i in range(90)]
        p=infer_series(rows)
        self.assertEqual(p['phi_step_deg'],4.)
        self.assertEqual(p['n_frames'],90)
        self.assertTrue(p['full_turn'])
        self.assertEqual(p['mode'],'stationary')
    def test_missing_frame_is_not_silently_interpolated(self):
        rows=[dict(phi=float(i*4),chi=38.,mode='SINGLE',sweep_deg=0.) for i in range(90) if i!=19]
        p=infer_series(rows)
        self.assertFalse(p['complete_grid'])
        self.assertEqual(len(p['missing_phi_deg']),1)
    def test_wraparound_360_duplicate(self):
        rows=[dict(phi=float(i*4),chi=38.,mode='SINGLE',sweep_deg=0.) for i in range(91)]
        p=infer_series(rows)
        self.assertTrue(p['full_turn'])
        self.assertEqual(p['duplicate_canonical_angles'],[0.])
    def test_different_xy_not_mixed(self):
        rows=[dict(phi=0.,chi=38.,x=0.,y=0.,path='a'),dict(phi=4.,chi=38.,x=5.,y=0.,path='b')]
        groups=select_point_groups(rows)
        self.assertEqual(len(groups),2)
    def test_flat_folder_position_suffixes_keep_39_points_separate(self):
        rows=[]
        for position in range(1,40):
            for i in range(90):
                rows.append(dict(phi=float(i*4),chi=38.,x=position*5.,y=0.,
                    mode='SINGLE',sweep_deg=0.,
                    path=f'/data/Ivan/batch/011_2D_LMO_011_batch_image_1_38.000_{i+1:03d}_{position:02d}.rasx'))
        groups=select_point_groups(rows)
        self.assertEqual(len(groups),39)
        self.assertTrue(all(len(group)==90 for group in groups.values()))
        self.assertEqual(len({key[3] for key in groups}),39)
        self.assertTrue(all(infer_series(group)['complete_grid'] for group in groups.values()))

    def test_original_filename_without_position_suffix_still_groups(self):
        rows=[dict(phi=float(i*4),chi=38.,x=0.,y=0.,
            path=f'/data/P001_LMO011_image_1_38.000_{i+1:03d}.rasx')
            for i in range(90)]
        groups=select_point_groups(rows)
        self.assertEqual(len(groups),1)
        self.assertIsNone(next(iter(groups))[3])

    def test_discover_progress_snapshots_include_failures(self):
        from pathlib import Path
        from tempfile import TemporaryDirectory
        from unittest.mock import patch
        from analysis.adaptive_pf2d import discover

        with TemporaryDirectory() as tmp:
            root=Path(tmp)
            for filename in ("sample_image_1_38.000_001_01.rasx",
                             "sample_image_1_38.000_002_01.rasx",
                             "bad.rasx"):
                (root/filename).touch()
            def fake_header(path):
                if path.name=="bad.rasx":
                    raise ValueError("bad header")
                return dict(x=0.,y=0.,chi=38.,phi=0.)
            snapshots=[]
            with patch("analysis.adaptive_pf2d.parse_rasx_header",side_effect=fake_header):
                groups,errors=discover(root,progress=snapshots.append)
        self.assertEqual(len(groups),1)
        self.assertEqual(len(errors),1)
        self.assertIn("enumerate",[s['stage'] for s in snapshots])
        self.assertIn("headers",[s['stage'] for s in snapshots])
        final=snapshots[-1]
        self.assertEqual(final['stage'],'group_done')
        self.assertEqual(final['total'],3)
        self.assertEqual(final['checked'],3)
        self.assertEqual(final['accepted'],2)
        self.assertEqual(final['rejected'],1)
        self.assertIn("sample_image",final['last_ok'])

    def test_phi_normalization(self):
        self.assertEqual(normalize_phi(-4.),356.)
        self.assertEqual(normalize_phi(360.),0.)

if __name__=='__main__':unittest.main()
