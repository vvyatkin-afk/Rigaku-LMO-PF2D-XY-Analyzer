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
    def test_phi_normalization(self):
        self.assertEqual(normalize_phi(-4.),356.)
        self.assertEqual(normalize_phi(360.),0.)

if __name__=='__main__':unittest.main()
