import unittest
import lmo_pf2d_gui as gui
class AdaptiveGUIWiring(unittest.TestCase):
 def test_new_auto_version(self):
  self.assertGreaterEqual(tuple(int(i) for i in gui.VERSION.split('.')),(1,1,0))
 def test_adaptive_dispatch(self):
  import inspect
  source=inspect.getsource(gui.run_analyze_mode)
  self.assertIn('adaptive_pf2d.py',source)
 def test_gui_has_adaptive_group_validation(self):
  import inspect
  source=inspect.getsource(gui.scan_groups)
  self.assertIn('adaptive_pf2d.py',source)
if __name__=='__main__':unittest.main()
