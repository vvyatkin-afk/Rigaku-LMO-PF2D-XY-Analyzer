import unittest
import lmo_pf2d_gui as gui

class ProgressMessageTests(unittest.TestCase):
 def test_file_progress_includes_stage_counts_last_checked_and_errors(self):
  state=dict(stage='headers',total=7020,checked=1450,accepted=1448,
             rejected=2,current='image_001_29.rasx',last_ok='image_001_28.rasx')
  msg=gui.format_scan_progress(state,65)
  self.assertIn('1450/7020',msg)
  self.assertIn('прошли 1448',msg)
  self.assertIn('ошибки 2',msg)
  self.assertIn('image_001_29.rasx',msg)
  self.assertIn('image_001_28.rasx',msg)
  self.assertIn('65 с',msg)

 def test_loading_and_group_validation_stages(self):
  self.assertIn('Подготовка',gui.format_scan_progress({'stage':'init'},5))
  self.assertIn('найдено 39',gui.format_scan_progress({'stage':'enumerate','found':39},10))
  state={'stage':'validate','groups_done':20,'groups_total':39,'groups_ok':19}
  self.assertIn('20/39',gui.format_scan_progress(state,80))

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
