import inspect
import unittest
import lmo_pf2d_gui as gui


class ProgressMessageTests(unittest.TestCase):
    def test_file_progress_includes_stage_counts_last_checked_and_errors(self):
        state=dict(stage="headers",total=7020,checked=1450,accepted=1448,
                   rejected=2,current="image_001_29.rasx",last_ok="image_001_28.rasx")
        msg=gui.format_scan_progress(state,65)
        self.assertIn("1450/7020",msg)
        self.assertIn("accepted 1448",msg)
        self.assertIn("errors 2",msg)
        self.assertIn("image_001_29.rasx",msg)
        self.assertIn("image_001_28.rasx",msg)
        self.assertIn("65 s",msg)

    def test_loading_and_validation_stages(self):
        self.assertIn("Preparing",gui.format_scan_progress({"stage":"init"},5))
        self.assertIn("39 found",gui.format_scan_progress({"stage":"enumerate","found":39},10))
        state={"stage":"validate","groups_done":20,"groups_total":39,"groups_ok":19}
        self.assertIn("20/39",gui.format_scan_progress(state,80))


class AdaptiveGUIWiring(unittest.TestCase):
    def test_new_one_pass_version(self):
        self.assertGreaterEqual(tuple(int(x) for x in gui.VERSION.split(".")),(1,3,0))

    def test_adaptive_dispatch(self):
        self.assertIn("adaptive_pf2d.py",inspect.getsource(gui.run_analyze_mode))

    def test_metadata_check_is_not_required_before_analysis(self):
        source=inspect.getsource(gui.make_gui)
        self.assertIn("Analyze & Generate Reports",source)
        self.assertIn("Rebuild PDFs from Saved Maps",source)
        self.assertNotIn('self.scan_button=',source)
        self.assertNotIn("scan_groups(",source)

    def test_gui_strings_are_english(self):
        source=inspect.getsource(gui.make_gui)
        self.assertFalse(any(1024 <= ord(c) <= 1279 for c in source))
        self.assertIn("One-pass workflow",source)


if __name__=="__main__":
    unittest.main()
