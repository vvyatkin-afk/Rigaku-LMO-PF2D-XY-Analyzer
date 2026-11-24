"""Regression checks for selecting nested folders and both Windows worker modes."""
import pathlib
import subprocess
import sys
import tempfile
import unittest

from lmo_pf2d_gui import VERSION, choose_existing_batch, choose_rasx_folder, ensure_output_parent


class FolderAndCliTests(unittest.TestCase):
    def test_select_only_nested_raw_folder(self):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp)
            raw=root/"batch_29"/"raw_scan"
            raw.mkdir(parents=True)
            (raw/"001_IMAGE_1_38.000_001_01.RASX").write_bytes(b"stub")
            other=root/"LMO_PF2D_BATCH_old"/"point_001"
            other.mkdir(parents=True)
            self.assertEqual(choose_rasx_folder(root),raw)
            self.assertEqual(choose_rasx_folder(raw),raw)

    def test_multi_sample_parent_folder(self):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp)
            for name in ("run1","run2"):
                d=root/name
                d.mkdir()
                (d/"frame.rasx").touch()
            self.assertEqual(choose_rasx_folder(root),root)

    def test_select_unique_previous_batch_subfolder(self):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp)
            batch=root/"LMO_PF2D_BATCH_20261123"
            point=batch/"point_001"
            point.mkdir(parents=True)
            (batch/"inspection.json").write_text('{"groups":[]}',encoding="utf-8")
            (point/"map.npz").write_bytes(b"stub")
            self.assertEqual(choose_existing_batch(root),batch)
            self.assertEqual(choose_existing_batch(batch),batch)

    def test_autocreate_output_with_missing_parents(self):
        with tempfile.TemporaryDirectory() as temp:
            p=pathlib.Path(temp)/"reports"/"new"/"nested"
            self.assertFalse(p.exists())
            self.assertEqual(ensure_output_parent(str(p)),p)
            self.assertTrue(p.is_dir())
            self.assertEqual(ensure_output_parent(str(p)),p)
            with self.assertRaises(ValueError):
                ensure_output_parent("")

    def test_both_modes_accepted_by_real_commandline_parser(self):
        app=pathlib.Path(__file__).resolve().parents[1]/"lmo_pf2d_gui.py"
        child=subprocess.run([sys.executable,str(app),"--analyze","--help"],
                             capture_output=True,text=True,errors="replace",timeout=40)
        self.assertEqual(child.returncode,0,child.stderr)
        self.assertIn("--reports-only",child.stdout)
        self.assertIn("--source",child.stdout)
        self.assertIn("--group-key",child.stdout)

    def test_version_bumped(self):
        self.assertEqual(VERSION,"1.5.1")


if __name__=="__main__":
    unittest.main()
