import unittest
import io
import inspect
from unittest.mock import Mock, patch, call
from io import StringIO
import os

from vasppy.summary import (
    Summary,
    md5sum,
    potcar_spec,
    find_vasp_calculations,
    load_vasp_summary,
)
from vasppy.vaspmeta import VASPMeta

from pymatgen.io.vasp.outputs import Vasprun

test_data_dir = "test_data"

mock_potcar_string = """foo
End of Dataset
bar
End of Dataset
sds
End of Dataset
"""

mock_potcar_data = {
    "PBE": {"A": "12", "B": "34"},
    "PBE_52": {"C": "01", "D": "23"},
    "PBE_54": {"E": "56", "F": "78"},
    "LDA": {"G": "89"},
    "LDA_52": {"H": "101"},
    "LDA_54": {"I": "202"},
    "GGA": {"J": "303"},
    "USPP_GGA": {"K": "404"},
    "USPP_LDA": {"L": "505"},
    "PBE_54r": {"M": "123"},
    "LDA_54r": {"N": "456"},
}


class SummaryInitTestCase(unittest.TestCase):
    @patch("vasppy.summary.VASPMeta")
    @patch("vasppy.summary.Summary.parse_vasprun")
    def test_summary_is_initialised(self, mock_parse_vasprun, MockVASPMeta):
        MockVASPMeta.from_file = Mock(return_value="foo")
        summary = Summary()
        self.assertEqual(mock_parse_vasprun.call_count, 1)
        expected_print_methods = [
            "title",
            "type",
            "status",
            "stoichiometry",
            "potcar",
            "eatom",
            "energy",
            "k-points",
            "functional",
            "encut",
            "plus_u",
            "ediffg",
            "ibrion",
            "converged",
            "version",
            "md5",
            "directory",
            "lreal",
            "vbm",
            "cbm",
        ]
        for key in expected_print_methods:
            self.assertTrue(key in summary.print_methods)
            self.assertTrue(inspect.ismethod(summary.print_methods[key]))

    @patch("vasppy.summary.VASPMeta")
    @patch("vasppy.summary.Summary.parse_vasprun")
    def test_summary_init_raises_filenotfounderror_if_file_is_not_found(
        self, mock_parse_vasprun, MockVASPMeta
    ):
        MockVASPMeta.from_file = Mock(side_effect=FileNotFoundError)
        with self.assertRaises(FileNotFoundError):
            Summary()


class SummaryTestCase(unittest.TestCase):
    @patch("vasppy.summary.VASPMeta")
    @patch("vasppy.summary.Summary.parse_vasprun")
    def setUp(self, mock_parse_vaspun, MockVASPMeta):
        MockVASPMeta.from_file = Mock(return_value="foo")
        self.summary = Summary()
        self.summary.vasprun = Mock(spec=Vasprun)
        self.summary.meta = Mock(spec=VASPMeta)
        self.summary.meta.notes = None

    def test_functional(self):
        self.summary.vasprun.run_type = "foo"
        self.assertEqual(self.summary.functional, "foo")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_cbm(self, mock_stdout):
        summary = self.summary
        summary.vasprun.eigenvalue_band_properties = ["null", "CBM", "VBM"]
        summary.print_cbm()
        self.assertEqual(mock_stdout.getvalue(), "cbm: CBM\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_vbm(self, mock_stdout):
        self.summary.vasprun.eigenvalue_band_properties = ["null", "CBM", "VBM"]
        self.summary.print_vbm()
        self.assertEqual(mock_stdout.getvalue(), "vbm: VBM\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_converged(self, mock_stdout):
        self.summary.vasprun.converged = "conv"
        self.summary.print_converged()
        self.assertEqual(mock_stdout.getvalue(), "converged: conv\n")

    def test_potcars_are_pbe_if_true(self):
        self.summary.vasprun.potcar_symbols = [
            "PAW_PBE Fe_pv 06Sep2000",
            "PAW_PBE O 08Apr2002",
        ]
        self.assertTrue(self.summary.potcars_are_pbe())

    def test_potcars_are_pbe_if_false(self):
        self.summary.vasprun.potcar_symbols = ["foo", "PAW_PBE O 08Apr2002"]
        self.assertFalse(self.summary.potcars_are_pbe())

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_type(self, mock_stdout):
        self.summary.meta.type = "TYPE"
        self.summary.print_type()
        self.assertEqual(mock_stdout.getvalue(), "type: TYPE\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_type_if_type_is_not_set(self, mock_stdout):
        self.summary.meta.type = None
        self.summary.print_type()
        self.assertEqual(mock_stdout.getvalue(), "")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_title(self, mock_stdout):
        self.summary.meta.title = "TITLE"
        self.summary.print_title()
        self.assertEqual(mock_stdout.getvalue(), "title: TITLE\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_notes(self, mock_stdout):
        self.summary.meta.notes = "NOTES"
        self.summary.print_notes()
        self.assertEqual(mock_stdout.getvalue(), "notes: NOTES\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_notes_handles_empty_notes_attribute(self, mock_stdout):
        self.summary.print_notes()
        self.assertEqual(mock_stdout.getvalue(), "notes: ~\n")


class SummaryHelperFunctionsTestCase(unittest.TestCase):
    def test_md5sum(self):
        self.assertEqual(md5sum("hello\n"), "b1946ac92492d2347c6235b4d2611184")

    def test_potcar_spec(self):
        mock_potcar_filename = "POTCAR"
        md5sum_return_values = ("12", "56", "23")
        with patch(
            "builtins.open", return_value=io.StringIO(mock_potcar_string)
        ) as mock_open:
            with patch(
                "vasppy.summary.md5sum", side_effect=md5sum_return_values
            ) as mock_md5sum:
                with patch.dict(
                    "vasppy.data.potcar_data.potcar_md5sum_data",
                    mock_potcar_data,
                    clear=True,
                ):
                    p_spec = potcar_spec(mock_potcar_filename)
                    mock_open.assert_called_with(mock_potcar_filename, "r")
                    mock_md5sum.assert_has_calls(
                        [
                            call("foo\nEnd of Dataset\n"),
                            call("bar\nEnd of Dataset\n"),
                            call("sds\nEnd of Dataset\n"),
                        ]
                    )
        self.assertEqual(p_spec, {"A": "PBE", "E": "PBE_54", "D": "PBE_52"})

    def test_potcar_spec_returns_hashes(self):
        mock_potcar_filename = "POTCAR"
        md5sum_return_values = ("12", "56", "23")
        with patch(
            "builtins.open", return_value=io.StringIO(mock_potcar_string)
        ) as mock_open:
            with patch(
                "vasppy.summary.md5sum", side_effect=md5sum_return_values
            ) as mock_md5sum:
                with patch.dict(
                    "vasppy.data.potcar_data.potcar_md5sum_data",
                    mock_potcar_data,
                    clear=True,
                ):
                    p_spec = potcar_spec(mock_potcar_filename, return_hashes=True)
                    mock_open.assert_called_with(mock_potcar_filename, "r")
                    mock_md5sum.assert_has_calls(
                        [
                            call("foo\nEnd of Dataset\n"),
                            call("bar\nEnd of Dataset\n"),
                            call("sds\nEnd of Dataset\n"),
                        ]
                    )
        self.assertEqual(p_spec, {"A": "12", "E": "56", "D": "23"})

    def test_potcar_spec_raises_valueerror_if_md5sum_not_matched(self):
        mock_potcar_filename = "POTCAR"
        md5sum_return_values = ("12", "56", "90")
        with patch("builtins.open", return_value=io.StringIO(mock_potcar_string)):
            with patch("vasppy.summary.md5sum", side_effect=md5sum_return_values):
                with patch.dict(
                    "vasppy.data.potcar_data.potcar_md5sum_data",
                    mock_potcar_data,
                    clear=True,
                ):
                    with self.assertRaises(ValueError):
                        potcar_spec(mock_potcar_filename)

    def test_find_vasp_calculations(self):
        mock_glob_output = ["dir_A/vasprun.xml", "dir_B/dir_C/vasprun.xml"]
        with patch("glob.iglob", side_effect=[mock_glob_output, []]):
            v = find_vasp_calculations()
        self.assertEqual(v, ["./dir_A/", "./dir_B/dir_C/"])

    def test_load_vasp_summary(self):
        vasp_summary_test_filename = os.path.join(
            os.path.dirname(__file__), test_data_dir, "vasp_summary_test.yaml"
        )
        expected_dict = {
            "foo": {"title": "foo", "data": "foo_data"},
            "bar": {"title": "bar", "data": "bar_data"},
        }
        vasp_summary = load_vasp_summary(vasp_summary_test_filename)
        self.assertEqual(vasp_summary, expected_dict)


if __name__ == "__main__":
    unittest.main()
