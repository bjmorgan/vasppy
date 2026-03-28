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
                    names, values = potcar_spec(mock_potcar_filename)
                    mock_open.assert_called_with(mock_potcar_filename, "r")
                    mock_md5sum.assert_has_calls(
                        [
                            call("foo\nEnd of Dataset\n"),
                            call("bar\nEnd of Dataset\n"),
                            call("sds\nEnd of Dataset\n"),
                        ]
                    )
        self.assertEqual(names, ["A", "E", "D"])
        self.assertEqual(values, ["PBE", "PBE_54", "PBE_52"])

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
                    names, hashes = potcar_spec(
                        mock_potcar_filename, return_hashes=True,
                    )
                    mock_open.assert_called_with(mock_potcar_filename, "r")
                    mock_md5sum.assert_has_calls(
                        [
                            call("foo\nEnd of Dataset\n"),
                            call("bar\nEnd of Dataset\n"),
                            call("sds\nEnd of Dataset\n"),
                        ]
                    )
        self.assertEqual(names, ["A", "E", "D"])
        self.assertEqual(hashes, ["12", "56", "23"])

    def test_potcar_spec_preserves_duplicate_species(self):
        """Duplicate species must appear in both returned lists."""
        mock_potcar_string_dupes = (
            "foo\nEnd of Dataset\nbar\nEnd of Dataset\n"
        )
        mock_potcar_data_dupes = {
            "PBE": {"O": "aa"},
            "PBE_52": {},
            "PBE_54": {},
            "PBE_54r": {},
            "LDA": {},
            "LDA_52": {},
            "LDA_54": {},
            "LDA_54r": {},
            "GGA": {},
            "USPP_GGA": {},
            "USPP_LDA": {},
        }
        with patch(
            "builtins.open",
            return_value=io.StringIO(mock_potcar_string_dupes),
        ):
            with patch(
                "vasppy.summary.md5sum", side_effect=("aa", "aa")
            ):
                with patch.dict(
                    "vasppy.data.potcar_data.potcar_md5sum_data",
                    mock_potcar_data_dupes,
                    clear=True,
                ):
                    names, values = potcar_spec("POTCAR")
        self.assertEqual(names, ["O", "O"])
        self.assertEqual(values, ["PBE", "PBE"])

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


class SummaryPrintMethodsTestCase(unittest.TestCase):
    """Tests for the remaining Summary print_* methods."""

    @patch("vasppy.summary.VASPMeta")
    @patch("vasppy.summary.Summary.parse_vasprun")
    def setUp(self, mock_parse_vasprun, MockVASPMeta):
        MockVASPMeta.from_file = Mock(return_value="foo")
        self.summary = Summary()
        self.summary.vasprun = Mock(spec=Vasprun)
        self.summary.meta = Mock(spec=VASPMeta)
        self.summary.meta.notes = None
        self.summary.directory = "."
        self.summary.vasprun_filename = "vasprun.xml"

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_status(self, mock_stdout):
        self.summary.meta.status = "converged"
        self.summary.print_status()
        self.assertEqual(mock_stdout.getvalue(), "status: converged\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_lreal(self, mock_stdout):
        self.summary.vasprun.parameters = {"LREAL": False}
        self.summary.print_lreal()
        self.assertEqual(mock_stdout.getvalue(), "lreal: False\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_description(self, mock_stdout):
        self.summary.meta.description = "  a description  "
        self.summary.print_description()
        self.assertEqual(mock_stdout.getvalue(), "description: a description\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_stoichiometry(self, mock_stdout):
        self.summary.vasprun.final_structure = Mock()
        self.summary.vasprun.final_structure.composition.get_el_amt_dict.return_value = {
            "Fe": 2.0, "O": 3.0
        }
        self.summary.print_stoichiometry()
        output = mock_stdout.getvalue()
        self.assertIn("stoichiometry:", output)
        self.assertIn("Fe: 2", output)
        self.assertIn("O: 3", output)

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_potcar(self, mock_stdout):
        self.summary.vasprun.final_structure = Mock()
        self.summary.vasprun.final_structure.composition.get_el_amt_dict.return_value = {
            "Fe": 2.0
        }
        self.summary.vasprun.potcar_symbols = ["PAW_PBE Fe_pv"]
        self.summary.print_potcar()
        output = mock_stdout.getvalue()
        self.assertIn("potcar:", output)
        self.assertIn("Fe: PAW_PBE Fe_pv", output)

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_energy_no_type(self, mock_stdout):
        self.summary.meta.type = None
        self.summary.vasprun.final_energy = -10.5
        self.summary.print_energy()
        self.assertEqual(mock_stdout.getvalue(), "energy: -10.5\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_energy_unsupported_type_raises(self, mock_stdout):
        self.summary.meta.type = "foo"
        with self.assertRaises(ValueError):
            self.summary.print_energy()

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_functional(self, mock_stdout):
        self.summary.vasprun.run_type = "GGA"
        self.summary.print_functional()
        self.assertEqual(mock_stdout.getvalue(), "functional: GGA\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_ibrion(self, mock_stdout):
        self.summary.vasprun.incar = {"IBRION": 2}
        self.summary.print_ibrion()
        self.assertEqual(mock_stdout.getvalue(), "ibrion: 2\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_ediffg(self, mock_stdout):
        self.summary.vasprun.incar = {"EDIFFG": -0.01}
        self.summary.print_ediffg()
        self.assertEqual(mock_stdout.getvalue(), "ediffg: -0.01\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_encut_with_encut(self, mock_stdout):
        self.summary.vasprun.incar = {"ENCUT": 520}
        self.summary.print_encut()
        self.assertEqual(mock_stdout.getvalue(), "encut: 520\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_encut_with_enmax(self, mock_stdout):
        self.summary.vasprun.incar = {"ENMAX": 500}
        self.summary.print_encut()
        self.assertEqual(mock_stdout.getvalue(), "encut: 500\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_encut_no_key(self, mock_stdout):
        self.summary.vasprun.incar = {}
        self.summary.print_encut()
        self.assertEqual(mock_stdout.getvalue(), "")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_directory(self, mock_stdout):
        self.summary.directory = "/some/path"
        self.summary.print_directory()
        self.assertEqual(mock_stdout.getvalue(), "directory: /some/path\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_nelect(self, mock_stdout):
        self.summary.vasprun.parameters = {"NELECT": 64.0}
        self.summary.print_nelect()
        self.assertEqual(mock_stdout.getvalue(), "nelect: 64.0\n")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_plus_u_no_ldauu(self, mock_stdout):
        self.summary.vasprun.incar = {}
        self.summary.print_plus_u()
        self.assertEqual(mock_stdout.getvalue(), "")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_plus_u_all_zero(self, mock_stdout):
        self.summary.vasprun.incar = {
            "LDAUU": [0, 0],
            "LDAUJ": [0, 0],
            "LDAUL": [2, 2],
        }
        self.summary.vasprun.final_structure = Mock()
        self.summary.vasprun.final_structure.composition.get_el_amt_dict.return_value = {
            "Fe": 2.0, "O": 3.0
        }
        self.summary.print_plus_u()
        self.assertEqual(mock_stdout.getvalue(), "")

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_plus_u_with_nonzero(self, mock_stdout):
        self.summary.vasprun.incar = {
            "LDAUU": [4, 0],
            "LDAUJ": [1, 0],
            "LDAUL": [2, 1],
        }
        self.summary.vasprun.final_structure = Mock()
        self.summary.vasprun.final_structure.composition.get_el_amt_dict.return_value = {
            "Fe": 2.0, "O": 3.0
        }
        self.summary.print_plus_u()
        output = mock_stdout.getvalue()
        self.assertIn("ldau:", output)
        self.assertIn("Fe: d 4 1", output)

    @patch("sys.stdout", new_callable=StringIO)
    @patch("vasppy.summary.file_md5", return_value="abc123")
    def test_print_vasprun_md5(self, mock_md5, mock_stdout):
        self.summary.print_vasprun_md5()
        self.assertIn("abc123", mock_stdout.getvalue())

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_file_tracking_no_track(self, mock_stdout):
        self.summary.meta.track = None
        self.summary.print_file_tracking()
        self.assertEqual(mock_stdout.getvalue(), "")

    @patch("sys.stdout", new_callable=StringIO)
    def test_output_skips_most_when_vasprun_is_none(self, mock_stdout):
        self.summary.vasprun = None
        self.summary.meta.title = "test"
        self.summary.meta.type = None
        self.summary.meta.status = "done"
        self.summary.output(["title", "type", "status", "energy"])
        output = mock_stdout.getvalue()
        self.assertIn("title: test", output)
        self.assertNotIn("energy", output)

    @patch("sys.stdout", new_callable=StringIO)
    def test_print_kpoints(self, mock_stdout):
        from unittest.mock import MagicMock
        self.summary.vasprun.kpoints = MagicMock()
        self.summary.vasprun.kpoints.style = "Gamma"
        self.summary.vasprun.kpoints.kpts = [[4, 4, 4]]
        self.summary.print_kpoints()
        output = mock_stdout.getvalue()
        self.assertIn("scheme: Gamma", output)
        self.assertIn("grid: 4 4 4", output)


if __name__ == "__main__":
    unittest.main()
