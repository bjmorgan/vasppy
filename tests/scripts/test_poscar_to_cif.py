"""Characterisation tests for the poscar_to_cif script."""

import os
import unittest
from unittest.mock import patch

TEST_POSCAR = os.path.join(
    os.path.dirname(__file__), "..", "test_data", "POSCAR_NaCl"
)


class PoscarToCifMainTestCase(unittest.TestCase):
    """Tests for the poscar_to_cif main function."""

    def test_main_produces_cif_output_without_symprec(self):
        """Test that main prints CIF output for a given POSCAR file."""
        with patch("sys.argv", ["poscar_to_cif", TEST_POSCAR]):
            from vasppy.scripts.poscar_to_cif import main
            with patch("builtins.print") as mock_print:
                main()
            mock_print.assert_called_once()
            output = str(mock_print.call_args[0][0])
            self.assertIn("_cell_length_a", output)

    def test_main_produces_cif_output_with_symprec(self):
        """Test that main prints symmetrised CIF output when symprec is given."""
        with patch("sys.argv", ["poscar_to_cif", TEST_POSCAR, "--symprec", "0.01"]):
            from vasppy.scripts.poscar_to_cif import main
            with patch("builtins.print") as mock_print:
                main()
            mock_print.assert_called_once()
            output = str(mock_print.call_args[0][0])
            self.assertIn("_cell_length_a", output)


if __name__ == "__main__":
    unittest.main()
