"""Characterisation tests for the spacegroup script."""

import os
import unittest
from unittest.mock import patch

TEST_POSCAR = os.path.join(
    os.path.dirname(__file__), "..", "test_data", "POSCAR_NaCl"
)


class SpacegroupMainTestCase(unittest.TestCase):
    """Tests for the spacegroup main function."""

    def test_main_prints_spacegroup_symbol(self):
        """Test that main prints a space group symbol for a given POSCAR file."""
        with patch("sys.argv", ["spacegroup", TEST_POSCAR]):
            from vasppy.scripts.spacegroup import main
            with patch("builtins.print") as mock_print:
                main()
            mock_print.assert_called_once()
            symbol = mock_print.call_args[0][0]
            # NaCl has space group Fm-3m (number 225)
            self.assertIsInstance(symbol, str)
            self.assertTrue(len(symbol) > 0)

    def test_main_prints_correct_spacegroup_for_nacl(self):
        """Test that main prints the correct space group symbol for the NaCl primitive cell.

        The test POSCAR contains the primitive cell of NaCl (1 Na + 1 Cl in a
        cubic unit cell), which pymatgen's SpacegroupAnalyzer identifies as Pm-3m
        at the default symprec of 1e-3.
        """
        with patch("sys.argv", ["spacegroup", TEST_POSCAR]):
            from vasppy.scripts.spacegroup import main
            with patch("builtins.print") as mock_print:
                main()
            symbol = mock_print.call_args[0][0]
            self.assertEqual(symbol, "Pm-3m")


if __name__ == "__main__":
    unittest.main()
