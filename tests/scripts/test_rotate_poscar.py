"""Characterisation tests for the rotate_poscar script."""

import os
import unittest
from unittest.mock import patch

import numpy as np

TEST_POSCAR = os.path.join(
    os.path.dirname(__file__), "..", "test_data", "POSCAR_NaCl"
)


class RotatePoscarMainTestCase(unittest.TestCase):
    """Integration tests for the rotate_poscar main function."""

    def test_main_produces_poscar_output(self):
        """Test that main prints POSCAR output for a given POSCAR file."""
        with patch("sys.argv", ["rotate_poscar", TEST_POSCAR, "-a", "0", "0", "1", "-d", "90"]):
            from vasppy.scripts.rotate_poscar import main
            with patch("builtins.print") as mock_print:
                main()
            mock_print.assert_called_once()
            output = str(mock_print.call_args[0][0])
            self.assertIn("Na", output)
            self.assertIn("Cl", output)

    def test_90_degree_rotation_about_z_swaps_x_and_y(self):
        """Test that a 90-degree rotation about z maps the a-vector to the b-direction.

        The NaCl test POSCAR has a cubic cell with vectors:
            a = [2.82, 0, 0]
            b = [0, 2.82, 0]
            c = [0, 0, 2.82]

        A 90-degree counterclockwise rotation about z should map:
            a -> [0, 2.82, 0]   (former b-direction)
            b -> [-2.82, 0, 0]  (former -a-direction)
            c -> [0, 0, 2.82]   (unchanged)
        """
        with patch("sys.argv", ["rotate_poscar", TEST_POSCAR, "-a", "0", "0", "1", "-d", "90"]):
            from vasppy.scripts.rotate_poscar import main
            with patch("builtins.print") as mock_print:
                main()
            output = str(mock_print.call_args[0][0])

        # Extract lattice vectors from the POSCAR output.
        # Standard POSCAR: line 0=comment, 1=scale, 2-4=lattice vectors.
        lines = output.strip().splitlines()
        scale = float(lines[1].strip())
        a_vec = np.array([float(x) for x in lines[2].split()]) * scale
        b_vec = np.array([float(x) for x in lines[3].split()]) * scale
        c_vec = np.array([float(x) for x in lines[4].split()]) * scale

        length = 2.82
        np.testing.assert_allclose(a_vec, [0.0, length, 0.0], atol=1e-5)
        np.testing.assert_allclose(b_vec, [-length, 0.0, 0.0], atol=1e-5)
        np.testing.assert_allclose(c_vec, [0.0, 0.0, length], atol=1e-5)

    def test_zero_degree_rotation_leaves_cell_unchanged(self):
        """Test that a 0-degree rotation leaves the lattice vectors unchanged."""
        with patch("sys.argv", ["rotate_poscar", TEST_POSCAR, "-a", "0", "0", "1", "-d", "0"]):
            from vasppy.scripts.rotate_poscar import main
            with patch("builtins.print") as mock_print:
                main()
            output = str(mock_print.call_args[0][0])

        lines = output.strip().splitlines()
        scale = float(lines[1].strip())
        a_vec = np.array([float(x) for x in lines[2].split()]) * scale
        b_vec = np.array([float(x) for x in lines[3].split()]) * scale
        c_vec = np.array([float(x) for x in lines[4].split()]) * scale

        length = 2.82
        np.testing.assert_allclose(a_vec, [length, 0.0, 0.0], atol=1e-5)
        np.testing.assert_allclose(b_vec, [0.0, length, 0.0], atol=1e-5)
        np.testing.assert_allclose(c_vec, [0.0, 0.0, length], atol=1e-5)


if __name__ == "__main__":
    unittest.main()
