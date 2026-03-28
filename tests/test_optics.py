import unittest
import numpy as np
from unittest.mock import patch, call

from vasppy.optics import (
    matrix_eigvals,
    to_matrix,
    parse_dielectric_data,
    absorption_coefficient,
    plot_dielectric_functions,
)


class TestOptics(unittest.TestCase):
    def test_matrix_eigvals(self):
        matrix = np.array([[2, 0, 3], [0, 3, 0], [0, 0, 3]])
        expected_eigenvalues = np.array([2, 3, 3])
        np.testing.assert_array_equal(matrix_eigvals(matrix), expected_eigenvalues)

    def test_to_matrix(self):
        expected_matrix = np.array([[1, 2, 3], [2, 4, 5], [3, 5, 6]])
        np.testing.assert_array_equal(
            to_matrix(xx=1, yy=4, zz=6, xy=2, yz=5, xz=3), expected_matrix
        )

    @patch("vasppy.optics.to_matrix")
    @patch("vasppy.optics.matrix_eigvals")
    def test_parse_dielectric_data(self, mock_matrix_eigvals, mock_to_matrix):
        input_data = [[0.1, 0.2, 0.3, 0.4, 0.5, 0.6], [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]]
        matrix_a = np.array([[0.1, 0.4, 0.6], [0.4, 0.2, 0.5], [0.6, 0.5, 0.3]])
        matrix_b = np.array([[1.1, 1.4, 1.6], [1.4, 1.2, 1.5], [1.6, 1.5, 1.3]])
        mock_to_matrix.side_effect = [matrix_a, matrix_b]
        mock_matrix_eigvals.side_effect = [np.array([1, 2, 3]), np.array([4, 5, 6])]
        expected_data = np.array([[1, 2, 3], [4, 5, 6]])
        np.testing.assert_array_equal(parse_dielectric_data(input_data), expected_data)
        mock_to_matrix.assert_has_calls([call(*input_data[0]), call(*input_data[1])])
        mock_matrix_eigvals.assert_has_calls([call(matrix_a), call(matrix_b)])

    # ------------------------------------------------------------------ #
    # absorption_coefficient
    # ------------------------------------------------------------------ #

    def test_absorption_coefficient_returns_array_of_correct_length(self):
        """Test that absorption_coefficient returns an array of the same
        length as the energy list."""
        n_points = 5
        energies = list(np.linspace(0.1, 5.0, n_points))
        # Diagonal tensors: xx=yy=zz=value, off-diagonal=0
        real_tensors = [[2.0, 2.0, 2.0, 0.0, 0.0, 0.0]] * n_points
        imag_tensors = [[1.0, 1.0, 1.0, 0.0, 0.0, 0.0]] * n_points
        dielectric = [energies, real_tensors, imag_tensors]
        result = absorption_coefficient(dielectric)
        self.assertEqual(result.shape, (n_points,))

    def test_absorption_coefficient_is_zero_at_zero_energy(self):
        """Test that the absorption coefficient is zero when energy is zero."""
        energies = [0.0, 1.0]
        real_tensors = [[2.0, 2.0, 2.0, 0.0, 0.0, 0.0]] * 2
        imag_tensors = [[1.0, 1.0, 1.0, 0.0, 0.0, 0.0]] * 2
        dielectric = [energies, real_tensors, imag_tensors]
        result = absorption_coefficient(dielectric)
        self.assertAlmostEqual(result[0], 0.0)

    def test_absorption_coefficient_is_non_negative(self):
        """Test that the absorption coefficient is non-negative for physical inputs."""
        energies = list(np.linspace(0.5, 5.0, 10))
        real_tensors = [[2.0, 2.0, 2.0, 0.0, 0.0, 0.0]] * 10
        imag_tensors = [[0.5, 0.5, 0.5, 0.0, 0.0, 0.0]] * 10
        dielectric = [energies, real_tensors, imag_tensors]
        result = absorption_coefficient(dielectric)
        self.assertTrue(np.all(result >= 0.0))

    def test_absorption_coefficient_numerical_value(self):
        """Test absorption_coefficient against a hand-calculated value."""
        from math import pi, sqrt
        from scipy.constants import physical_constants, speed_of_light

        eV_to_recip_cm = 1.0 / (
            physical_constants["Planck constant in eV s"][0] * speed_of_light * 1e2
        )
        # Single-point isotropic case: eps1=2, eps2=1
        energy = 2.0
        eps1, eps2 = 2.0, 1.0
        expected = (
            2.0 * sqrt(2.0) * pi * eV_to_recip_cm * energy
            * sqrt(-eps1 + sqrt(eps1**2 + eps2**2))
        )
        energies = [energy]
        real_tensors = [[eps1, eps1, eps1, 0.0, 0.0, 0.0]]
        imag_tensors = [[eps2, eps2, eps2, 0.0, 0.0, 0.0]]
        dielectric = [energies, real_tensors, imag_tensors]
        result = absorption_coefficient(dielectric)
        self.assertAlmostEqual(result[0], expected, places=6)

    # ------------------------------------------------------------------ #
    # plot_dielectric_functions
    # ------------------------------------------------------------------ #

    def test_plot_dielectric_functions_creates_new_figure_when_ax_is_none(self):
        """Test that a new Figure is returned when no axes are supplied."""
        energies = list(np.linspace(0.1, 8.0, 10))
        real_tensors = [[2.0, 2.0, 2.0, 0.0, 0.0, 0.0]] * 10
        imag_tensors = [[0.5, 0.5, 0.5, 0.0, 0.0, 0.0]] * 10
        dielectric = [energies, real_tensors, imag_tensors]

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = plot_dielectric_functions(dielectric)
        self.assertIsNotNone(fig)
        plt.close("all")

    def test_plot_dielectric_functions_returns_none_when_ax_is_supplied(self):
        """Test that None is returned when an existing Axes is supplied."""
        energies = list(np.linspace(0.1, 8.0, 10))
        real_tensors = [[2.0, 2.0, 2.0, 0.0, 0.0, 0.0]] * 10
        imag_tensors = [[0.5, 0.5, 0.5, 0.0, 0.0, 0.0]] * 10
        dielectric = [energies, real_tensors, imag_tensors]

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        _, ax = plt.subplots()
        result = plot_dielectric_functions(dielectric, ax=ax)
        self.assertIsNone(result)
        plt.close("all")

    def test_plot_dielectric_functions_plots_two_lines(self):
        """Test that two lines are drawn (real and imaginary)."""
        energies = list(np.linspace(0.1, 8.0, 10))
        real_tensors = [[2.0, 2.0, 2.0, 0.0, 0.0, 0.0]] * 10
        imag_tensors = [[0.5, 0.5, 0.5, 0.0, 0.0, 0.0]] * 10
        dielectric = [energies, real_tensors, imag_tensors]

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        _, ax = plt.subplots()
        plot_dielectric_functions(dielectric, ax=ax)
        self.assertEqual(len(ax.lines), 2)
        plt.close("all")


if __name__ == "__main__":
    unittest.main()
