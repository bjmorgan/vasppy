import unittest
import os
from vasppy import procar
import numpy as np
from unittest.mock import patch
import warnings
from copy import deepcopy

test_data_dir = "test_data"
test_procar_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "PROCAR_test"
)
test_procar_spin_polarised_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "PROCAR_spin_polarised_test"
)


class KPointTestCase(unittest.TestCase):
    """Test for procar.KPoint class"""

    def setUp(self):
        index = 1
        frac_coords = np.array([0.1, 0.2, 0.3])
        weight = 0.1
        self.k_point = procar.KPoint(
            index=index, frac_coords=frac_coords, weight=weight
        )

    def test_kpoint_is_initialised(self):
        """Test KPoint object is initialised"""
        index = 1
        frac_coords = np.array([0.1, 0.2, 0.3])
        weight = 0.1
        k_point = procar.KPoint(index=index, frac_coords=frac_coords, weight=weight)
        self.assertEqual(k_point.index, index)
        np.testing.assert_equal(k_point.frac_coords, frac_coords)
        self.assertEqual(k_point.index, index)

    def test_cart_coords(self):
        reciprocal_lattice = np.array(
            [[10.0, 10.0, 0.0], [10.0, 0.0, 10.0], [0.0, 10.0, 10.0]]
        )
        np.testing.assert_equal(
            self.k_point.cart_coords(reciprocal_lattice=reciprocal_lattice),
            np.dot(self.k_point.frac_coords, reciprocal_lattice),
        )

    def test___eq___equal_kpoints(self):
        other_k_point = deepcopy(self.k_point)
        self.assertTrue(self.k_point == other_k_point)

    def test___eq___unequal_index(self):
        other_k_point = deepcopy(self.k_point)
        other_k_point.index = 99
        self.assertFalse(self.k_point == other_k_point)

    def test___eq___unequal_frac_coords(self):
        other_k_point = deepcopy(self.k_point)
        other_k_point.frac_coords = np.array([0.9, 0.8, 0.7])
        self.assertFalse(self.k_point == other_k_point)

    def test___eq___unequal_weight(self):
        other_k_point = deepcopy(self.k_point)
        other_k_point.weight = 0.99
        self.assertFalse(self.k_point == other_k_point)


def _load_procar(filename=test_procar_filename, **kwargs):
    """Helper to load a Procar from file, suppressing warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return procar.Procar.from_file(filename, **kwargs)


class ProcarTestCase(unittest.TestCase):
    """Test for Procar class"""

    def test_procar_is_read_from_file(self):
        """Checking that `PROCAR_test` is read"""
        pcar = _load_procar()
        self.assertEqual(pcar.spin_channels, 4)
        self.assertEqual(pcar.number_of_ions, 22)
        self.assertEqual(pcar.number_of_bands, 4)
        self.assertEqual(pcar.number_of_k_points, 2)

    def test_procar_from_file_correctly_parses_bands(self):
        pcar = _load_procar()
        np.testing.assert_equal(
            [b.index for b in pcar._bands], [1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0]
        )
        np.testing.assert_equal(
            [b.energy for b in pcar._bands],
            [
                -13.17934476,
                -13.17934476,
                -13.16936722,
                -13.16936722,
                -13.1849117,
                -13.1849117,
                -13.16621473,
                -13.16621472,
            ],
        )
        np.testing.assert_equal(
            [b.occupancy for b in pcar._bands],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -0.03191968],
        )

    def test_spin_polarised_procar_is_read_from_file(self):
        """Checking that `PROCAR_spin_polarised_test` is read"""
        pcar = _load_procar(test_procar_spin_polarised_filename)
        self.assertEqual(pcar.spin_channels, 2)
        self.assertEqual(pcar.number_of_ions, 25)
        self.assertEqual(pcar.number_of_bands, 112)
        self.assertEqual(pcar.number_of_k_points, 8)

    def test___add___(self):
        pcar1 = _load_procar()
        pcar2 = _load_procar()
        combined_pcar = pcar1 + pcar2
        self.assertEqual(combined_pcar.spin_channels, 4)
        self.assertEqual(combined_pcar.number_of_ions, 22)
        self.assertEqual(combined_pcar.number_of_bands, 4)
        self.assertEqual(combined_pcar.number_of_k_points, 4)
        expected_bands = np.ravel(np.concatenate([pcar1.bands, pcar2.bands], axis=1))
        np.testing.assert_equal(combined_pcar._bands, expected_bands)
        for k1, k2 in zip(
            combined_pcar.k_points, pcar1.k_points + pcar2.k_points
        ):
            np.testing.assert_equal(k1.frac_coords, k2.frac_coords)
            self.assertEqual(k1.weight, k2.weight)
        self.assertEqual([k.index for k in combined_pcar.k_points], [1, 2, 3, 4])

    def test___add___spin_polarised_procars(self):
        pcar1 = _load_procar(test_procar_spin_polarised_filename)
        pcar2 = _load_procar(test_procar_spin_polarised_filename)
        combined_pcar = pcar1 + pcar2
        self.assertEqual(combined_pcar.spin_channels, 2)
        self.assertEqual(combined_pcar.number_of_ions, 25)
        self.assertEqual(combined_pcar.number_of_bands, 112)
        self.assertEqual(combined_pcar.number_of_k_points, 16)


class ParserTestCase(unittest.TestCase):
    """Tests for VASP PROCAR parsers"""

    def test_k_points_are_parsed(self):
        """Checking that k-points are parsed from PROCAR format strings"""
        procar_string = " k-point    1 :    0.50000000 0.25000000 0.75000000     weight = 0.00806452\nk-point    2 :    0.50000000 0.25735294 0.74264706     weight = 0.00806452"
        k_points = procar.k_point_parser(procar_string)
        np.testing.assert_array_equal(
            [k.frac_coords for k in k_points],
            [
                [0.50000000, 0.25000000, 0.75000000],
                [0.50000000, 0.25735294, 0.74264706],
            ],
        )
        self.assertEqual([k.weight for k in k_points], [0.00806452, 0.00806452])

    def test_negative_k_points_are_parsed(self):
        """Checking that negative k-points are parsed from PROCAR format strings"""
        procar_string = "  k-point  119 :   -0.01282051 0.00000000 0.00000000     weight = 0.00500000\n k-point  122 :    0.00000000-0.01282051 0.01282051     weight = 0.00500000\n k-point    1 :   -0.50000000 0.00000000-0.50000000     weight = 0.00500000"
        k_points = procar.k_point_parser(procar_string)
        np.testing.assert_array_equal(
            [k.frac_coords for k in k_points],
            [
                [-0.01282051, 0.00000000, 0.00000000],
                [0.00000000, -0.01282051, 0.01282051],
                [-0.50000000, 0.00000000, -0.50000000],
            ],
        )
        self.assertEqual([k.weight for k in k_points], [0.005, 0.005, 0.005])

    def test_k_point_parser_extra_space_after_weight_equals(self):
        """Regression test for issue #11: extra space after 'weight =' should be parsed."""
        procar_string = " k-point    1 :    0.50000000 0.25000000 0.75000000     weight =  0.00806452"
        k_points = procar.k_point_parser(procar_string)
        self.assertEqual(len(k_points), 1)
        np.testing.assert_array_equal(k_points[0].frac_coords, [0.5, 0.25, 0.75])
        self.assertAlmostEqual(k_points[0].weight, 0.00806452)

    def test_get_numbers_from_string(self):
        """Checking function for extracting numbers from a string"""
        self.assertEqual(
            procar.get_numbers_from_string("asd834asd2.11 -23as"), [834.0, 2.11, -23.0]
        )

    def test_projections_are_parsed(self):
        """Checking that projections are parsed from PROCAR format strings"""
        procar_string = "ion      s      p      d    tot\n  1  0.006  0.000  0.000  0.006\n  2  0.009  0.000  0.000  0.009\ntot  0.835  0.021  0.012  0.868\nion      s      p      d    tot\n  1  0.006  0.000  0.000  0.006\n  2  0.009  0.000  0.000  0.009\ntot  0.835  0.021  0.012  0.868\n"
        self.assertEqual(
            procar.projections_parser(procar_string).tolist(),
            np.array(
                [
                    [
                        1.0,
                        0.006,
                        0.0,
                        0.0,
                        0.006,
                        2.0,
                        0.009,
                        0.0,
                        0.0,
                        0.009,
                        0.0,
                        0.835,
                        0.021,
                        0.012,
                        0.868,
                    ],
                    [
                        1.0,
                        0.006,
                        0.0,
                        0.0,
                        0.006,
                        2.0,
                        0.009,
                        0.0,
                        0.0,
                        0.009,
                        0.0,
                        0.835,
                        0.021,
                        0.012,
                        0.868,
                    ],
                ]
            ).tolist(),
        )


class ProcarSupportFunctionsTestCase(unittest.TestCase):
    """Test for the support functions in procar.py"""

    def test_area_of_a_triangle_in_cartesian_space(self):
        a = np.array([0.0, 0.0, 0.0])
        b = np.array([4.0, 0.0, 0.0])
        c = np.array([0.0, 3.0, 0.0])
        self.assertEqual(procar.area_of_a_triangle_in_cartesian_space(a, b, c), 6.0)

    def test_points_are_in_a_straight_line(self):
        a = np.array([0.0, 0.0, 0.0])
        b = np.array([2.0, 0.0, 0.0])
        c = np.array([3.0, 0.0, 0.0])
        points = [a, b, c]
        tolerance = 1e-7
        self.assertEqual(procar.points_are_in_a_straight_line(points, tolerance), True)

    def test_points_are_not_in_a_straight_line(self):
        a = np.array([0.0, 0.0, 0.0])
        b = np.array([2.0, 1.0, 0.0])
        c = np.array([3.0, 0.0, 0.0])
        points = [a, b, c]
        tolerance = 1e-7
        self.assertEqual(procar.points_are_in_a_straight_line(points, tolerance), False)

    def test_least_squares_effective_mass(self):
        k_points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        eigenvalues = np.array([0.0, 1.0, 4.0])
        with patch(
            "vasppy.procar.points_are_in_a_straight_line"
        ) as mock_straight_line_test:
            mock_straight_line_test.return_value = True
            self.assertAlmostEqual(
                procar.least_squares_effective_mass(k_points, eigenvalues), 13.605693123
            )

    def test_least_squares_effective_mass_with_nonzero_offset(self):
        """Test that effective mass is independent of energy offset."""
        k_points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        # Same parabola as test_least_squares_effective_mass but shifted by +10.0 eV
        eigenvalues = np.array([10.0, 11.0, 14.0])
        with patch(
            "vasppy.procar.points_are_in_a_straight_line"
        ) as mock_straight_line_test:
            mock_straight_line_test.return_value = True
            self.assertAlmostEqual(
                procar.least_squares_effective_mass(k_points, eigenvalues), 13.605693123
            )

    def test_least_squares_effective_mass_raises_valueerror_if_points_are_not_collinear(
        self,
    ):
        k_points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 1.0], [2.0, 0.0, 0.0]])
        eigenvalues = np.array([0.0, 1.0, 4.0])
        with patch(
            "vasppy.procar.points_are_in_a_straight_line"
        ) as mock_straight_line_test:
            mock_straight_line_test.return_value = False
            with self.assertRaises(ValueError):
                procar.least_squares_effective_mass(k_points, eigenvalues)

    def test_two_point_effective_mass(self):
        k_points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        eigenvalues = np.array([0.0, 1.0])
        self.assertAlmostEqual(
            procar.two_point_effective_mass(k_points, eigenvalues), 13.605693123
        )


class KPointReprTestCase(unittest.TestCase):
    """Tests for KPoint.__repr__."""

    def test_repr_contains_index(self):
        kp = procar.KPoint(index=3, frac_coords=np.array([0.1, 0.2, 0.3]), weight=0.5)
        self.assertIn("3", repr(kp))

    def test_repr_contains_weight(self):
        kp = procar.KPoint(index=1, frac_coords=np.array([0.0, 0.0, 0.0]), weight=0.25)
        self.assertIn("0.25", repr(kp))


class TwoPointEffectiveMassRaisesTestCase(unittest.TestCase):
    """Tests that two_point_effective_mass raises for bad inputs."""

    def test_raises_valueerror_for_wrong_kpoint_count(self):
        k_points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        eigenvalues = np.array([0.0, 1.0, 4.0])
        with self.assertRaises(ValueError):
            procar.two_point_effective_mass(k_points, eigenvalues)

    def test_raises_valueerror_for_wrong_eigenvalue_count(self):
        k_points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        eigenvalues = np.array([0.0, 1.0, 2.0])
        with self.assertRaises(ValueError):
            procar.two_point_effective_mass(k_points, eigenvalues)


class ProcarInitTestCase(unittest.TestCase):
    """Tests for Procar.__init__."""

    def _make_minimal_procar(self, **overrides):
        """Build a minimal valid Procar for testing __init__ validation."""
        defaults = dict(
            data=np.zeros((1, 1, 1, 2, 3)),
            bands=np.array([]),
            k_points=[],
            number_of_k_points=1,
            number_of_bands=1,
            number_of_ions=1,
            number_of_projections=3,
            spin_channels=1,
            k_point_blocks=1,
            calculation={
                "non_spin_polarised": True,
                "non_collinear": False,
                "spin_polarised": False,
            },
            negative_occupancies="warn",
        )
        defaults.update(overrides)
        return procar.Procar(**defaults)

    def test_invalid_negative_occupancies_raises(self):
        with self.assertRaises(ValueError):
            self._make_minimal_procar(negative_occupancies="ignore")

    def test_valid_negative_occupancies_warn(self):
        pcar = self._make_minimal_procar(negative_occupancies="warn")
        self.assertEqual(pcar.negative_occupancies, "warn")

    def test_valid_negative_occupancies_zero(self):
        pcar = self._make_minimal_procar(negative_occupancies="zero")
        self.assertEqual(pcar.negative_occupancies, "zero")

    def test_attributes_are_set(self):
        pcar = self._make_minimal_procar()
        self.assertEqual(pcar._number_of_k_points, 1)
        self.assertEqual(pcar._number_of_bands, 1)
        self.assertEqual(pcar._spin_channels, 1)
        self.assertEqual(pcar._number_of_ions, 1)
        self.assertEqual(pcar._k_point_blocks, 1)
        self.assertIsNotNone(pcar._data)


class ProcarSanityCheckTestCase(unittest.TestCase):
    """Tests for Procar.sanity_check."""

    def test_sanity_check_passes_for_valid_data(self):
        pcar = _load_procar()
        # Should not raise
        pcar.sanity_check()

    def test_sanity_check_raises_for_kpoint_mismatch(self):
        pcar = _load_procar()
        pcar._number_of_k_points = 999
        with self.assertRaises(ValueError):
            pcar.sanity_check()

    def test_sanity_check_raises_for_band_mismatch(self):
        pcar = _load_procar()
        pcar._number_of_bands = 999
        with self.assertRaises(ValueError):
            pcar.sanity_check()


class ProcarPropertyMismatchTestCase(unittest.TestCase):
    """Tests that Procar properties raise ValueError on metadata mismatch."""

    def test_number_of_k_points_raises_on_mismatch(self):
        pcar = _load_procar()
        pcar._number_of_k_points = 999
        with self.assertRaises(ValueError):
            _ = pcar.number_of_k_points

    def test_number_of_bands_raises_on_mismatch(self):
        pcar = _load_procar()
        pcar._number_of_bands = 999
        with self.assertRaises(ValueError):
            _ = pcar.number_of_bands

    def test_spin_channels_raises_on_mismatch(self):
        pcar = _load_procar()
        pcar._spin_channels = 999
        with self.assertRaises(ValueError):
            _ = pcar.spin_channels

    def test_number_of_ions_raises_on_mismatch(self):
        pcar = _load_procar()
        pcar._number_of_ions = 999
        with self.assertRaises(ValueError):
            _ = pcar.number_of_ions

    def test_number_of_projections_raises_on_mismatch(self):
        pcar = _load_procar()
        pcar._number_of_projections = 999
        with self.assertRaises(ValueError):
            _ = pcar.number_of_projections


class ProcarAddRaisesTestCase(unittest.TestCase):
    """Tests that __add__ raises ValueError for incompatible Procars."""

    def test_add_raises_for_mismatched_spin_channels(self):
        pcar1 = _load_procar(test_procar_filename)
        pcar2 = _load_procar(test_procar_spin_polarised_filename)
        with self.assertRaises(ValueError):
            _ = pcar1 + pcar2

    def test_add_raises_for_mismatched_ions(self):
        pcar1 = _load_procar()
        pcar2 = deepcopy(pcar1)
        pcar2._number_of_ions = 999
        with self.assertRaises(ValueError):
            _ = pcar1 + pcar2

    def test_add_raises_for_mismatched_bands(self):
        pcar1 = _load_procar()
        pcar2 = deepcopy(pcar1)
        pcar2._number_of_bands = 999
        pcar2._data = pcar2._data[:, :2, :, :, :]
        with self.assertRaises(ValueError):
            _ = pcar1 + pcar2


class ProcarXAxisTestCase(unittest.TestCase):
    """Tests for Procar.x_axis."""

    def test_x_axis_without_reciprocal_lattice(self):
        pcar = _load_procar()
        x = pcar.x_axis()
        self.assertEqual(len(x), pcar.number_of_k_points)
        np.testing.assert_array_equal(x, np.arange(pcar.number_of_k_points))

    def test_x_axis_with_reciprocal_lattice(self):
        pcar = _load_procar()
        recip_lattice = np.eye(3) * 2 * np.pi
        x = pcar.x_axis(recip_lattice)
        self.assertEqual(len(x), pcar.number_of_k_points)
        self.assertEqual(x[0], 0.0)

    def test_x_axis_is_monotonically_increasing(self):
        pcar = _load_procar()
        recip_lattice = np.eye(3) * 2 * np.pi
        x = pcar.x_axis(recip_lattice)
        self.assertTrue(np.all(np.diff(x) >= 0))


class ProcarSelectKPointsTestCase(unittest.TestCase):
    """Tests for Procar.select_k_points."""

    def test_select_k_points_reduces_count(self):
        pcar = _load_procar()
        new_pcar = pcar.select_k_points([0])
        self.assertEqual(new_pcar.number_of_k_points, 1)

    def test_select_k_points_renumbers_indices(self):
        pcar = _load_procar()
        new_pcar = pcar.select_k_points([0, 1])
        self.assertEqual([kp.index for kp in new_pcar._k_points], [1, 2])


class ProcarWeightedBandStructureTestCase(unittest.TestCase):
    """Tests for Procar.weighted_band_structure."""

    def test_weighted_band_structure_shape(self):
        pcar = _load_procar()
        bs = pcar.weighted_band_structure()
        self.assertEqual(bs.shape[0], pcar.number_of_bands)
        self.assertEqual(bs.shape[1], pcar.number_of_k_points)
        self.assertEqual(bs.shape[2], 3)

    def test_weighted_band_structure_with_e_fermi(self):
        pcar = _load_procar()
        bs_no_efermi = pcar.weighted_band_structure(e_fermi=0.0)
        bs_with_efermi = pcar.weighted_band_structure(e_fermi=5.0)
        # Energies should be shifted
        np.testing.assert_array_almost_equal(
            bs_no_efermi[:, :, 1] - 5.0,
            bs_with_efermi[:, :, 1],
        )


if __name__ == "__main__":
    unittest.main()
