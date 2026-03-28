import unittest
import numpy as np
from pymatgen.core import Structure, Lattice
from unittest.mock import Mock, patch, call

from vasppy.rdf import (
    RadialDistributionFunction,
    VanHoveAnalysis,
    shell_volumes,
)


class TestShellVolumes(unittest.TestCase):
    """Tests for shell_volumes."""

    def test_shell_volumes_single_shell(self):
        """Test volume of a single spherical shell."""
        intervals = np.array([0.0, 1.0])
        expected = 4.0 / 3.0 * np.pi * (1.0**3 - 0.0**3)
        np.testing.assert_allclose(shell_volumes(intervals), [expected])

    def test_shell_volumes_multiple_shells(self):
        """Test that volumes are strictly positive for increasing intervals."""
        intervals = np.linspace(0.0, 10.0, 11)
        vols = shell_volumes(intervals)
        self.assertEqual(len(vols), 10)
        self.assertTrue(np.all(vols > 0))

    def test_shell_volumes_are_increasing(self):
        """Test that outer shells have greater volume than inner shells."""
        intervals = np.linspace(0.0, 5.0, 6)
        vols = shell_volumes(intervals)
        self.assertTrue(np.all(np.diff(vols) > 0))

    def test_shell_volumes_inner_sphere_zero(self):
        """Total volume of all shells equals sphere of radius r_max."""
        intervals = np.linspace(0.0, 3.0, 4)
        total = np.sum(shell_volumes(intervals))
        expected = 4.0 / 3.0 * np.pi * 3.0**3
        self.assertAlmostEqual(total, expected)


class TestRadialDistributionFunction(unittest.TestCase):

    def test_RadialDistributionFunction_raises_ValueError_if_weights_doesnt_match_structures(
        self,
    ):
        mock_structures = [Mock(spec=Structure), Mock(spec=Structure)]
        weights = [1, 2, 3]
        indices_i = [0, 1]
        with self.assertRaises(ValueError):
            RadialDistributionFunction(
                structures=mock_structures, indices_i=indices_i, weights=weights
            )

    def test_RadialDistributionFunction_init(self):
        mock_structures = [Mock(spec=Structure), Mock(spec=Structure)]
        for s in mock_structures:
            s.lattice = Mock(spec=Lattice)
            s.lattice.volume = 1.0
        indices_i = [0, 1]
        with patch("vasppy.rdf.dr_ij") as mock_dr_ij:
            mock_dr_ij.side_effect = [np.array([5.0, 6.0]), np.array([6.0, 7.0])]
            with patch("vasppy.rdf.shell_volumes") as mock_shell_volumes:
                mock_shell_volumes.return_value = np.ones(500)
                rdf = RadialDistributionFunction(
                    structures=mock_structures, indices_i=indices_i
                )
        self.assertEqual(rdf.indices_i, [0, 1])
        self.assertEqual(rdf.indices_j, [0, 1])
        self.assertEqual(rdf.nbins, 500)
        self.assertEqual(rdf.range, (0.0, 10.0))
        np.testing.assert_array_equal(rdf.intervals, np.linspace(0, 10, 501))
        self.assertEqual(rdf.dr, 0.02)
        np.testing.assert_array_equal(rdf.r, np.linspace(0.01, 9.99, 500))
        expected_rdf = np.zeros_like(rdf.r)
        expected_rdf[250] = 0.125
        expected_rdf[300] = 0.25
        expected_rdf[350] = 0.125
        expected_coordination_number = np.cumsum(expected_rdf) * 2.0
        np.testing.assert_array_almost_equal(rdf.rdf, expected_rdf)
        np.testing.assert_array_almost_equal(
            rdf.coordination_number, expected_coordination_number
        )
        expected_calls = [
            call(
                structure=mock_structures[0],
                indices_i=[0, 1],
                indices_j=[0, 1],
                self_reference=False,
            ),
            call(
                structure=mock_structures[1],
                indices_i=[0, 1],
                indices_j=[0, 1],
                self_reference=False,
            ),
        ]
        mock_dr_ij.assert_has_calls(expected_calls)

    def test_RadialDistributionFunction_accepts_numpy_array_indices(self):
        mock_structures = [Mock(spec=Structure)]
        mock_structures[0].lattice = Mock(spec=Lattice)
        mock_structures[0].lattice.volume = 1.0
        indices_i = np.array([0, 1])
        indices_j = np.array([0, 1])
        with patch("vasppy.rdf.dr_ij") as mock_dr_ij:
            mock_dr_ij.return_value = np.array([5.0, 6.0])
            with patch("vasppy.rdf.shell_volumes") as mock_shell_volumes:
                mock_shell_volumes.return_value = np.ones(500)
                # Should not raise ValueError: ambiguous truth value of an array
                rdf = RadialDistributionFunction(
                    structures=mock_structures,
                    indices_i=indices_i,
                    indices_j=indices_j,
                )
        self.assertIsInstance(rdf.indices_i, list)
        self.assertIsInstance(rdf.indices_j, list)
        self.assertEqual(rdf.indices_i, [0, 1])
        self.assertEqual(rdf.indices_j, [0, 1])

    def test_RadialDistributionFunction_accepts_numpy_array_indices_j_none(self):
        mock_structures = [Mock(spec=Structure)]
        mock_structures[0].lattice = Mock(spec=Lattice)
        mock_structures[0].lattice.volume = 1.0
        indices_i = np.array([0, 1])
        with patch("vasppy.rdf.dr_ij") as mock_dr_ij:
            mock_dr_ij.return_value = np.array([5.0, 6.0])
            with patch("vasppy.rdf.shell_volumes") as mock_shell_volumes:
                mock_shell_volumes.return_value = np.ones(500)
                rdf = RadialDistributionFunction(
                    structures=mock_structures,
                    indices_i=indices_i,
                )
        self.assertIsInstance(rdf.indices_i, list)
        self.assertIsInstance(rdf.indices_j, list)
        self.assertEqual(rdf.indices_i, [0, 1])
        self.assertEqual(rdf.indices_j, [0, 1])

    def test_smeared_rdf_returns_array_of_same_length(self):
        """Test that smeared_rdf returns an array of the same length as rdf."""
        mock_structures = [Mock(spec=Structure)]
        mock_structures[0].lattice = Mock(spec=Lattice)
        mock_structures[0].lattice.volume = 1.0
        indices_i = [0, 1]
        with patch("vasppy.rdf.dr_ij") as mock_dr_ij:
            mock_dr_ij.return_value = np.array([5.0, 6.0])
            with patch("vasppy.rdf.shell_volumes") as mock_shell_volumes:
                mock_shell_volumes.return_value = np.ones(500)
                rdf = RadialDistributionFunction(
                    structures=mock_structures, indices_i=indices_i
                )
        smeared = rdf.smeared_rdf(sigma=0.1)
        self.assertEqual(smeared.shape, rdf.rdf.shape)

    def test_smeared_rdf_broadens_peaks(self):
        """Test that smearing reduces the maximum value of the RDF."""
        mock_structures = [Mock(spec=Structure)]
        mock_structures[0].lattice = Mock(spec=Lattice)
        mock_structures[0].lattice.volume = 1.0
        indices_i = [0, 1]
        with patch("vasppy.rdf.dr_ij") as mock_dr_ij:
            mock_dr_ij.return_value = np.array([5.0])
            with patch("vasppy.rdf.shell_volumes") as mock_shell_volumes:
                mock_shell_volumes.return_value = np.ones(500)
                rdf = RadialDistributionFunction(
                    structures=mock_structures, indices_i=indices_i
                )
        smeared = rdf.smeared_rdf(sigma=0.2)
        self.assertLess(smeared.max(), rdf.rdf.max())

    def test_from_species_strings_raises_if_species_i_not_found(self):
        """Test that from_species_strings raises ValueError for absent species_i."""
        lattice = Lattice.cubic(4.0)
        structure = Structure(lattice, ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        with self.assertRaises(ValueError):
            RadialDistributionFunction.from_species_strings(
                structures=[structure], species_i="Fe", species_j="Na"
            )

    def test_from_species_strings_raises_if_species_j_not_found(self):
        """Test that from_species_strings raises ValueError for absent species_j."""
        lattice = Lattice.cubic(4.0)
        structure = Structure(lattice, ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        with self.assertRaises(ValueError):
            RadialDistributionFunction.from_species_strings(
                structures=[structure], species_i="Na", species_j="Fe"
            )

    def test_from_species_strings_returns_rdf_instance(self):
        """Test that from_species_strings constructs an RDF correctly."""
        lattice = Lattice.cubic(4.0)
        structure = Structure(lattice, ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        with patch("vasppy.rdf.dr_ij") as mock_dr_ij:
            mock_dr_ij.return_value = np.array([3.46])
            rdf = RadialDistributionFunction.from_species_strings(
                structures=[structure], species_i="Na", species_j="Cl"
            )
        self.assertIsInstance(rdf, RadialDistributionFunction)
        self.assertEqual(rdf.indices_i, [0])
        self.assertEqual(rdf.indices_j, [1])


class TestVanHoveAnalysis(unittest.TestCase):
    """Tests for VanHoveAnalysis."""

    def _make_mock_structures(self, n: int = 3) -> list:
        """Return a list of n mock Structure objects sharing a common lattice."""
        lattice = Mock(spec=Lattice)
        lattice.volume = 8.0
        frac_coords = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        n_atoms = len(frac_coords)
        structures = []
        for _ in range(n):
            s = Mock(spec=Structure)
            s.lattice = lattice
            s.frac_coords = frac_coords
            structures.append(s)
        # All pairwise distances set to a large value so no bin is hit
        lattice.get_all_distances.return_value = np.full((n_atoms, n_atoms), 20.0)
        return structures

    def test_vanhove_init_returns_arrays_of_correct_shape(self):
        """Test that gsrt and gdrt are initialised to arrays of length nbins."""
        structures = self._make_mock_structures(n=4)
        indices = [0, 1]
        vha = VanHoveAnalysis(structures=structures, indices=indices, d_steps=1, nbins=100)
        self.assertEqual(vha.gsrt.shape, (100,))
        self.assertEqual(vha.gdrt.shape, (100,))

    def test_vanhove_init_stores_correct_metadata(self):
        """Test that metadata attributes are stored correctly."""
        structures = self._make_mock_structures(n=3)
        indices = [0, 1]
        vha = VanHoveAnalysis(
            structures=structures, indices=indices, d_steps=1,
            nbins=200, r_min=0.0, r_max=5.0,
        )
        self.assertEqual(vha.nbins, 200)
        self.assertEqual(vha.range, (0.0, 5.0))
        self.assertAlmostEqual(vha.dr, 5.0 / 200)

    def test_vanhove_self_returns_gsrt_without_sigma(self):
        """Test that self() without sigma returns gsrt."""
        structures = self._make_mock_structures(n=3)
        vha = VanHoveAnalysis(structures=structures, indices=[0], d_steps=1)
        np.testing.assert_array_equal(vha.self(), vha.gsrt)

    def test_vanhove_distinct_returns_gdrt_without_sigma(self):
        """Test that distinct() without sigma returns gdrt."""
        structures = self._make_mock_structures(n=3)
        vha = VanHoveAnalysis(structures=structures, indices=[0], d_steps=1)
        np.testing.assert_array_equal(vha.distinct(), vha.gdrt)

    def test_vanhove_self_returns_smeared_with_sigma(self):
        """Test that self() with sigma returns a smeared array."""
        structures = self._make_mock_structures(n=3)
        vha = VanHoveAnalysis(structures=structures, indices=[0, 1], d_steps=1)
        # Introduce a non-zero spike so smearing produces a different result
        vha.gsrt[50] = 1.0
        smeared = vha.self(sigma=0.5)
        with self.assertRaises(AssertionError):
            np.testing.assert_array_equal(smeared, vha.gsrt)

    def test_vanhove_distinct_returns_smeared_with_sigma(self):
        """Test that distinct() with sigma returns a smeared array."""
        structures = self._make_mock_structures(n=3)
        vha = VanHoveAnalysis(structures=structures, indices=[0, 1], d_steps=1)
        vha.gdrt[50] = 1.0
        smeared = vha.distinct(sigma=0.5)
        with self.assertRaises(AssertionError):
            np.testing.assert_array_equal(smeared, vha.gdrt)


if __name__ == "__main__":
    unittest.main()
