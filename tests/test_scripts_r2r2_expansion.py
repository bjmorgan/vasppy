"""Characterisation tests for r2r2_expansion script."""

import numpy as np
from numpy.testing import assert_allclose
from pymatgen.core import Lattice, Structure

from vasppy.scripts.r2r2_expansion import sqrt2_by_sqrt2_expansion


def _cubic_structure() -> Structure:
    """Create a simple cubic structure with one atom for testing."""
    lattice = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    return Structure(lattice, ["Na"], [[0.0, 0.0, 0.0]])


class TestSqrt2BySqrt2Expansion:
    """Characterisation tests for the sqrt(2) x sqrt(2) expansion."""

    def test_returns_structure(self):
        """The function should return a pymatgen Structure."""
        structure = _cubic_structure()
        result = sqrt2_by_sqrt2_expansion(structure, axis="z")
        assert isinstance(result, Structure)

    def test_lattice_is_orthorhombic(self):
        """The resulting lattice should be orthorhombic (diagonal matrix)."""
        structure = _cubic_structure()
        result = sqrt2_by_sqrt2_expansion(structure, axis="z")
        matrix = result.lattice.matrix
        # Off-diagonal elements should be zero (orthorhombic).
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0.0)
        assert_allclose(off_diag, 0.0, atol=1e-10)

    def test_supercell_atom_count(self):
        """A 2x2x1 supercell should have 4x the original atom count."""
        structure = _cubic_structure()
        result = sqrt2_by_sqrt2_expansion(structure, axis="z")
        assert len(result) == 4 * len(structure)

    def test_lattice_parameters_z_axis(self):
        """For a cubic cell with a=4, rotating about z and expanding 2x2x1
        should give a ~ 4*sqrt(2), b ~ 4*sqrt(2), c ~ 4."""
        structure = _cubic_structure()
        result = sqrt2_by_sqrt2_expansion(structure, axis="z")
        a, b, c = result.lattice.abc
        expected_ab = 4.0 * np.sqrt(2)
        assert_allclose(a, expected_ab, atol=1e-10)
        assert_allclose(b, expected_ab, atol=1e-10)
        assert_allclose(c, 4.0, atol=1e-10)

    def test_fractional_coords_are_finite(self):
        """All fractional coordinates should be finite."""
        structure = _cubic_structure()
        result = sqrt2_by_sqrt2_expansion(structure, axis="z")
        frac = result.frac_coords
        assert np.all(np.isfinite(frac))

    def test_species_preserved(self):
        """Species should be preserved through the expansion."""
        lattice = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
        structure = Structure(lattice, ["Na", "Cl"], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        result = sqrt2_by_sqrt2_expansion(structure, axis="z")
        species_strings = sorted(set(str(s) for s in result.species))
        assert species_strings == ["Cl", "Na"]
        assert len(result) == 8

    def test_x_axis_rotation(self):
        """Rotation about x axis: a is doubled by 2x supercell, b and c
        reflect the rotated-and-diagonalised lattice extents."""
        structure = _cubic_structure()
        result = sqrt2_by_sqrt2_expansion(structure, axis="x")
        a, b, c = result.lattice.abc
        assert_allclose(a, 8.0, atol=1e-10)
        assert_allclose(b, 4.0 * np.sqrt(2), atol=1e-10)
        assert_allclose(c, 2.0 * np.sqrt(2), atol=1e-10)
