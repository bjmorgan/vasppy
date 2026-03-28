import numpy as np
import pytest
from pathlib import Path
from pymatgen.core import Lattice, Structure
from vasppy.grid import Grid, interpolate, trilinear_interpolation

CHGCAR_MINIMAL = str(Path(__file__).parent / "test_data" / "CHGCAR_minimal")


def _dummy_structure() -> Structure:
    """Return a simple cubic structure for testing."""
    return Structure(
        Lattice.cubic(2.82),
        ["Na", "Cl"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )


class TestInterpolate:
    def test_interpolate_at_zero(self):
        assert interpolate(3.0, 7.0, 0.0) == 3.0

    def test_interpolate_at_one(self):
        assert interpolate(3.0, 7.0, 1.0) == 7.0

    def test_interpolate_at_half(self):
        assert interpolate(3.0, 7.0, 0.5) == 5.0


class TestTrilinearInterpolation:
    def test_at_origin(self):
        cube = np.arange(8).reshape((2, 2, 2), order="C").astype(float)
        result = trilinear_interpolation(cube, [0, 0, 0])
        assert result == cube[0, 0, 0]

    def test_at_centre(self):
        cube = np.ones((2, 2, 2))
        result = trilinear_interpolation(cube, [0.5, 0.5, 0.5])
        assert result == pytest.approx(1.0)


class TestGridInit:
    def test_requires_structure(self):
        with pytest.raises(TypeError):
            Grid(dimensions=(2, 2, 2), grid=np.zeros((2, 2, 2)))

    def test_requires_dimensions(self):
        with pytest.raises(TypeError):
            Grid(structure=_dummy_structure(), grid=np.zeros((2, 2, 2)))

    def test_requires_grid(self):
        with pytest.raises(TypeError):
            Grid(structure=_dummy_structure(), dimensions=(2, 2, 2))

    def test_stores_structure(self):
        s = _dummy_structure()
        g = Grid(structure=s, dimensions=(2, 2, 2), grid=np.zeros((2, 2, 2)))
        assert g.structure is s

    def test_stores_dimensions(self):
        g = Grid(
            structure=_dummy_structure(),
            dimensions=(3, 4, 5),
            grid=np.zeros((3, 4, 5)),
        )
        assert g.dimensions == (3, 4, 5)

    def test_stores_grid(self):
        data = np.ones((3, 4, 5))
        g = Grid(structure=_dummy_structure(), dimensions=(3, 4, 5), grid=data)
        np.testing.assert_array_equal(g.grid, data)

    def test_spacing(self):
        g = Grid(
            structure=_dummy_structure(),
            dimensions=(2, 4, 5),
            grid=np.zeros((2, 4, 5)),
        )
        np.testing.assert_array_almost_equal(g.spacing, [0.5, 0.25, 0.2])

    def test_no_filename_attribute(self):
        g = Grid(
            structure=_dummy_structure(),
            dimensions=(2, 2, 2),
            grid=np.zeros((2, 2, 2)),
        )
        assert not hasattr(g, "filename")


class TestGridFromFile:
    def test_reads_dimensions(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        assert g.dimensions == (2, 5, 2)

    def test_reads_grid_data(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        assert g.grid.shape == (2, 5, 2)
        assert g.grid[0, 0, 0] == pytest.approx(1.0)
        assert g.grid[1, 0, 0] == pytest.approx(2.0)
        assert g.grid[0, 0, 1] == pytest.approx(11.0)
        assert g.grid[1, 4, 1] == pytest.approx(20.0)

    def test_reads_structure(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        assert g.structure is not None
        np.testing.assert_array_almost_equal(
            g.structure.lattice.matrix,
            [[2.82, 0.0, 0.0], [0.0, 2.82, 0.0], [0.0, 0.0, 2.82]],
        )
        assert len(g.structure) == 2

    def test_returns_grid_instance(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        assert isinstance(g, Grid)


class TestGridWriteToFilename:
    def test_roundtrip(self, tmp_path):
        """Read a CHGCAR, write it, read it back, and compare grid data."""
        g = Grid.from_file(CHGCAR_MINIMAL)
        out_file = str(tmp_path / "CHGCAR_out")
        g.write_to_filename(out_file)

        g2 = Grid.from_file(out_file)
        assert g2.dimensions == g.dimensions
        np.testing.assert_array_almost_equal(g2.grid, g.grid)

    def test_output_contains_lattice(self, tmp_path):
        g = Grid.from_file(CHGCAR_MINIMAL)
        out_file = str(tmp_path / "CHGCAR_out")
        g.write_to_filename(out_file)

        with open(out_file) as f:
            content = f.read()
        assert "2.8199" in content

    def test_output_contains_dimensions(self, tmp_path):
        g = Grid.from_file(CHGCAR_MINIMAL)
        out_file = str(tmp_path / "CHGCAR_out")
        g.write_to_filename(out_file)

        with open(out_file) as f:
            lines = f.readlines()
        dim_line = None
        for line in lines:
            parts = line.split()
            if parts == ["2", "5", "2"]:
                dim_line = line
                break
        assert dim_line is not None


class TestGridAverage:
    def test_average_z(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        avg = g.average("z")
        np.testing.assert_array_almost_equal(avg, [5.5, 15.5])

    def test_average_x(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        avg = g.average("x")
        np.testing.assert_array_almost_equal(avg, [10.0, 11.0])


class TestGridCoordinateMethods:
    def test_fractional_coordinate_at_index(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        frac = g.fractional_coordinate_at_index([1, 0, 0])
        np.testing.assert_array_almost_equal(frac, [0.5, 0.0, 0.0])

    def test_cartesian_coordinate_at_index(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        cart = g.cartesian_coordinate_at_index([1, 0, 0])
        np.testing.assert_array_almost_equal(cart, [1.41, 0.0, 0.0])


class TestGridInterpolation:
    def test_interpolated_value_at_fractional_coordinate(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        val = g.interpolated_value_at_fractional_coordinate([0.25, 0.25, 0.25])
        assert val == pytest.approx(9.0)


class TestGridCubeSlice:
    def test_cube_slice_at_origin(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        cube = g.cube_slice(0, 0, 0)
        assert cube.shape == (2, 2, 2)
        assert cube[0, 0, 0] == g.grid[0, 0, 0]
        assert cube[1, 0, 0] == g.grid[1, 0, 0]


class TestGridByIndex:
    def test_by_index_returns_correct_value(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        assert g.by_index([0, 0, 0]) == pytest.approx(1.0)
        assert g.by_index([1, 0, 0]) == pytest.approx(2.0)

    def test_by_index_with_numpy_array(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        idx = np.array([0, 0, 1])
        assert g.by_index(idx) == pytest.approx(11.0)


class TestGridAverageAxes:
    """Explicit tests verifying the planar average divisor is correct
    for non-cubic grids where nx != ny != nz."""

    def test_average_y_axis(self):
        """Average along y should divide by nx*nz, not nx*ny."""
        g = Grid.from_file(CHGCAR_MINIMAL)
        # Grid shape (2, 5, 2); average along y -> sum over axes 0 and 2, divide by 2*2=4
        avg = g.average("y")
        assert avg.shape == (5,)
        # Manual check: grid[:, 0, :].sum() / 4 = (1+11+2+12)/4 = 26/4 = 6.5
        assert avg[0] == pytest.approx(6.5)

    def test_average_invalid_axis_raises(self):
        g = Grid.from_file(CHGCAR_MINIMAL)
        with pytest.raises(KeyError):
            g.average("w")
