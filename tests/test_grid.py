import numpy as np
import pytest
from pathlib import Path
from vasppy.grid import Grid, interpolate, trilinear_interpolation

CHGCAR_MINIMAL = str(Path(__file__).parent / "test_data" / "CHGCAR_minimal")


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
    def test_default_dimensions(self):
        g = Grid()
        assert g.dimensions == (1, 1, 1)
        assert g.grid.shape == (1, 1, 1)

    def test_custom_dimensions(self):
        g = Grid(dimensions=(3, 4, 5))
        assert g.dimensions == (3, 4, 5)
        assert g.grid.shape == (3, 4, 5)
        np.testing.assert_array_equal(g.grid, np.zeros((3, 4, 5)))

    def test_spacing(self):
        g = Grid(dimensions=(2, 4, 5))
        np.testing.assert_array_almost_equal(g.spacing, [0.5, 0.25, 0.2])

    def test_structure_is_none_by_default(self):
        g = Grid()
        assert g.structure is None

    def test_no_filename_by_default(self):
        g = Grid()
        assert g.filename is None


class TestGridReadFromFilename:
    def test_reads_dimensions(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        assert g.dimensions == (2, 5, 2)

    def test_reads_grid_data(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        assert g.grid.shape == (2, 5, 2)
        assert g.grid[0, 0, 0] == pytest.approx(1.0)
        assert g.grid[1, 0, 0] == pytest.approx(2.0)
        assert g.grid[0, 0, 1] == pytest.approx(11.0)
        assert g.grid[1, 4, 1] == pytest.approx(20.0)

    def test_reads_structure(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        assert g.structure is not None
        np.testing.assert_array_almost_equal(
            g.structure.lattice.matrix,
            [[2.82, 0.0, 0.0], [0.0, 2.82, 0.0], [0.0, 0.0, 2.82]],
        )
        assert len(g.structure) == 2

    def test_returns_self(self):
        g = Grid()
        result = g.read_from_filename(CHGCAR_MINIMAL)
        assert result is g


class TestGridWriteToFilename:
    def test_roundtrip(self, tmp_path):
        """Read a CHGCAR, write it, read it back, and compare grid data."""
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        out_file = str(tmp_path / "CHGCAR_out")
        g.write_to_filename(out_file)

        g2 = Grid().read_from_filename(out_file)
        assert g2.dimensions == g.dimensions
        np.testing.assert_array_almost_equal(g2.grid, g.grid)

    def test_output_contains_lattice(self, tmp_path):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        out_file = str(tmp_path / "CHGCAR_out")
        g.write_to_filename(out_file)

        with open(out_file) as f:
            content = f.read()
        assert "2.8199" in content

    def test_output_contains_dimensions(self, tmp_path):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        out_file = str(tmp_path / "CHGCAR_out")
        g.write_to_filename(out_file)

        with open(out_file) as f:
            lines = f.readlines()
        # Find the dimensions line (after POSCAR block + blank line)
        dim_line = None
        for line in lines:
            parts = line.split()
            if parts == ["2", "5", "2"]:
                dim_line = line
                break
        assert dim_line is not None


class TestGridAverage:
    def test_average_z(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        avg = g.average("z")
        np.testing.assert_array_almost_equal(avg, [5.5, 15.5])

    def test_average_x(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        avg = g.average("x")
        np.testing.assert_array_almost_equal(avg, [10.0, 11.0])


class TestGridCoordinateMethods:
    def test_fractional_coordinate_at_index(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        frac = g.fractional_coordinate_at_index([1, 0, 0])
        np.testing.assert_array_almost_equal(frac, [0.5, 0.0, 0.0])

    def test_cartesian_coordinate_at_index(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        cart = g.cartesian_coordinate_at_index([1, 0, 0])
        np.testing.assert_array_almost_equal(cart, [1.41, 0.0, 0.0])


class TestGridInterpolation:
    def test_interpolated_value_at_fractional_coordinate(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        val = g.interpolated_value_at_fractional_coordinate([0.25, 0.25, 0.25])
        assert val == pytest.approx(9.0)


class TestGridCubeSlice:
    def test_cube_slice_at_origin(self):
        g = Grid().read_from_filename(CHGCAR_MINIMAL)
        cube = g.cube_slice(0, 0, 0)
        assert cube.shape == (2, 2, 2)
        assert cube[0, 0, 0] == g.grid[0, 0, 0]
        assert cube[1, 0, 0] == g.grid[1, 0, 0]
