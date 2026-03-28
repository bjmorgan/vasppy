"""Module for reading and manipulating VASP volumetric data (CHGCAR/LOCPOT)."""

from __future__ import annotations

import math
from typing import ClassVar, cast

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp.inputs import Poscar as PmgPoscar


def interpolate(i: float, j: float, x: float) -> float:
    """Linearly interpolate between two values.

    Args:
        i: Value at x=0.
        j: Value at x=1.
        x: Interpolation parameter in [0, 1].

    Returns:
        The interpolated value.
    """
    return (i * (1.0 - x)) + (j * x)


def trilinear_interpolation(cube: np.ndarray, r: np.ndarray) -> float:
    """Trilinearly interpolate within a 2x2x2 cube of values.

    Args:
        cube: A (2, 2, 2) array of values at cube vertices.
        r: Fractional position within the cube, shape (3,).

    Returns:
        The interpolated value.
    """
    return interpolate(
        interpolate(
            interpolate(cube[0, 0, 0], cube[1, 0, 0], r[0]),
            interpolate(cube[0, 1, 0], cube[1, 1, 0], r[0]),
            r[1],
        ),
        interpolate(
            interpolate(cube[0, 0, 1], cube[1, 0, 1], r[0]),
            interpolate(cube[0, 1, 1], cube[1, 1, 1], r[0]),
            r[1],
        ),
        r[2],
    )


def _read_poscar_header(filename: str) -> tuple[Structure, int]:
    """Read the POSCAR header from a VASP volumetric file.

    Parses the header to determine the number of header lines and
    constructs a pymatgen Structure from the POSCAR portion.

    Args:
        filename: Path to the VASP volumetric file.

    Returns:
        A tuple of (structure, n_header_lines) where n_header_lines
        is the line index of the grid dimensions line (i.e. the
        POSCAR block plus the blank separator line).
    """
    with open(filename) as f:
        lines = f.readlines()
    # Line 0: title
    # Line 1: scaling
    # Lines 2-4: lattice vectors
    # Line 5: species names
    # Line 6: species counts
    n_atoms = sum(int(x) for x in lines[6].split())
    offset = 7
    # Check for selective dynamics
    if lines[offset].strip()[0] in "sS":
        offset += 1
    # Coordinate type line
    offset += 1
    # Atom coordinate lines
    offset += n_atoms
    poscar_str = "".join(lines[:offset])
    pmg_poscar = PmgPoscar.from_str(poscar_str)
    # Add 1 for the blank line between POSCAR block and grid dimensions
    return pmg_poscar.structure, offset + 1


def _read_dimensions(filename: str, n_header_lines: int) -> tuple[int, int, int]:
    """Read grid dimensions from the volumetric file.

    Args:
        filename: Path to the VASP volumetric file.
        n_header_lines: Number of header lines before the dimensions line.

    Returns:
        Grid dimensions as (nx, ny, nz).
    """
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == n_header_lines:
                nx, ny, nz = (int(x) for x in line.split())
                return (nx, ny, nz)
    raise ValueError(f"Could not read dimensions from {filename}")


def _read_grid(
    filename: str,
    n_header_lines: int,
    dimensions: tuple[int, int, int],
) -> np.ndarray:
    """Read grid data values from the volumetric file.

    Args:
        filename: Path to the VASP volumetric file.
        n_header_lines: Number of header lines before the dimensions line.
        dimensions: Grid dimensions (nx, ny, nz).

    Returns:
        3D numpy array of grid data in shape ``dimensions``.
    """
    total_points = math.prod(dimensions)
    grid_data_lines = math.ceil(total_points / 5)
    lines: list[str] = []
    with open(filename) as f:
        for i, line in enumerate(f):
            if (i > n_header_lines) and (
                i <= n_header_lines + grid_data_lines
            ):
                lines.append(line.strip())
    values = np.array(" ".join(lines).split(), dtype=float)
    return values.reshape(dimensions, order="F")


class Grid:
    """Represents volumetric data on a regular grid from VASP CHGCAR/LOCPOT files.

    Every ``Grid`` instance is fully initialised at construction: it always
    has a valid ``structure``, ``dimensions``, and ``grid`` array.

    Attributes:
        projections: Mapping from axis labels to indices.
        structure: pymatgen Structure from the POSCAR header.
        dimensions: Grid dimensions (nx, ny, nz).
        spacing: Fractional spacing along each axis.
        grid: 3D numpy array of grid data.
    """

    projections: ClassVar[dict[str, int]] = {"x": 0, "y": 1, "z": 2}

    def __init__(
        self,
        structure: Structure,
        dimensions: tuple[int, int, int],
        grid: np.ndarray,
    ) -> None:
        """Initialise a Grid object.

        Args:
            structure: pymatgen Structure describing the unit cell.
            dimensions: Grid dimensions as (nx, ny, nz).
            grid: 3D numpy array of volumetric data with shape
                matching ``dimensions``.
        """
        self.structure = structure
        self.dimensions = dimensions
        self.spacing = np.array([1.0 / n for n in self.dimensions])
        self.grid = grid

    @classmethod
    def from_file(cls, filename: str) -> Grid:
        """Read volumetric data from a VASP CHGCAR/LOCPOT file.

        Args:
            filename: Path to the file.

        Returns:
            A fully-constructed Grid instance.
        """
        structure, n_header_lines = _read_poscar_header(filename)
        dimensions = _read_dimensions(filename, n_header_lines)
        grid_data = _read_grid(filename, n_header_lines, dimensions)
        return cls(structure=structure, dimensions=dimensions, grid=grid_data)

    def write_to_filename(self, filename: str) -> None:
        """Write the volumetric data to a file in VASP CHGCAR format.

        Grid values are written in Fortran column-major order with five
        values per line, matching the VASP convention.

        Args:
            filename: Path to the output file.
        """
        with open(filename, "w") as f:
            poscar_str = PmgPoscar(self.structure).get_str()
            f.write(poscar_str)
            f.write(f"\n{' '.join(str(i) for i in self.dimensions)}\n")
            np.savetxt(
                f,
                np.swapaxes(self.grid, 0, 2).reshape(-1, 5),
                fmt="%.11E",
            )

    def average(self, normal_axis_label: str) -> np.ndarray:
        """Calculate the planar average perpendicular to a given axis.

        Args:
            normal_axis_label: Axis label (``'x'``, ``'y'``, or ``'z'``).

        Returns:
            1D array of averaged values along the specified axis.

        Raises:
            KeyError: If ``normal_axis_label`` is not ``'x'``, ``'y'``,
                or ``'z'``.
        """
        axes = [0, 1, 2]
        axes.remove(Grid.projections[normal_axis_label])
        n_plane = self.dimensions[axes[0]] * self.dimensions[axes[1]]
        return cast(np.ndarray, np.sum(np.sum(self.grid, axis=axes[1]), axis=axes[0]) / n_plane)

    def by_index(self, index: list[int] | np.ndarray) -> float:
        """Return the grid value at a given index.

        Args:
            index: Three-element list or array ``[i, j, k]``.

        Returns:
            The grid value at that index.
        """
        return float(self.grid[index[0], index[1], index[2]])

    def fractional_coordinate_at_index(self, index: np.ndarray | list) -> np.ndarray:
        """Convert a grid index to fractional coordinates.

        Args:
            index: Three-element array or list of grid indices.

        Returns:
            Fractional coordinates as a numpy array.
        """
        return cast(np.ndarray, np.multiply(self.spacing, index))

    def cartesian_coordinate_at_index(self, index: np.ndarray | list) -> np.ndarray:
        """Convert a grid index to Cartesian coordinates.

        Args:
            index: Three-element array or list of grid indices.

        Returns:
            Cartesian coordinates as a numpy array.
        """
        return cast(np.ndarray, self.fractional_coordinate_at_index(index).dot(
            self.structure.lattice.matrix
        ))

    def cube_slice(self, x0: int, y0: int, z0: int) -> np.ndarray:
        """Extract a 2x2x2 cube of grid values around a point.

        Wraps around periodic boundaries.

        Args:
            x0: x index of the lower corner.
            y0: y index of the lower corner.
            z0: z index of the lower corner.

        Returns:
            A (2, 2, 2) numpy array of grid values.
        """
        x1 = (x0 + 1) % self.dimensions[0]
        y1 = (y0 + 1) % self.dimensions[1]
        z1 = (z0 + 1) % self.dimensions[2]
        return np.array([
            self.grid[x0, y0, z0], self.grid[x0, y0, z1],
            self.grid[x0, y1, z0], self.grid[x0, y1, z1],
            self.grid[x1, y0, z0], self.grid[x1, y0, z1],
            self.grid[x1, y1, z0], self.grid[x1, y1, z1],
        ]).reshape((2, 2, 2))

    def interpolated_value_at_fractional_coordinate(
        self, coord: np.ndarray | list,
    ) -> float:
        """Interpolate the grid value at an arbitrary fractional coordinate.

        Uses trilinear interpolation from the eight surrounding grid points.

        Args:
            coord: Fractional coordinates [x, y, z].

        Returns:
            The interpolated value.
        """
        point = np.array(self.dimensions, dtype=float) * np.asarray(coord)
        origin = point.astype(int)
        delta = point - origin
        cube = self.cube_slice(*origin)
        return trilinear_interpolation(cube, delta)

    def interpolate_to_orthorhombic_grid(
        self, dimensions: tuple[int, int, int],
    ) -> Grid:
        """Interpolate grid data onto an orthorhombic grid.

        Creates a new Grid with a diagonal lattice matrix (using only the
        diagonal elements of the current lattice) and interpolates all
        values from the original grid.

        Args:
            dimensions: Dimensions for the new orthorhombic grid.

        Returns:
            A new Grid with the interpolated data.

        Note:
            This may need a more robust minimum image function for highly
            non-orthorhombic cells.
        """
        old_lattice = self.structure.lattice
        new_matrix = np.diag(np.diag(old_lattice.matrix))
        new_lattice = Lattice(new_matrix)
        # Create a dummy structure for the new grid
        new_structure = Structure(
            new_lattice, ["X"], [[0, 0, 0]],
        )
        new_grid = Grid(
            structure=new_structure,
            dimensions=dimensions,
            grid=np.zeros(dimensions, dtype=float),
        )
        index_grid = np.array(
            [[i, j, k] for (i, j, k), _ in np.ndenumerate(new_grid.grid)]
        )
        cart_coord_grid = np.array(
            [new_grid.cartesian_coordinate_at_index(idx) for idx in index_grid]
        )
        init_frac_coord_grid = old_lattice.get_fractional_coords(cart_coord_grid)
        frac_coord_grid = init_frac_coord_grid % 1.0
        new_grid_data = np.array([
            self.interpolated_value_at_fractional_coordinate(r)
            for r in frac_coord_grid
        ])
        new_grid.grid = new_grid_data.reshape(new_grid.dimensions)
        return new_grid
