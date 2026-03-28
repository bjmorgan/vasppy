#! /usr/bin/env python3

"""Script to generate a sqrt(2) x sqrt(2) supercell from a VASP POSCAR."""

import argparse
from typing import Literal

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp.inputs import Poscar

from vasppy.cell import Cell


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace with poscar filename and axis choice.
    """
    parser = argparse.ArgumentParser(
        description="Generate a sqrt(2) x sqrt(2) supercell from a VASP POSCAR"
    )
    parser.add_argument(
        "poscar", help="filename of the VASP POSCAR to be processed"
    )
    parser.add_argument(
        "-a",
        "--axis",
        choices=["x", "y", "z"],
        type=str,
        help="normal vector for sqrt(2) x sqrt(2) expansion",
        required=True,
    )
    return parser.parse_args()


def sqrt2_by_sqrt2_expansion(
    structure: Structure,
    axis: Literal["x", "y", "z"],
) -> Structure:
    """Perform a sqrt(2) x sqrt(2) supercell expansion.

    Rotates the cell by 45 degrees about the specified axis, replicates
    to a 2x2x1 supercell, then adjusts the lattice to an orthorhombic
    cell based on the Cartesian coordinate extents.

    Args:
        structure: The input pymatgen Structure.
        axis: The rotation axis ('x', 'y', or 'z').

    Returns:
        A new Structure with the expanded orthorhombic cell.
    """
    axis_vectors: dict[str, np.ndarray] = {
        "x": np.array([1, 0, 0]),
        "y": np.array([0, 1, 0]),
        "z": np.array([0, 0, 1]),
    }
    # Rotate the cell by 45 degrees about the chosen axis.
    cell = Cell(structure.lattice.matrix.copy())
    cell.rotate(axis=axis_vectors[axis], theta=np.pi / 4)
    rotated_lattice = Lattice(cell.matrix)
    rotated_structure = Structure(
        rotated_lattice,
        structure.species,
        structure.frac_coords,
    )
    # Make a 2x2x1 supercell.
    rotated_structure.make_supercell([2, 2, 1])
    # Extract Cartesian coordinates and build an orthorhombic cell
    # from the diagonal of the current cell matrix.
    cart_coords = rotated_structure.cart_coords
    ortho_matrix = np.diag(rotated_structure.lattice.matrix.diagonal())
    ortho_lattice = Lattice(ortho_matrix)
    # Convert Cartesian coordinates to fractional in the new lattice.
    frac_coords = cart_coords @ np.linalg.inv(ortho_matrix)
    return Structure(
        ortho_lattice,
        rotated_structure.species,
        frac_coords,
    )


def main() -> None:
    """Read a POSCAR, perform sqrt(2) x sqrt(2) expansion, and print the result."""
    args = parse_command_line_arguments()
    structure = Structure.from_file(args.poscar)
    result = sqrt2_by_sqrt2_expansion(structure=structure, axis=args.axis)
    print(Poscar(result))


if __name__ == "__main__":
    main()
