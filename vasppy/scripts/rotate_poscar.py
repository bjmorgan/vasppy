#! /usr/bin/env python3

"""Script to rotate the cell lattice of a VASP POSCAR file."""

import argparse
import math

from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp.inputs import Poscar

from vasppy.cell import Cell


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace with poscar path, axis vector, and
        rotation angle in degrees.
    """
    parser = argparse.ArgumentParser(
        description="Rotates the cell lattice in VASP POSCAR files"
    )
    parser.add_argument("poscar", help="filename of the VASP POSCAR to be processed")
    parser.add_argument(
        "-a",
        "--axis",
        nargs=3,
        type=float,
        help="vector for rotation axis",
        required=True,
    )
    parser.add_argument(
        "-d", "--degrees", type=int, help="rotation angle in degrees", required=True
    )
    return parser.parse_args()


def main() -> None:
    """Read a POSCAR, rotate its cell lattice, and print the result."""
    args = parse_command_line_arguments()
    structure = Structure.from_file(args.poscar)
    theta = math.pi * args.degrees / 180.0

    cell = Cell(structure.lattice.matrix.copy())
    cell.rotate(args.axis, theta)

    new_lattice = Lattice(cell.matrix)
    rotated_structure = Structure(
        new_lattice,
        structure.species,
        structure.frac_coords,
    )
    print(Poscar(rotated_structure))


if __name__ == "__main__":
    main()
