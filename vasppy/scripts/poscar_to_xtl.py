#! /usr/bin/env python3
"""Convert a VASP POSCAR file to the XTL file format."""

import argparse

from pymatgen.core import Structure


def poscar_to_xtl_output(filename: str) -> str:
    """Generate XTL-format output from a VASP POSCAR file.

    Reads a POSCAR using pymatgen and produces a string in the XTL crystal
    structure format, including the unit cell parameters, symmetry label, and
    fractional atomic coordinates.

    Args:
        filename: Path to the VASP POSCAR file.

    Returns:
        A string containing the XTL-formatted structure.
    """
    structure = Structure.from_file(filename)
    lattice = structure.lattice

    lines = []
    lines.append(structure.formula)
    lines.append("CELL")
    lines.append(
        "".join(
            [
                "   {: .8f}".format(v)
                for v in (
                    lattice.a,
                    lattice.b,
                    lattice.c,
                    lattice.alpha,
                    lattice.beta,
                    lattice.gamma,
                )
            ]
        )
    )
    lines.append(" Symmetry label P1\n\nATOMS\nNAME      X       Y     Z")

    for site in structure:
        label = site.species_string
        x, y, z = site.frac_coords
        lines.append("{:<6}  {: .10f}  {: .10f}  {: .10f}".format(label, x, y, z))

    return "\n".join(lines) + "\n"


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command-line arguments for poscar_to_xtl.

    Returns:
        Parsed argument namespace with a ``poscar`` attribute.
    """
    parser = argparse.ArgumentParser(
        description="Converts a VASP POSCAR file to the .xtl file format"
    )
    parser.add_argument("poscar", help="filename of the VASP POSCAR to be processed")
    return parser.parse_args()


def main() -> None:
    """Entry point: read a POSCAR file and print XTL output to stdout."""
    args = parse_command_line_arguments()
    print(poscar_to_xtl_output(args.poscar), end="")


if __name__ == "__main__":
    main()
