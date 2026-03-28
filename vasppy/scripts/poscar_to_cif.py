#! /usr/bin/env python3

import sys
if sys.platform != "win32":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

import argparse

from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description="Converts a VASP POSCAR file to the .cif file format"
    )
    parser.add_argument("poscar", help="filename of the VASP POSCAR to be processed")
    parser.add_argument(
        "-s",
        "--symprec",
        type=float,
        help="Symmetry precision for a symmetrised .cif output",
    )
    return parser.parse_args()


def main() -> None:
    """Read a POSCAR file and print its contents in CIF format."""
    args = parse_command_line_arguments()
    structure = Structure.from_file(args.poscar)
    cif_writer = CifWriter(structure, symprec=args.symprec)
    print(cif_writer)


if __name__ == "__main__":
    main()
