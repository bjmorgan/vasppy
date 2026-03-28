#! /usr/bin/env python3

import sys
if sys.platform != "win32":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

import argparse

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description="Finds the spacegroup for a VASP POSCAR file"
    )
    parser.add_argument("poscar", help="filename of the VASP POSCAR to be processed")
    parser.add_argument(
        "-s",
        "--symprec",
        type=float,
        help="Precision for symmetry analysis (default=1e-3)",
        default=1e-3,
    )
    return parser.parse_args()


def main() -> None:
    """Read a POSCAR file and print the space group symbol."""
    args = parse_command_line_arguments()
    structure = Structure.from_file(args.poscar)
    symmetry_analyzer = SpacegroupAnalyzer(structure, symprec=args.symprec)
    print(symmetry_analyzer.get_space_group_symbol())


if __name__ == "__main__":
    main()
