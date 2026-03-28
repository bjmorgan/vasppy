#! /usr/bin/env python3
"""Generate a weighted (fat) band structure from a VASP PROCAR."""

import sys
if sys.platform != "win32":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

from vasppy import procar
from vasppy.outcar import reciprocal_lattice_from_outcar
import argparse


def orbitals_with_l(l: str) -> list[int] | None:
    """Return orbital indices for a given angular momentum label.

    Args:
        l: Angular momentum label; one of ``'s'``, ``'p'``, ``'d'``,
           ``'f'``, or ``'all'``.

    Returns:
        List of orbital indices, or ``None`` for ``'all'``.
    """
    to_return: dict[str, list[int] | None] = {
        "s": [0],
        "p": [1, 2, 3],
        "d": [4, 5, 6, 7, 8],
        "f": [9, 10, 11, 12, 13],
        "all": None,
    }
    return to_return[l]


def main() -> None:
    """Entry point for the fat_bands command-line script."""
    parser = argparse.ArgumentParser(
        description="Generate a weighted (fat) band structure from a VASP PROCAR"
    )
    parser.add_argument(
        "-i",
        "--ions",
        help="ion indices for band projection (default: sum over all ions)",
        nargs="+",
        type=int,
    )
    parser.add_argument(
        "-s",
        "--spins",
        help="spin indices for band projection (default [ 1 ])",
        nargs="+",
        type=int,
        default=[1],
    )
    parser.add_argument(
        "-o",
        "--orbitals",
        help="orbital indices for band projection (default: sum over all orbitals)",
        nargs="+",
        type=int,
    )
    parser.add_argument(
        "-e",
        "--efermi",
        help="set fermi energy as reference for energy scale",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "-l",
        "--l-angular-momentum",
        help="select all orbitals with angular momentum L for band projection. This supercedes the --orbitals option",
        choices=["s", "p", "d", "f", "all"],
    )
    parser.add_argument(
        "-f",
        "--procar",
        help="PROCAR filename (default PROCAR)",
        type=str,
        default="PROCAR",
    )
    parser.add_argument(
        "--scaling",
        help="Energy scaling for band widths (default 0.2 eV)",
        type=float,
        default=0.2,
    )
    parser.add_argument(
        "-x",
        "--xscaling",
        help="Automatic scaling of x-axis using reciprocal lattice vectors read from OUTCAR",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    if args.l_angular_momentum:
        args.orbitals = orbitals_with_l(args.l_angular_momentum)

    reciprocal_lattice = reciprocal_lattice_from_outcar("OUTCAR") if args.xscaling else None

    pcar = procar.Procar.from_file(args.procar)
    pcar.print_weighted_band_structure(
        spins=args.spins,
        ions=args.ions,
        orbitals=args.orbitals,
        scaling=args.scaling,
        e_fermi=args.efermi,
        reciprocal_lattice=reciprocal_lattice,
    )


if __name__ == "__main__":
    main()
