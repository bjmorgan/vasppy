#! /usr/bin/env python3
"""Calculate an effective mass from a VASP PROCAR using a fitted quadratic."""

import sys
if sys.platform != "win32":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

from typing import Any, Sequence

from vasppy import procar
from vasppy.outcar import reciprocal_lattice_from_outcar
import argparse


def minimum_length(nmin: int) -> type[argparse.Action]:
    """Return an argparse Action class that enforces a minimum argument count.

    Args:
        nmin: Minimum number of arguments required.

    Returns:
        An argparse Action subclass that raises an error if fewer than
        ``nmin`` arguments are supplied.
    """
    class MinimumLength(argparse.Action):
        def __call__(
            self,
            parser: argparse.ArgumentParser,
            args: argparse.Namespace,
            values: str | Sequence[Any] | None,
            option_string: str | None = None,
        ) -> None:
            if values is None or not nmin <= len(values):  # type: ignore[arg-type]
                msg = f'argument "{self.dest}" requires at least {nmin} arguments'
                raise argparse.ArgumentError(self, msg)
            setattr(args, self.dest, values)

    return MinimumLength


def main() -> None:
    """Entry point for the effective_mass command-line script."""
    parser = argparse.ArgumentParser(
        description="Calculate an effective mass from a VASP PROCAR using a fitted quadratic"
    )
    parser.add_argument(
        "-k",
        "--k-points",
        help="index of k-points for calculating effective mass",
        nargs="+",
        type=int,
        required=True,
        action=minimum_length(2),
    )
    parser.add_argument(
        "-b",
        "--band-index",
        help="index of band for calculating effective mass",
        type=int,
        required=True,
    )
    parser.add_argument(
        "-f",
        "--procar",
        help="PROCAR filename (default PROCAR)",
        type=str,
        default="PROCAR",
    )
    parser.add_argument("-v", "--verbose", help="Verbose output", action="store_true")
    parser.add_argument(
        "-o",
        "--outcar",
        help="OUTCAR filename (default OUTCAR)",
        type=str,
        default="OUTCAR",
    )
    parser.add_argument(
        "-s",
        "--spin",
        help="select spin channel (default 1 / non-spin-polarised)",
        type=int,
        default="1",
    )
    args = parser.parse_args()

    reciprocal_lattice = reciprocal_lattice_from_outcar(args.outcar)

    pcar = procar.Procar.from_file(args.procar)
    effective_mass = pcar.effective_mass_calc(
        k_point_indices=args.k_points,
        band_index=args.band_index,
        reciprocal_lattice=reciprocal_lattice,
        spin=args.spin,
        printing=args.verbose,
    )
    print(effective_mass)


if __name__ == "__main__":
    main()
