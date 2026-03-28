#! /usr/bin/env python3

"""Manipulate VASP POSCAR files.

Reads a POSCAR file and outputs it with optional transformations
such as supercell generation, Bohr conversion, and coordinate
type selection.
"""

from __future__ import annotations

import sys
if sys.platform != "win32":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

import argparse

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar as PmgPoscar

from vasppy.units import angstrom_to_bohr


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command-line arguments for POSCAR manipulation.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Manipulates VASP POSCAR files")
    parser.add_argument("poscar", help="filename of the VASP POSCAR to be processed")
    parser.add_argument(
        "-l",
        "--label",
        type=int,
        choices=[1, 4],
        help="label coordinates with atom name at position {1,4}",
    )
    parser.add_argument(
        "-c", "--coordinates-only", help="only output coordinates", action="store_true"
    )
    parser.add_argument(
        "-t",
        "--coordinate-type",
        type=str,
        choices=["c", "cartesian", "d", "direct"],
        default="direct",
        help="specify coordinate type for output {(c)artesian|(d)irect} [default = (d)irect]",
    )
    parser.add_argument(
        "-g", "--group", help="group atoms within supercell", action="store_true"
    )
    parser.add_argument(
        "-s",
        "--supercell",
        type=int,
        nargs=3,
        metavar=("h", "k", "l"),
        help="construct supercell by replicating (h,k,l) times along [a b c]",
    )
    parser.add_argument(
        "-b",
        "--bohr",
        action="store_true",
        help="assumes the input file is in Angstrom, and converts everything to bohr",
    )
    parser.add_argument(
        "-n",
        "--number-atoms",
        action="store_true",
        help="label coordinates with atom number",
    )
    parser.add_argument(
        "-o",
        "--orthorhombic",
        action="store_true",
        help="force orthorhombic cell matrix (set off-diagonal elements to zero)",
    )
    parser.add_argument(
        "--scale",
        action="store_true",
        help="scale the lattice parameters by the scaling factor",
    )
    parser.add_argument(
        "--selective",
        choices=["T", "F"],
        help="generate Selective Dynamics POSCAR with all values set to T / F",
    )
    args = parser.parse_args()
    return args


def _species_groups(structure: Structure) -> tuple[list[str], list[int]]:
    """Return species names and counts grouped by contiguous runs.

    Args:
        structure: A pymatgen Structure.

    Returns:
        Tuple of (species names, atom counts) matching POSCAR format.
    """
    from itertools import groupby

    labels = [site.species_string for site in structure]
    species = []
    counts = []
    for key, group in groupby(labels):
        species.append(key)
        counts.append(len(list(group)))
    return species, counts


def _get_coordinates(structure: Structure, coordinate_type: str) -> np.ndarray:
    """Return coordinates for the given coordinate type.

    Args:
        structure: A pymatgen Structure.
        coordinate_type: Either 'Direct' or 'Cartesian'.

    Returns:
        Array of coordinates with shape (n_sites, 3).
    """
    return structure.frac_coords if coordinate_type == "Direct" else structure.cart_coords


def _output_header(
    structure: Structure,
    title: str,
    scaling: float,
    coordinate_type: str,
    opts: dict[str, object],
) -> None:
    """Print the POSCAR header to stdout.

    Args:
        structure: A pymatgen Structure.
        title: Title line for the POSCAR.
        scaling: Scaling factor for the lattice.
        coordinate_type: Either 'Direct' or 'Cartesian'.
        opts: Output options dictionary.
    """
    print(title)
    print(scaling)
    matrix = np.array(structure.lattice.matrix)
    if opts.get("orthorhombic"):
        matrix = matrix * np.eye(3)
    for row in matrix:
        print("".join(f"   {v: .10f}" for v in row))
    species, counts = _species_groups(structure)
    print(" ".join(species))
    print(" ".join(str(n) for n in counts))
    if opts.get("selective"):
        print("Selective Dynamics")
    print(coordinate_type)


def _output_coordinates(
    structure: Structure,
    coordinate_type: str,
    opts: dict[str, object],
) -> None:
    """Print coordinate lines to stdout.

    Args:
        structure: A pymatgen Structure.
        coordinate_type: Either 'Direct' or 'Cartesian'.
        opts: Output options dictionary.
    """
    coords = _get_coordinates(structure, coordinate_type)
    labels = [site.species_string for site in structure]
    for i, (coord, label) in enumerate(zip(coords, labels)):
        prefix = label.ljust(6) if opts.get("label") == 1 else ""
        coord_str = "".join(f"  {v: .10f}" for v in coord)
        suffix_parts: list[str] = []
        if opts.get("selective"):
            suffix_parts.append(f" {opts['selective']} {opts['selective']} {opts['selective']}")
        if opts.get("numbered"):
            suffix_parts.append(f" {i + 1}")
        if opts.get("label") == 4:
            suffix_parts.append(f" {label}")
        print(f"{prefix}{coord_str}{''.join(suffix_parts)}")


def _convert_to_bohr(structure: Structure) -> tuple[Structure, float]:
    """Convert a structure from Angstrom to Bohr units.

    Args:
        structure: A pymatgen Structure in Angstrom.

    Returns:
        Tuple of (new Structure with lattice in Bohr, scaling factor).
    """
    new_lattice = Lattice(structure.lattice.matrix / angstrom_to_bohr)
    return Structure(new_lattice, structure.species, structure.frac_coords), angstrom_to_bohr


def main() -> None:
    """Main entry point for the proc_poscar script."""
    args = parse_command_line_arguments()
    coordinate_type = "Cartesian" if args.coordinate_type[0].lower() == "c" else "Direct"

    poscar_data = PmgPoscar.from_file(args.poscar)
    structure = poscar_data.structure
    title = poscar_data.comment
    # pymatgen already applies the scaling factor when reading,
    # so we track it separately for output formatting.
    scaling = 1.0

    if args.supercell:
        if args.group:
            for i in args.supercell:
                if i % 2 == 1 and i > 1:
                    raise Exception(
                        "odd supercell expansions != 1 are incompatible with automatic grouping"
                    )
        structure.make_supercell(args.supercell)

    if args.bohr:
        structure, scaling = _convert_to_bohr(structure)

    output_opts: dict[str, object] = {
        "label": args.label,
        "numbered": args.number_atoms,
        "coordinates_only": args.coordinates_only,
        "selective": args.selective,
        "orthorhombic": args.orthorhombic,
    }
    if not args.coordinates_only:
        _output_header(
            structure,
            title=title,
            scaling=scaling,
            coordinate_type=coordinate_type,
            opts=output_opts,
        )
    _output_coordinates(structure, coordinate_type=coordinate_type, opts=output_opts)


if __name__ == "__main__":
    main()
