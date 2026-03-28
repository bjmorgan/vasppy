#! /usr/bin/env python3
"""Check species consistency between a VASP POSCAR and POTCAR pair.

Species are considered consistent if the species labels in the POSCAR match
the start of the pseudopotential labels in the POTCAR, in order. For example,
a POSCAR containing ``Na Cl`` will match a POTCAR containing ``Na_pv Cl``.
If any species labels do not match, an ``AttributeError`` is raised.

The ``-p`` flag additionally checks that all pseudopotentials belong to a
specific pseudopotential set.
"""

import sys
if sys.platform != "win32":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

import argparse

from pymatgen.core import Structure

from vasppy.summary import potcar_spec, potcar_sets


def unique_species_from_structure(filename: str) -> list[str]:
    """Return unique species labels from a POSCAR file, in order of appearance.

    Args:
        filename: Path to the VASP POSCAR file.

    Returns:
        A list of element symbol strings, in the order they first appear.
    """
    structure = Structure.from_file(filename)
    return list(dict.fromkeys(site.species_string for site in structure))


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command-line arguments for check_species.

    Returns:
        Parsed argument namespace with ``poscar``, ``potcar``, and
        optional ``ppset`` attributes.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Check species consistency between a VASP POSCAR file and a POTCAR file."
        )
    )
    parser.add_argument(
        "poscar",
        help="filename of the VASP POSCAR to be processed",
        nargs="?",
        default="POSCAR",
    )
    parser.add_argument(
        "potcar",
        help="filename of the VASP POTCAR to be processed",
        nargs="?",
        default="POTCAR",
    )
    parser.add_argument(
        "-p",
        "--ppset",
        help=(
            "check whether the POTCAR pseudopotentials belong to a specific "
            "pseudopotential set"
        ),
        choices=potcar_sets,
    )
    return parser.parse_args()


def main() -> None:
    """Entry point: check that POSCAR and POTCAR species are consistent."""
    args = parse_command_line_arguments()
    species = unique_species_from_structure(args.poscar)
    potcar_names, potcar_datasets = potcar_spec(args.potcar)
    for i, (sp, name, dataset) in enumerate(
        zip(species, potcar_names, potcar_datasets, strict=True), 1
    ):
        if not name.startswith(sp):
            raise AttributeError(
                f"Species {i} mismatch:\nPOSCAR contains {sp}\nPOTCAR contains {name}"
            )
        if args.ppset and args.ppset != dataset:
            raise AttributeError(f"Pseudopotential set mismatch: {potcar_datasets}")


if __name__ == "__main__":
    main()
