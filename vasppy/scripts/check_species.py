#! /usr/bin/env python3

from vasppy.poscar import Poscar
from vasppy.summary import potcar_spec, potcar_sets
import argparse

"""
A command line utility for testing species consistency between a VASP POSCAR and POTCAR pair of files. Species are considered consistent if the species labels used in the POSCAR file match the start of the pseudopotential labels in the POTCAR file, in order. e.g. a POSCAR that contains `Ti O` will match a POTCAR that contains `Ti_pv O`. If any species labels do not match the script raises an AttributeError.

The `-p` flag will check that all the pseudopotentials in the POTCAR file belong to a specific pseudopotential set.
"""


def parse_command_line_arguments():
    parser = argparse.ArgumentParser(
        description="Check species consistency between a VASP POSCAR file and a POTCAR file."
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
        help="check whether the POTCAR pseudopotentials belong to a specific pseudopotential set",
        choices=potcar_sets,
    )
    return parser.parse_args()


def main():
    args = parse_command_line_arguments()
    poscar = Poscar.from_file(args.poscar)
    potcar_names, potcar_datasets = potcar_spec(args.potcar)
    for i, (species, name, dataset) in enumerate(
        zip(poscar.atoms, potcar_names, potcar_datasets, strict=True), 1,
    ):
        if not name.startswith(species):
            raise AttributeError(
                "Species {} mismatch:\nPOSCAR contains {}\nPOTCAR contains {}".format(
                    i, species, name
                )
            )
        if args.ppset and args.ppset != dataset:
            raise AttributeError(
                "Pseudopotential set mismatch: {}".format(potcar_datasets)
            )


if __name__ == "__main__":
    main()
