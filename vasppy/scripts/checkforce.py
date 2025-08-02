#! /usr/bin/env python3

import argparse
import re
import numpy as np
from dataclasses import dataclass, field
from typing import Union

forces_re = re.compile(
    r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)"
)


@dataclass
class ForcesData:
    forces: Union[list, np.ndarray] = field(repr=False)
    convergence: float
    number_of_ions: int = field(init=False)
    mean_excess_force: float = field(init=False)
    max_force: float = field(init=False)
    n_not_converged: int = field(init=False)

    def __post_init__(self):
        if isinstance(self.forces, list):
            self.forces = np.array(self.forces)
        self.number_of_ions = self.forces.shape[0]
        force_magnitude = np.sqrt(np.sum(self.forces**2, axis=1))
        self.mean_excess_force = (
            np.where(
                force_magnitude > self.convergence,
                force_magnitude - self.convergence,
                0,
            ).sum()
            / self.number_of_ions
        )
        self.n_not_converged = np.sum(force_magnitude > self.convergence)
        self.max_force = force_magnitude.max()

    def print_forces(self):
        for theseforces in self.forces:
            output_str = [f"  {force: .6f}" for force in theseforces]
            for force in theseforces:
                if abs(force) > self.convergence:
                    output_str.append("   ****")
                else:
                    output_str.append("   conv")
            force_norm = np.linalg.norm(theseforces)
            output_str.append(f"   {force_norm:.6f}")
            if np.linalg.norm(force_norm) > self.convergence:
                output_str.append("   ****")
            else:
                output_str.append("   conv")
            print("".join(output_str))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Checks for convergence of geometry optimisations in VASP"
    )
    parser.add_argument(
        "-o",
        "--outcar",
        type=str,
        default="OUTCAR",
        help='The filepath of the OUTCAR file to be parsed. Defaule is "OUTCAR")',
    )
    parser.add_argument(
        "-c",
        "--convergence",
        type=float,
        help="Set force convergence. Default is to read the convergence from the OUTCAR file.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose output. Show convergence status for all atoms.",
    )
    parser.add_argument("-w", "--warn", action="store_true", help="Print warnings.")
    args = parser.parse_args()
    return args


def get_forces_data(outcar_filename="OUTCAR", convergence=None, warn=False):
    """Parse an OUTCAR file and return forces data, includig various summary statistics.

    args:
        outcar_filename (optional, `str`): OUTCAR filename. Default is "OUTCAR".
        convergence (optional, `str`): Optionally set the forces convergence threshold.
            Default is to read from the OUTCAR file.
        warn (optional, `bool`): Optionally warn if the OUTCAR file appears to be incomplete.

    returns:
        ForcesData

    """
    with open(outcar_filename, "r") as f:
        outcar = f.read()

    number_of_ions_re = re.compile(r"NIONS = \s+([\d]+)")
    try:
        number_of_ions = int(number_of_ions_re.findall(outcar)[0])
    except IndexError as exc:
        raise ValueError("Unable to read NIONS.") from exc
    except:
        raise
    if not convergence:
        convergence_re = re.compile(r"EDIFFG = -([\d\.E-]+)")
        try:
            convergence = float(convergence_re.findall(outcar)[0])
        except IndexError as exc:
            raise ValueError("Unable to read EDIFFG.") from exc
        except:
            raise

    # find force output block positions
    forces_block_start = []
    forces_header_re = re.compile(r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)")
    for i, line in enumerate(outcar.split("\n")):
        if forces_header_re.search(line):
            forces_block_start.append(i)
    if not forces_block_start:
        raise ValueError("No FORCES blocks found.")
    # Ideally we want the most recent complete block of forces data.
    most_recent_block_start = forces_block_start[-1] + 2
    forces_lines = outcar.split("\n")[
        most_recent_block_start : most_recent_block_start + number_of_ions
    ]
    if not forces_block_is_well_formed(forces_lines):
        # If the most recent forces block is ill-formed, try to parse the previous block.
        if warn:
            print(
                "The last FORCES block is not well-formed. Trying to parse the preceeding block."
            )
        next_most_recent_block_start = forces_block_start[-2] + 2
        forces_lines = outcar.split("\n")[
            next_most_recent_block_start : next_most_recent_block_start + number_of_ions
        ]
        if not forces_block_is_well_formed(forces_lines):
            # If the last two forces blocks are ill-formed, we assume the input file is mis-formatted.
            raise Exception(
                "The last two FORCES blocks are not well-formed. Your input might be mis-formatted."
            )
    forces = []
    for line in forces_lines:
        forces.append([float(s) for s in line.split()[-3:]])
    forces_data = ForcesData(forces=forces, convergence=convergence)
    return forces_data


def forces_block_is_well_formed(forces_lines):
    return all(forces_re.match(line) for line in forces_lines)


def main():
    args = parse_args()
    forces_data = get_forces_data(
        outcar_filename=args.outcar, convergence=args.convergence, warn=args.warn
    )
    if args.verbose:
        forces_data.print_forces()
        print()
    print(f"remainder:  {forces_data.mean_excess_force:.6f}")
    print(f"maximum:    {forces_data.max_force:.6f}")
    print(f"non-opt:    {forces_data.n_not_converged} / {forces_data.number_of_ions}")


if __name__ == "__main__":
    main()
