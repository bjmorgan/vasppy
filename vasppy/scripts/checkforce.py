#! /usr/bin/env python3

import argparse
import re
import numpy as np
from dataclasses import dataclass, field

from vasppy.outcar import forces_from_outcar


@dataclass
class ForcesData:
    forces: np.ndarray = field(repr=False)
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

    def print_forces(self) -> None:
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Checks for convergence of geometry optimisations in VASP"
    )
    parser.add_argument(
        "-o",
        "--outcar",
        type=str,
        default="OUTCAR",
        help='The filepath of the OUTCAR file to be parsed. Default is "OUTCAR")',
    )
    parser.add_argument(
        "-c",
        "--convergence",
        type=float,
        help="Set force convergence. Default is to read the convergence from the OUTCAR file.",
    )
    # Create mutually exclusive group for --verbose and --all
    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose output. Show convergence status for all atoms.",
    )
    output_group.add_argument(
        "-a",
        "--all",
        action="store_true",
        help="Print summary data for every ionic step.",
    )
    args = parser.parse_args()
    return args


def get_forces_data(
    outcar_filename: str = "OUTCAR",
    convergence: float | None = None,
) -> ForcesData:
    """Parse an OUTCAR file and return forces data, including various summary statistics.

    Args:
        outcar_filename: OUTCAR filename. Default is "OUTCAR".
        convergence: Optionally set the forces convergence threshold.
            Default is to read from the OUTCAR file.

    Returns:
        ForcesData object with forces and summary statistics.

    """
    if not convergence:
        convergence = read_ediffg_from_outcar(outcar_filename)
    forces = forces_from_outcar(outcar_filename, last_one_only=True)
    forces_data = ForcesData(forces=forces, convergence=convergence)
    return forces_data

def get_all_forces_data(
    outcar_filename: str = "OUTCAR",
    convergence: float | None = None,
):
    """Parse an OUTCAR file and yield forces data for all ionic steps.

    Args:
        outcar_filename: OUTCAR filename. Default is "OUTCAR".
        convergence: Optionally set the forces convergence threshold.
            Default is to read from the OUTCAR file.

    Yields:
        ForcesData object for each ionic step.

    """
    if not convergence:
        convergence = read_ediffg_from_outcar(outcar_filename)
    all_forces = forces_from_outcar(outcar_filename, last_one_only=False)
    for step_forces in all_forces:
        yield ForcesData(forces=step_forces, convergence=convergence)

def read_ediffg_from_outcar(outcar_filename: str) -> float:
    """Read EDIFFG convergence criterion from OUTCAR file.

    Args:
        outcar_filename: OUTCAR filename.

    Returns:
        Convergence threshold (absolute value of EDIFFG).

    Raises:
        ValueError: If EDIFFG cannot be found in the file.

    """
    with open(outcar_filename, "r") as f:
        outcar = f.read()
    convergence_re = re.compile(r"EDIFFG = -([\d\.E-]+)")
    try:
        convergence = float(convergence_re.findall(outcar)[0])
    except IndexError as exc:
        raise ValueError("Unable to read EDIFFG.") from exc
    return convergence

def main() -> None:
    args = parse_args()

    if args.all:
        # Print summary for all ionic steps
        print(f"{'Max Force':<12}{'Non-opt':<10}{'Mean Excess':<10}")
        for step, forces_data in enumerate(
            get_all_forces_data(
                outcar_filename=args.outcar, convergence=args.convergence
            ),
            start=1,
        ):
            print(
                f"{forces_data.max_force:<12.6f}"
                f"{forces_data.n_not_converged:<10}"
                f"{forces_data.mean_excess_force:<10.6f}"
            )
    else:
        # Original behavior: print only last step
        forces_data = get_forces_data(
            outcar_filename=args.outcar, convergence=args.convergence
        )
        if args.verbose:
            forces_data.print_forces()
            print()
        print(f"remainder:  {forces_data.mean_excess_force:.6f}")
        print(f"maximum:    {forces_data.max_force:.6f}")
        print(f"non-opt:    {forces_data.n_not_converged} / {forces_data.number_of_ions}")

if __name__ == "__main__":
    main()
