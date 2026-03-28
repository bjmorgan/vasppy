#! /usr/bin/env python3

# Adapted from http://kitchingroup.cheme.cmu.edu/blog/2013/02/18/Nonlinear-curve-fitting/

import argparse
import warnings
import xml.etree.ElementTree as ET

import matplotlib  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
from numpy.typing import NDArray
from pymatgen.core import Structure
from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.outputs import UnconvergedVASPWarning
from scipy.optimize import leastsq

from vasppy.summary import find_vasp_calculations
from vasppy.utils import match_filename

matplotlib.use("agg")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Murnaghan fit script."""
    parser = argparse.ArgumentParser(
        description="Perform a Murnaghan equation of state fit across VASP subdirectories"
    )
    parser.add_argument(
        "-p", "--plot", action="store_true", help="generate murn.pdf plot of fit"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose output")
    args = parser.parse_args()
    return args


def read_vasprun(filename: str) -> Vasprun:
    """Read a vasprun.xml file with minimal parsing.

    Args:
        filename: Path to the vasprun.xml file.

    Returns:
        A Vasprun object.
    """
    return Vasprun(
        filename, parse_potcar_file=False, parse_dos=False, parse_eigen=False
    )


def read_data(verbose: bool = True) -> pd.DataFrame:
    """Read volume and energy data from VASP calculation subdirectories.

    Scans subdirectories for vasprun.xml and POSCAR files, collecting
    initial POSCAR volumes, relaxed volumes, and energies. Checks that
    the ratio of relaxed to initial volume is consistent across all
    calculations (indicating the same base cell geometry).

    Args:
        verbose: If True, print the collected data table.

    Returns:
        A DataFrame with columns: poscar_volume, volume, energy,
        converged, and volume_ratio.

    Raises:
        ValueError: If no VASP calculations are found, or if the
            volume ratios are inconsistent.
    """
    dir_list = find_vasp_calculations()
    if not dir_list:
        raise ValueError(
            "Did not find any subdirectories containing vasprun.xml or vasprun.xml.gz files"
        )
    data = []
    for d in dir_list:
        converged = True
        try:
            filename = match_filename(d + "vasprun.xml")
            if filename is None:
                continue
            with warnings.catch_warnings(record=True) as w:
                vasprun = read_vasprun(filename)
                for warning in w:
                    if isinstance(warning.message, UnconvergedVASPWarning):
                        converged = False
                    else:
                        print(warning.message)
        except (ET.ParseError, FileNotFoundError) as e:
            warnings.warn(f"Skipping {d}: {e}", stacklevel=2)
            continue
        poscar_structure = Structure.from_file(d + "POSCAR")
        data.append(
            [
                poscar_structure.volume,
                vasprun.final_structure.volume,
                vasprun.final_energy,
                converged,
            ]
        )
    column_titles = ["poscar_volume", "volume", "energy", "converged"]
    df = pd.DataFrame(data, columns=column_titles).sort_values(by="poscar_volume")
    df = df.reset_index(drop=True)
    df["volume_ratio"] = df.volume / df.poscar_volume
    volume_ratio_round = 4
    if verbose:
        print(df.to_string(index=False))
    if len(set(df.volume_ratio.round(volume_ratio_round))) != 1:
        raise ValueError("POSCAR volumes and relaxed volumes are inconsistent")
    return df


def murnaghan(
    vol: float | NDArray[np.floating],
    e0: float,
    b0: float,
    bp: float,
    v0: float,
) -> float | NDArray[np.floating]:
    """Calculate energy using the Murnaghan equation of state.

    Murnaghan, Proc. Nat. Acad. Sci. 30, 244 (1944).
    https://en.wikipedia.org/wiki/Murnaghan_equation_of_state
    cf. Fu and Ho, Phys. Rev. B 28, 5480 (1983).

    Args:
        vol: Volume(s) at which to evaluate the energy.
        e0: Energy at the equilibrium volume, E0.
        b0: Bulk modulus at the equilibrium volume, B0.
        bp: Pressure derivative of the bulk modulus, B0'.
        v0: Equilibrium volume, V0.

    Returns:
        The energy at the given volume(s).
    """
    energy = (
        e0 + b0 * vol / bp * (((v0 / vol) ** bp) / (bp - 1) + 1) - v0 * b0 / (bp - 1.0)
    )
    return energy


def objective(
    pars: tuple[float, ...],
    x: NDArray[np.floating],
    y: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Compute residuals between observed energies and the Murnaghan model.

    Args:
        pars: Murnaghan parameters (e0, b0, bp, v0).
        x: Observed volumes.
        y: Observed energies.

    Returns:
        Array of residuals (observed minus predicted).
    """
    err = y - murnaghan(x, *pars)
    return err


def lstsq_fit(
    volumes: NDArray[np.floating],
    energies: NDArray[np.floating],
) -> tuple[NDArray[np.floating], int]:
    """Fit the Murnaghan equation of state to volume-energy data.

    Args:
        volumes: Array of volumes.
        energies: Array of corresponding energies.

    Returns:
        A tuple of (fitted parameters, info flag) from
        scipy.optimize.leastsq.
    """
    e_min = energies.min()
    v_min = volumes[np.argwhere(energies == e_min)[0][0]]
    x0 = [e_min, 2.0, 10.0, v_min]  # initial guess of parameters
    plsq = leastsq(objective, x0, args=(volumes, energies))  # type: ignore[arg-type]
    return plsq


def make_plot(df: pd.DataFrame, fit_params: tuple[float, ...]) -> None:
    """Generate a plot of the equation of state fit.

    Args:
        df: DataFrame containing volume, energy, and converged columns.
        fit_params: Fitted Murnaghan parameters (e0, b0, bp, v0).
    """
    v_min = df.volume.min() * 0.99
    v_max = df.volume.max() * 1.01
    v_fitting = np.linspace(v_min, v_max, num=50)
    e_fitting = murnaghan(v_fitting, *fit_params)
    plt.figure(figsize=(8.0, 6.0))
    # plot converged data points
    loc = df.converged
    plt.plot(df[loc].volume, df[loc].energy, "o")
    # plot unconverged data points
    loc = [not b for b in df.converged]
    plt.plot(df[loc].volume, df[loc].energy, "o", c="grey")
    # plot fitted equation of state curve
    plt.plot(v_fitting, e_fitting, "--")
    plt.xlabel(r"volume [$\mathrm{\AA}^3$]")
    plt.ylabel(r"energy [eV]")
    plt.tight_layout()
    plt.savefig("murn.pdf")


def fit(verbose: bool = False, plot: bool = False) -> None:
    """Perform a Murnaghan equation of state fit and print results.

    Args:
        verbose: If True, print the raw data table.
        plot: If True, save a plot to murn.pdf.
    """
    df = read_data(verbose=verbose)
    e0, b0, bp, v0 = lstsq_fit(np.array(df.volume), np.array(df.energy))[0]
    if plot:
        make_plot(df, (e0, b0, bp, v0))
    print(f"E0: {e0:.4f}")
    print(f"V0: {v0:.4f}")
    print(f"opt. POSCAR volume: {v0 / df.volume_ratio.mean():.4f}")


def main() -> None:
    """Entry point for the murnfit command-line script."""
    args = parse_args()
    fit(verbose=args.verbose, plot=args.plot)


if __name__ == "__main__":
    main()
