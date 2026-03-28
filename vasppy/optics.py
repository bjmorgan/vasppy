"""Functions for working with optical properties from vasprun.xml."""

from math import pi
from typing import cast

import numpy as np
from scipy.constants import physical_constants, speed_of_light  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.axes import Axes
from matplotlib.figure import Figure

eV_to_recip_cm = 1.0 / (
    physical_constants["Planck constant in eV s"][0] * speed_of_light * 1e2
)


def matrix_eigvals(matrix: np.ndarray) -> np.ndarray:
    """Calculate the eigenvalues of a matrix.

    Args:
        matrix: The matrix to diagonalise.

    Returns:
        Array of the matrix eigenvalues.
    """
    return np.linalg.eigvals(matrix)


def to_matrix(
    xx: float,
    yy: float,
    zz: float,
    xy: float,
    yz: float,
    xz: float,
) -> np.ndarray:
    """Convert a list of matrix components to a symmetric 3x3 matrix.

    Inputs should be in the order xx, yy, zz, xy, yz, xz.

    Args:
        xx: xx component of the matrix.
        yy: yy component of the matrix.
        zz: zz component of the matrix.
        xy: xy component of the matrix.
        yz: yz component of the matrix.
        xz: xz component of the matrix.

    Returns:
        The matrix as a 3x3 numpy array.
    """
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def plot_dielectric_functions(
    dielectric: list,
    ax: Axes | None = None,
) -> Figure | None:
    """Plot the real and imaginary dielectric functions.

    Args:
        dielectric: Dielectric data in pymatgen vasprun format.
            Element 0 is a list of energies; element 1 contains the real
            dielectric tensors; element 2 contains the imaginary dielectric
            tensors.
        ax: Optional matplotlib Axes to plot on. If ``None``, a new figure and
            axes are created.

    Returns:
        The matplotlib Figure if a new one was created, otherwise ``None``.
    """
    real_dielectric = parse_dielectric_data(dielectric[1])
    imag_dielectric = parse_dielectric_data(dielectric[2])
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6.0, 3.0))
    else:
        fig = None
    ax.plot(
        dielectric[0], np.mean(real_dielectric, axis=1), "-", zorder=2
    )  # better to pass in v.dielectric
    ax.plot(dielectric[0], np.mean(imag_dielectric, axis=1), "-", zorder=2)
    ax.set_xlim((0, 8))
    ax.set_ylim((0, 5))
    return fig


def parse_dielectric_data(data: list) -> np.ndarray:
    """Convert a set of 2D vasprun-formatted dielectric data to eigenvalues.

    Converts each entry to a 3x3 symmetric numpy matrix and returns the
    eigenvalues of each matrix.

    Args:
        data: Length-N list of dielectric data. Each entry should be a list of
            ``[xx, yy, zz, xy, yz, xz]`` dielectric tensor elements.

    Returns:
        An Nx3 numpy array. Each row contains the eigenvalues for the
        corresponding row in ``data``.
    """
    return np.array([matrix_eigvals(to_matrix(*e)) for e in data])


def absorption_coefficient(dielectric: list) -> np.ndarray:
    """Calculate the optical absorption coefficient from dielectric data.

    Args:
        dielectric: A list containing the dielectric response function in
            pymatgen vasprun format.

            - Element 0: list of energies.
            - Element 1: real dielectric tensors in ``[xx, yy, zz, xy, yz, xz]``
              format.
            - Element 2: imaginary dielectric tensors in
              ``[xx, yy, zz, xy, yz, xz]`` format.

    Returns:
        Absorption coefficient using eV as frequency units (cm :sup:`-1`).

    Note:
        The absorption coefficient is calculated as

        .. math::

            \\alpha = \\frac{2\\sqrt{2}\\pi}{\\lambda}
            \\sqrt{-\\epsilon_1 + \\sqrt{\\epsilon_1^2 + \\epsilon_2^2}}
    """
    energies_in_eV = np.array(dielectric[0])
    real_dielectric = parse_dielectric_data(dielectric[1])
    imag_dielectric = parse_dielectric_data(dielectric[2])
    epsilon_1 = np.mean(real_dielectric, axis=1)
    epsilon_2 = np.mean(imag_dielectric, axis=1)
    return cast(
        np.ndarray,
        2.0
        * np.sqrt(2.0)
        * pi
        * eV_to_recip_cm
        * energies_in_eV
        * np.sqrt(-epsilon_1 + np.sqrt(epsilon_1**2 + epsilon_2**2)),
    )
