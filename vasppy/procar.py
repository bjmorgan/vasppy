from functools import reduce
from copy import deepcopy
from typing import cast
import math
import re

import numpy as np
from .units import angstrom_to_bohr, ev_to_hartree
from .band import Band


class KPoint:
    """Represents a single k-point from a VASP PROCAR file.

    Attributes:
        index: 1-based k-point index.
        frac_coords: Fractional reciprocal coordinates as a 1D numpy array.
        weight: Integration weight of this k-point.
    """

    def __init__(
        self,
        index: int,
        frac_coords: np.ndarray,
        weight: float,
    ) -> None:
        """Initialise a KPoint object.

        Args:
            index: 1-based k-point index.
            frac_coords: Fractional reciprocal coordinates, shape (3,).
            weight: Integration weight of this k-point.
        """
        self.index = index
        self.frac_coords = frac_coords
        self.weight = weight

    def cart_coords(self, reciprocal_lattice: np.ndarray) -> np.ndarray:
        """Convert the reciprocal fractional coordinates for this k-point to
        reciprocal Cartesian coordinates.

        Args:
            reciprocal_lattice: 3x3 numpy array containing the Cartesian
                reciprocal lattice.

        Returns:
            The reciprocal Cartesian coordinates of this k-point.
        """
        return cast(np.ndarray, np.dot(self.frac_coords, reciprocal_lattice))

    def __eq__(self, other: object) -> bool:
        """Check equality based on index, fractional coordinates, and weight."""
        if not isinstance(other, KPoint):
            return NotImplemented
        return (
            (self.index == other.index)
            and (self.frac_coords == other.frac_coords).all()
            and (self.weight == other.weight)
        )

    def __repr__(self) -> str:
        coords_str = " ".join(str(c) for c in self.frac_coords)
        return f"k-point {self.index}: {coords_str} weight = {self.weight}"


def get_numbers_from_string(string: str) -> list[float]:
    """Extract all numbers from a string.

    Args:
        string: Input string containing numbers.

    Returns:
        List of floats parsed from the string.
    """
    p = re.compile(r"-?\d+[.\d]*")
    return [float(s) for s in p.findall(string)]


def k_point_parser(string: str) -> list[KPoint]:
    """Parse k-point data from a PROCAR string.

    Finds all lines of the form::

         k-point    1 :    0.50000000 0.25000000 0.75000000     weight = 0.00806452

    and extracts the k-point index, reciprocal fractional coordinates, and weight
    into a :obj:`procar.KPoint` object.

    Args:
        string: String containing a full PROCAR file.

    Returns:
        List of :obj:`procar.KPoint` objects.
    """
    regex = re.compile(
        r"k-point\s+(\d+)\s*:\s+([- ][01].\d+)([- ][01].\d+)([- ][01].\d+)\s+weight = *(-*[01].\d+)"
    )
    captured = regex.findall(string)
    k_points = []
    for kp in captured:
        index = int(kp[0])
        frac_coords = np.array([float(s) for s in kp[1:4]])
        weight = float(kp[4])
        k_points.append(KPoint(index=index, frac_coords=frac_coords, weight=weight))
    return k_points


def projections_parser(string: str) -> np.ndarray:
    """Parse projection data from a PROCAR string.

    Args:
        string: String containing PROCAR file data.

    Returns:
        2D numpy array of projection values.
    """
    regex = re.compile(r"([-.\d\se]+tot.+)\n")
    data = regex.findall(string)
    data = [x.replace("tot", "0") for x in data]
    return np.array([x.split() for x in data], dtype=float)


def area_of_a_triangle_in_cartesian_space(
    a: np.ndarray, b: np.ndarray, c: np.ndarray
) -> float:
    """Return the area of a triangle defined by three points in Cartesian space.

    Args:
        a: Cartesian coordinates of point A.
        b: Cartesian coordinates of point B.
        c: Cartesian coordinates of point C.

    Returns:
        The area of the triangle.
    """
    return float(0.5 * np.linalg.norm(np.cross(b - a, c - a)))


def points_are_in_a_straight_line(
    points: np.ndarray | list[np.ndarray],
    tolerance: float = 1e-7,
) -> bool:
    """Check whether a set of points fall on a straight line.

    Calculates the areas of triangles formed by triplets of the points.
    Returns False if any of these areas are larger than the tolerance.

    Args:
        points: List of Cartesian coordinates for each point.
        tolerance: Maximum triangle size for points to be considered
            collinear. Default is 1e-7.

    Returns:
        True if all points fall on a straight line (within the allowed
        tolerance).
    """
    a = points[0]
    b = points[1]
    for c in points[2:]:
        if area_of_a_triangle_in_cartesian_space(a, b, c) > tolerance:
            return False
    return True


def two_point_effective_mass(
    cartesian_k_points: np.ndarray,
    eigenvalues: np.ndarray
) -> float:
    """Calculate the effective mass given eigenvalues at two k-points.

    Reimplemented from Aron Walsh's original effective mass Fortran code.

    Args:
        cartesian_k_points: 2D numpy array containing the k-points in
            (reciprocal) Cartesian coordinates.
        eigenvalues: numpy array containing the eigenvalues at each k-point.

    Returns:
        The effective mass.

    Raises:
        ValueError: If ``cartesian_k_points`` does not contain exactly
            two k-points, or if ``eigenvalues`` does not have exactly
            two elements.
    """
    if cartesian_k_points.shape[0] != 2:
        raise ValueError(
            f"two_point_effective_mass requires exactly 2 k-points, got {cartesian_k_points.shape[0]}"
        )
    if eigenvalues.size != 2:
        raise ValueError(
            f"two_point_effective_mass requires exactly 2 eigenvalues, got {eigenvalues.size}"
        )
    dk = cartesian_k_points[1] - cartesian_k_points[0]
    mod_dk = np.sqrt(np.dot(dk, dk))
    delta_e = (eigenvalues[1] - eigenvalues[0]) * ev_to_hartree * 2.0
    effective_mass = mod_dk * mod_dk / delta_e
    return float(effective_mass)


def least_squares_effective_mass(
    cartesian_k_points: np.ndarray,
    eigenvalues: np.ndarray
) -> float:
    """Calculate the effective mass using a least squares quadratic fit.

    Args:
        cartesian_k_points: Cartesian reciprocal coordinates for the k-points.
        eigenvalues: Energy eigenvalues at each k-point to be used in the fit.

    Returns:
        The fitted effective mass.

    Raises:
        ValueError: If the k-points do not sit on a straight line.
    """
    if not points_are_in_a_straight_line(cartesian_k_points):
        raise ValueError("k-points are not collinear")
    dk = cartesian_k_points - cartesian_k_points[0]
    mod_dk = np.linalg.norm(dk, axis=1)
    effective_mass = 1.0 / (np.polyfit(mod_dk, eigenvalues, 2)[0] * ev_to_hartree * 2.0)
    return float(effective_mass)


class Procar:
    """Object for working with VASP PROCAR data.

    Attributes:
        data: A 5D numpy array that stores the projection data.
            Axes are k-points, bands, spin-channels, ions and sum over
            ions, lm-projections.
        bands: A numpy array of :obj:`Band` objects containing band index,
            energy, and occupancy data.
        k_points: A numpy array of :obj:`KPoint` objects containing
            fractional coordinates and weights for each k-point.
        number_of_k_points: The number of k-points.
        number_of_bands: The number of bands.
        spin_channels: Number of spin channels in the PROCAR data:

            - 1 for non-spin-polarised calculations.
            - 2 for spin-polarised calculations.
            - 4 for non-collinear calculations.

        number_of_ions: The number of ions.
        number_of_projections: The number of projections.
        calculation: Dictionary of True/False values describing the
            calculation type. Keys are ``'non_spin_polarised'``,
            ``'non_collinear'``, and ``'spin_polarised'``.
    """

    def __init__(
        self,
        data: np.ndarray,
        bands: np.ndarray,
        k_points: list[KPoint],
        number_of_k_points: int,
        number_of_bands: int,
        number_of_ions: int,
        number_of_projections: int,
        spin_channels: int,
        k_point_blocks: int,
        calculation: dict[str, bool],
        negative_occupancies: str = "warn",
    ) -> None:
        """Initialise a Procar object with fully parsed data.

        Args:
            data: 5D numpy array of projection data with axes
                (k-points, bands, spin-channels, ions+tot, projections).
            bands: 1D numpy array of :obj:`Band` objects.
            k_points: List of :obj:`KPoint` objects.
            number_of_k_points: Number of k-points.
            number_of_bands: Number of bands.
            number_of_ions: Number of ions.
            number_of_projections: Number of lm-projections.
            spin_channels: Number of spin channels (1, 2, or 4).
            k_point_blocks: Number of k-point blocks (1 or 2).
            calculation: Dictionary describing the calculation type with
                keys ``'non_spin_polarised'``, ``'non_collinear'``, and
                ``'spin_polarised'``.
            negative_occupancies: How to handle negative occupancies.
                Accepted values are ``'warn'``, ``'raise'``, or ``'zero'``.

        Raises:
            ValueError: If ``negative_occupancies`` is not one of the
                accepted values.
        """
        if negative_occupancies not in ["warn", "raise", "zero"]:
            raise ValueError(
                "negative_occupancies can be one of [ 'warn', 'raise', 'zero' ]"
            )
        self._data = data
        self._bands = bands
        self._k_points = k_points
        self._number_of_k_points = number_of_k_points
        self._number_of_bands = number_of_bands
        self._number_of_ions = number_of_ions
        self._number_of_projections = number_of_projections
        self._spin_channels = spin_channels
        self._k_point_blocks = k_point_blocks
        self.calculation = calculation
        self.negative_occupancies = negative_occupancies

    @property
    def occupancy(self) -> np.ndarray:
        """Band index and occupancy for all bands.

        Returns:
            2D array with columns [band_index, occupancy].
        """
        return np.array([[band.index, band.occupancy] for band in self._bands])

    def __add__(self, other: "Procar") -> "Procar":
        """Concatenate two Procar objects along the k-point axis.

        Args:
            other: Another :obj:`Procar` instance to concatenate with.

        Returns:
            A new :obj:`Procar` containing data from both objects.

        Raises:
            ValueError: If the two objects are incompatible (mismatched
                spin channels, number of ions, bands, projections,
                k-point blocks, or calculation type).
        """
        if self.spin_channels != other.spin_channels:
            raise ValueError(
                f"Can only concatenate Procars with equal spin_channels: "
                f"{self.spin_channels}, {other.spin_channels}"
            )
        if self.number_of_ions != other.number_of_ions:
            raise ValueError(
                f"Can only concatenate Procars with equal number_of_ions: "
                f"{self.number_of_ions}, {other.number_of_ions}"
            )
        if self.number_of_bands != other.number_of_bands:
            raise ValueError(
                f"Can only concatenate Procars with equal number_of_bands: "
                f"{self.number_of_bands}, {other.number_of_bands}"
            )
        if self.number_of_projections != other.number_of_projections:
            raise ValueError(
                f"Can only concatenate Procars with equal number_of_projections: "
                f"{self.number_of_projections}, {other.number_of_projections}"
            )
        if self._k_point_blocks != other._k_point_blocks:
            raise ValueError(
                f"Can only concatenate Procars with equal k_point_blocks: "
                f"{self._k_point_blocks}, {other._k_point_blocks}"
            )
        if self.calculation != other.calculation:
            raise ValueError(
                f"Can only concatenate Procars from equal calculations: "
                f"{self.calculation}, {other.calculation}"
            )
        new_procar = deepcopy(self)
        new_procar._data = np.concatenate((self._data, other._data), axis=0)
        new_procar._number_of_k_points = (
            self.number_of_k_points + other.number_of_k_points
        )
        new_procar._bands = np.ravel(np.concatenate([self.bands, other.bands], axis=1))
        new_procar._k_points = self._k_points + other._k_points
        for i, kp in enumerate(new_procar._k_points, 1):
            kp.index = i
        new_procar.sanity_check()
        return new_procar

    def sanity_check(self) -> None:
        """Verify that the parsed data is internally consistent.

        Raises:
            ValueError: If the number of k-points or bands in the header
                does not match what was found in the file body.
        """
        if self._number_of_k_points != len(self._k_points):
            raise ValueError(
                f"k-point number mismatch: {self._number_of_k_points} in header; "
                f"{len(self._k_points)} in file"
            )
        read_bands = len(self._bands) / self._number_of_k_points / self._k_point_blocks
        if self._number_of_bands != read_bands:
            raise ValueError(
                f"band mismatch: {self._number_of_bands} in header; {read_bands} in file"
            )

    @classmethod
    def from_files(cls, filenames: list[str], **kwargs) -> "Procar":
        """Create a :obj:`Procar` object from a series of VASP ``PROCAR`` files.

        Useful when a band-structure calculation has been split over multiple
        VASP calculations, for example when using hybrid functionals.

        Args:
            filenames: List of ``PROCAR`` filenames to read and concatenate.
            **kwargs: See :meth:`from_file` for a description of keyword
                arguments.

        Returns:
            A combined :obj:`Procar` instance.
        """
        pcars = [cls.from_file(f, **kwargs) for f in filenames]
        return reduce(cls.__add__, pcars)

    @classmethod
    def from_file(
        cls,
        filename: str,
        negative_occupancies: str = "warn",
        select_zero_weighted_k_points: bool = False,
    ) -> "Procar":
        """Create a :obj:`Procar` object from a VASP ``PROCAR`` file.

        Args:
            filename: Filename of the ``PROCAR`` file.
            negative_occupancies: How to handle negative occupancies.
                Options are:

                - ``'warn'`` (default): Warn that some partial occupancies
                  are negative.
                - ``'raise'``: Raise a ``ValueError``.
                - ``'zero'``: Set negative partial occupancies to zero.

            select_zero_weighted_k_points: Set to True to only read
                zero-weighted k-points from the file. Default is False.

        Returns:
            A :obj:`Procar` instance.
        """
        with open(filename, "r") as file_in:
            file_in.readline()
            number_of_k_points, number_of_bands, number_of_ions = [
                int(f) for f in get_numbers_from_string(file_in.readline())
            ]
            read_in = file_in.read()

        # Parse k-points
        k_points = k_point_parser(read_in)[:number_of_k_points]

        # Parse bands
        band_data = re.findall(
            r"band\s*(\d+)\s*#\s*energy\s*([-.\d]+)\s?\s*#\s*occ.\s*([-.\d]+)",
            read_in,
        )
        bands = np.array([
            Band(
                index=int(index),
                energy=float(energy),
                occupancy=float(occupancy),
                negative_occupancies=negative_occupancies,
            )
            for index, energy, occupancy in band_data
        ])

        # Parse projections and determine calculation type
        projection_data = projections_parser(read_in)
        calculation = {
            "non_spin_polarised": False,
            "non_collinear": False,
            "spin_polarised": False,
        }
        n_proj_rows = len(projection_data)
        expected = number_of_bands * number_of_k_points
        if n_proj_rows == expected:
            spin_channels = 1
            k_point_blocks = 1
            calculation["non_spin_polarised"] = True
        elif n_proj_rows == expected * 4:
            spin_channels = 4
            k_point_blocks = 1
            calculation["non_collinear"] = True
        elif n_proj_rows == expected * 2:
            spin_channels = 2
            k_point_blocks = 2
            calculation["spin_polarised"] = True
        else:
            raise ValueError(
                f"Cannot determine calculation type: {n_proj_rows} projection "
                f"rows for {number_of_k_points} k-points and "
                f"{number_of_bands} bands"
            )

        number_of_projections = (
            int(projection_data.shape[1] / (number_of_ions + 1)) - 1
        )

        # Reshape projection data into the 5D array
        if calculation["spin_polarised"]:
            data = (
                projection_data.reshape(
                    (
                        spin_channels,
                        number_of_k_points,
                        number_of_bands,
                        number_of_ions + 1,
                        number_of_projections + 1,
                    )
                )[:, :, :, :, 1:]
                .swapaxes(0, 1)
                .swapaxes(1, 2)
            )
        else:
            data = projection_data.reshape(
                (
                    number_of_k_points,
                    number_of_bands,
                    spin_channels,
                    number_of_ions + 1,
                    number_of_projections + 1,
                )
            )[:, :, :, :, 1:]

        pcar = cls(
            data=data,
            bands=bands,
            k_points=k_points,
            number_of_k_points=number_of_k_points,
            number_of_bands=number_of_bands,
            number_of_ions=number_of_ions,
            number_of_projections=number_of_projections,
            spin_channels=spin_channels,
            k_point_blocks=k_point_blocks,
            calculation=calculation,
            negative_occupancies=negative_occupancies,
        )
        pcar.sanity_check()
        if select_zero_weighted_k_points:
            k_point_indices = [
                i for i, kp in enumerate(pcar.k_points) if kp.weight == 0.0
            ]
            pcar = pcar.select_k_points(k_point_indices)
        return pcar

    @property
    def number_of_k_points(self) -> int:
        """The number of k-points described by this :obj:`Procar` object.

        Raises:
            ValueError: If the metadata count does not match the data array.
        """
        if self._number_of_k_points != self._data.shape[0]:
            raise ValueError(
                f"Number of k-points in metadata ({self._number_of_k_points}) "
                f"not equal to number in PROCAR data ({self._data.shape[0]})"
            )
        return self._number_of_k_points

    @property
    def number_of_bands(self) -> int:
        """The number of bands described by this :obj:`Procar` object.

        Raises:
            ValueError: If the metadata count does not match the data array.
        """
        if self._number_of_bands != self._data.shape[1]:
            raise ValueError(
                f"Number of bands in metadata ({self._number_of_bands}) "
                f"not equal to number in PROCAR data ({self._data.shape[1]})"
            )
        return self._number_of_bands

    @property
    def spin_channels(self) -> int:
        """The number of spin-channels described by this :obj:`Procar` object.

        Raises:
            ValueError: If the metadata count does not match the data array.
        """
        if self._spin_channels != self._data.shape[2]:
            raise ValueError(
                f"Number of spin channels in metadata ({self._spin_channels}) "
                f"not equal to number in PROCAR data ({self._data.shape[2]})"
            )
        return self._spin_channels

    @property
    def number_of_ions(self) -> int:
        """The number of ions described by this :obj:`Procar` object.

        Raises:
            ValueError: If the metadata count does not match the data array.
        """
        if self._number_of_ions != self._data.shape[3] - 1:
            raise ValueError(
                f"Number of ions in metadata ({self._number_of_ions}) "
                f"not equal to number in PROCAR data ({self._data.shape[3] - 1})"
            )
        return self._number_of_ions

    @property
    def number_of_projections(self) -> int:
        """The number of lm-projections described by this :obj:`Procar` object.

        Raises:
            ValueError: If the metadata count does not match the data array.
        """
        if self._number_of_projections != self._data.shape[4]:
            raise ValueError(
                f"Number of projections in metadata ({self._number_of_projections}) "
                f"not equal to number in PROCAR data ({self._data.shape[4]})"
            )
        return self._number_of_projections

    def print_weighted_band_structure(
        self,
        spins: list[int] | None = None,
        ions: list[int] | None = None,
        orbitals: list[int] | None = None,
        scaling: float = 1.0,
        e_fermi: float = 0.0,
        reciprocal_lattice: np.ndarray | None = None,
    ) -> None:
        """Print the weighted band structure to stdout.

        Args:
            spins: List of 1-based spin channel indices to include. Default
                is all spin channels.
            ions: List of ion indices to include. Default is the ``tot`` row
                (sum over all ions).
            orbitals: List of orbital projection indices to include. Default
                is all projections.
            scaling: Multiplicative scaling factor for the projection weights.
                Default is 1.0.
            e_fermi: Fermi energy in eV to subtract from all eigenvalues.
                Default is 0.0.
            reciprocal_lattice: 3x3 Cartesian reciprocal lattice used to
                compute real k-point spacings for the x-axis. If None,
                sequential integers are used.
        """
        band_structure_data = self.weighted_band_structure(
            spins=spins,
            ions=ions,
            orbitals=orbitals,
            scaling=scaling,
            e_fermi=e_fermi,
            reciprocal_lattice=reciprocal_lattice,
        )
        for i, band_data in enumerate(band_structure_data, 1):
            print(f"# band: {i}")
            for k_point_data in band_data:
                print(" ".join(str(f) for f in k_point_data))
            print()

    def weighted_band_structure(
        self,
        spins: list[int] | None = None,
        ions: list[int] | None = None,
        orbitals: list[int] | None = None,
        scaling: float = 1.0,
        e_fermi: float = 0.0,
        reciprocal_lattice: np.ndarray | None = None,
    ) -> np.ndarray:
        """Compute the weighted band structure array.

        Args:
            spins: List of 1-based spin channel indices to include. Default
                is all spin channels.
            ions: List of ion indices to include. Default is the ``tot`` row.
            orbitals: List of orbital projection indices to include. Default
                is all projections.
            scaling: Multiplicative scaling factor for the projection weights.
                Default is 1.0.
            e_fermi: Fermi energy in eV to subtract from all eigenvalues.
                Default is 0.0.
            reciprocal_lattice: 3x3 Cartesian reciprocal lattice for x-axis
                generation. If None, sequential integers are used.

        Returns:
            3D numpy array of shape ``(n_bands, n_k_points, 3)`` where the
            last axis is ``[x, energy, projection_weight]``.
        """
        if spins:
            spins = [s - 1 for s in spins]
        else:
            spins = list(range(self.spin_channels))
        if not ions:
            ions = [self.number_of_ions]  # nions+1 is the `tot` index
        if not orbitals:
            orbitals = list(range(self.number_of_projections))
        if self.calculation["spin_polarised"]:
            band_energies = (
                np.array([band.energy for band in self._bands])
                .reshape(
                    (self.spin_channels, self.number_of_k_points, self.number_of_bands)
                )[spins[0]]
                .T
            )
        else:
            band_energies = (
                np.array([band.energy for band in self._bands])
                .reshape((self.number_of_k_points, self.number_of_bands))
                .T
            )
        orbital_projection = np.sum(self._data[:, :, :, :, orbitals], axis=4)
        ion_projection = np.sum(orbital_projection[:, :, :, ions], axis=3)
        spin_projection = np.sum(ion_projection[:, :, spins], axis=2)
        x_axis = self.x_axis(reciprocal_lattice)
        rows = []
        for i in range(self.number_of_bands):
            for k, (e, p) in enumerate(
                zip(band_energies[i], spin_projection.T[i], strict=False)
            ):
                rows.append([x_axis[k], e - e_fermi, p * scaling])
        return np.array(rows).reshape((self.number_of_bands, -1, 3))

    def effective_mass_calc(
        self,
        k_point_indices: list[int],
        band_index: int,
        reciprocal_lattice: np.ndarray,
        spin: int = 1,
        printing: bool = False,
    ) -> float:
        """Calculate the effective mass at a band extremum.

        Args:
            k_point_indices: List of 1-based k-point indices to use.
            band_index: 1-based band index.
            reciprocal_lattice: 3x3 Cartesian reciprocal lattice in inverse Angstroms.
            spin: 1-based spin channel index. Default is 1.
            printing: If True, print k-point and eigenvalue data to stdout.
                Default is False.

        Returns:
            The effective mass in units of the free electron mass.

        Raises:
            ValueError: If ``spin`` exceeds the number of k-point blocks.
            ValueError: If fewer than 2 k-point indices are provided.
        """
        if spin > self._k_point_blocks:
            raise ValueError(
                f"spin index {spin} exceeds number of k-point blocks {self._k_point_blocks}"
            )
        if len(k_point_indices) < 2:
            raise ValueError("at least 2 k-point indices are required for effective mass calculation")
        band_energies = self._bands[:, 1:].reshape(
            self._k_point_blocks, self.number_of_k_points, self.number_of_bands
        )
        frac_k_point_coords = np.array(
            [self._k_points[k - 1].frac_coords for k in k_point_indices]
        )
        eigenvalues = np.array(
            [band_energies[spin - 1][k - 1][band_index - 1] for k in k_point_indices]
        )
        if printing:
            print("# h k l e")
            for row in np.concatenate(
                (frac_k_point_coords, np.array([eigenvalues]).T), axis=1
            ):
                print(" ".join(str(f) for f in row))
        reciprocal_lattice = reciprocal_lattice * 2 * math.pi * angstrom_to_bohr
        cart_k_point_coords = np.array(
            [k.cart_coords(reciprocal_lattice) for k in self._k_points]
        )  # convert k-points to cartesian
        if len(k_point_indices) == 2:
            effective_mass_function = two_point_effective_mass
        else:
            effective_mass_function = least_squares_effective_mass
        return effective_mass_function(cart_k_point_coords, eigenvalues)

    def x_axis(self, reciprocal_lattice: np.ndarray | None = None) -> np.ndarray:
        """Generate the x-axis values for a band-structure plot.

        Returns an array of cumulative distances in reciprocal space between
        sequential k-points.

        Args:
            reciprocal_lattice: 3x3 Cartesian reciprocal lattice. If None,
                the returned x-axis values will be sequential integers, giving
                even spacings between sequential k-points.

        Returns:
            Array of x-axis values.
        """
        if reciprocal_lattice is not None:
            cartesian_k_points = np.array(
                [k.cart_coords(reciprocal_lattice) for k in self._k_points]
            )
            x_axis = [0.0]
            for i in range(1, len(cartesian_k_points)):
                dk = cartesian_k_points[i - 1] - cartesian_k_points[i]
                mod_dk = np.sqrt(np.dot(dk, dk))
                x_axis.append(mod_dk + x_axis[-1])
            x_axis_array = np.array(x_axis)
        else:
            x_axis_array = np.arange(len(self._k_points))
        return x_axis_array

    @property
    def bands(self) -> np.ndarray:
        """3D array of Band objects shaped (k_point_blocks, n_k_points, n_bands)."""
        return self._bands.reshape(
            self._k_point_blocks, self._number_of_k_points, self.number_of_bands
        )

    @property
    def k_points(self) -> list[KPoint]:
        """List of KPoint objects for this Procar."""
        return self._k_points

    def select_bands_by_kpoint(self, band_indices: list[int]) -> np.ndarray:
        """Return a flattened array of Band objects at the specified k-points.

        Args:
            band_indices: List of k-point indices (0-based) to select.

        Returns:
            Flattened numpy array of Band objects.
        """
        return np.ravel(self.bands[:, band_indices, :])

    def select_k_points(self, band_indices: list[int]) -> "Procar":
        """Return a new Procar containing only the specified k-points.

        Args:
            band_indices: List of k-point indices (0-based) to keep.

        Returns:
            A new :obj:`Procar` instance with the selected k-points.
        """
        new_procar = deepcopy(self)
        new_procar._bands = np.ravel(new_procar.bands[:, band_indices, :])
        new_procar._data = np.array(
            [kp for i, kp in enumerate(new_procar._data) if i in band_indices]
        )
        new_procar._number_of_k_points = len(band_indices)
        new_procar._k_points = [
            kp for i, kp in enumerate(new_procar._k_points) if i in band_indices
        ]
        for i, kp in enumerate(new_procar._k_points, 1):
            kp.index = i
        new_procar.sanity_check()
        return new_procar

