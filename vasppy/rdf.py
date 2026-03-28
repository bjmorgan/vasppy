from __future__ import annotations

from typing import Any
import numpy as np
from scipy.ndimage import gaussian_filter1d
from pymatgen.core import Structure
from vasppy.utils import dr_ij
from numpy.typing import NDArray


class RadialDistributionFunction:
    """Class for computing radial distribution functions.

    Attributes:
        nbins: Number of bins.
        range: Minimum and maximum values of r.
        intervals: r values of the bin edges.
        dr: Bin width.
        r: Mid-points of each bin.
        rdf: RDF values.
        coordination_number: Volume integral of the RDF.
    """

    def __init__(
        self,
        structures: list[Structure],
        indices_i: list[int],
        indices_j: list[int] | None = None,
        nbins: int = 500,
        r_min: float = 0.0,
        r_max: float = 10.0,
        weights: list[float] | None = None,
    ) -> None:
        """Initialise a RadialDistributionFunction instance.

        Args:
            structures: List of pymatgen Structure objects.
            indices_i: List of indices for species i.
            indices_j: List of indices for species j. Defaults to ``None``
                (uses the same indices as species i).
            nbins: Number of bins used for the RDF. Defaults to ``500``.
            r_min: Minimum r value. Defaults to ``0.0``.
            r_max: Maximum r value. Defaults to ``10.0``.
            weights: List of weights for each structure. Defaults to ``None``
                (equal weights of 1.0 are used).

        Raises:
            ValueError: If the length of ``weights`` does not match the number
                of structures.
        """
        if isinstance(indices_i, np.ndarray):
            indices_i = indices_i.tolist()
        if indices_j is not None and isinstance(indices_j, np.ndarray):
            indices_j = indices_j.tolist()
        if weights is not None:
            if len(weights) != len(structures):
                raise ValueError(
                    "List of structure weights needs to be the same length"
                    " as the list of structures."
                )
        else:
            weights = [1.0] * len(structures)
        if indices_j is None:
            indices_j = indices_i
        self.indices_i = indices_i
        self.indices_j = indices_j
        self.nbins = nbins
        self.range = (r_min, r_max)
        self.intervals = np.linspace(r_min, r_max, nbins + 1)
        self.dr = (r_max - r_min) / nbins
        self.r = self.intervals[:-1] + self.dr / 2.0
        ff = shell_volumes(self.intervals)
        self.coordination_number = np.zeros(nbins)
        self.rdf = np.zeros((nbins,), dtype=np.double)
        for structure, weight in zip(structures, weights):
            all_dr_ij = dr_ij(
                structure=structure,
                indices_i=self.indices_i,
                indices_j=self.indices_j,
                self_reference=False,
            ).flatten()
            hist = np.histogram(
                all_dr_ij, bins=nbins, range=(r_min, r_max), density=False
            )[0]
            rho = float(len(self.indices_i)) / structure.lattice.volume
            self.rdf += hist * weight / rho
            self.coordination_number += np.cumsum(hist)
        self.rdf /= ff * sum(weights) * float(len(indices_j))
        self.coordination_number /= sum(weights) * float(len(self.indices_j))

    def smeared_rdf(self, sigma: float = 0.1) -> NDArray[np.floating[Any]]:
        """Smear the RDF with a Gaussian kernel.

        Args:
            sigma: Standard deviation for the Gaussian kernel. Defaults to
                ``0.1``.

        Returns:
            Smeared RDF data.
        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.rdf, sigma=sigma_n_bins)

    @classmethod
    def from_species_strings(
        cls: type[RadialDistributionFunction],
        structures: list[Structure],
        species_i: str,
        species_j: str | None = None,
        **kwargs,
    ) -> RadialDistributionFunction:
        """Initialise a RadialDistributionFunction by specifying species strings.

        Args:
            structures: List of pymatgen Structure objects.
            species_i: String for species i, e.g. ``"Na"``.
            species_j: String for species j, e.g. ``"Cl"``. Defaults to
                ``None``.
            **kwargs: Additional keyword arguments passed to
                :class:`RadialDistributionFunction`.

        Returns:
            A RadialDistributionFunction instance.

        Raises:
            ValueError: If species i is not found in the first structure.
            ValueError: If species j is not found in the first structure.
        """
        indices_i: list[int]
        indices_j: list[int] | None

        indices_i = [
            i
            for i, site in enumerate(structures[0])
            if site.species_string == species_i
        ]
        if species_j:
            indices_j = [
                j
                for j, site in enumerate(structures[0])
                if site.species_string == species_j
            ]
        else:
            indices_j = None

        if not indices_i:
            raise ValueError("Species i not found.")
        if indices_j is not None and not indices_j:
            raise ValueError("Species j not found.")

        return cls(
            structures=structures, indices_i=indices_i, indices_j=indices_j, **kwargs
        )


class VanHoveAnalysis:
    """Class for computing Van Hove correlation functions.

    Attributes:
        nbins: Number of bins.
        range: Minimum and maximum values of r.
        intervals: r values of the bin edges.
        dr: Bin width.
        r: Mid-points of each bin.
        gsrt: Self part of the Van Hove correlation function.
        gdrt: Distinct part of the Van Hove correlation function.
    """

    def __init__(
        self,
        structures: list[Structure],
        indices: list[int],
        d_steps: int,
        nbins: int = 500,
        r_min: float = 0.0,
        r_max: float = 10.0,
    ) -> None:
        """Initialise a VanHoveAnalysis instance.

        Args:
            structures: List of pymatgen Structure objects.
            indices: List of site indices for the species to consider.
            d_steps: Number of steps between structures at dt=0 and dt=t.
            nbins: Number of bins used for the correlation functions. Defaults
                to ``500``.
            r_min: Minimum r value. Defaults to ``0.0``.
            r_max: Maximum r value. Defaults to ``10.0``.
        """
        self.nbins = nbins
        self.range = (r_min, r_max)
        self.intervals = np.linspace(r_min, r_max, nbins + 1)
        self.dr = (r_max - r_min) / nbins
        self.r = self.intervals[:-1] + self.dr / 2.0
        self.gdrt = np.zeros((nbins), dtype=np.double)
        self.gsrt = np.zeros((nbins), dtype=np.double)
        lattice = structures[0].lattice
        ff = shell_volumes(self.intervals)
        rho = len(indices) / lattice.volume
        for struc_i, struc_j in zip(
            structures[: len(structures) - d_steps], structures[d_steps:]
        ):
            i_frac_coords = struc_i.frac_coords[indices]
            j_frac_coords = struc_j.frac_coords[indices]
            all_dr_ij = lattice.get_all_distances(i_frac_coords, j_frac_coords)
            mask = np.ones(all_dr_ij.shape, dtype=bool)
            np.fill_diagonal(mask, 0)
            distinct_dr_ij = np.ndarray.flatten(all_dr_ij[mask])
            hist = np.histogram(
                distinct_dr_ij, bins=nbins, range=(0.0, r_max), density=False
            )[0]
            self.gdrt += hist / rho
            self_dr_ij = np.ndarray.flatten(all_dr_ij[np.invert(mask)])
            hist = np.histogram(
                self_dr_ij, bins=nbins, range=(0.0, r_max), density=False
            )[0]
            self.gsrt += hist / rho
        self.gdrt /= ff * (len(structures) - d_steps) * float(len(indices))
        self.gsrt /= (len(structures) - d_steps) * float(len(indices))

    def self(self, sigma: float | None = None) -> np.ndarray:
        """Return the self part of the Van Hove correlation function.

        Args:
            sigma: Optional smearing width.

        Returns:
            The (optionally smeared) self correlation function.
        """
        if sigma:
            return self.smeared_gsrt(sigma=sigma)
        return self.gsrt

    def distinct(self, sigma: float | None = None) -> np.ndarray:
        """Return the distinct part of the Van Hove correlation function.

        Args:
            sigma: Optional smearing width.

        Returns:
            The (optionally smeared) distinct correlation function.
        """
        if sigma:
            return self.smeared_gdrt(sigma=sigma)
        return self.gdrt

    def smeared_gsrt(self, sigma: float = 0.1) -> NDArray[np.floating[Any]]:
        """Smear the self part of the Van Hove correlation function.

        Args:
            sigma: Standard deviation for the Gaussian kernel. Defaults to
                ``0.1``.

        Returns:
            Smeared self correlation function data.
        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.gsrt, sigma=sigma_n_bins)

    def smeared_gdrt(self, sigma: float = 0.1) -> NDArray[np.floating[Any]]:
        """Smear the distinct part of the Van Hove correlation function.

        Args:
            sigma: Standard deviation for the Gaussian kernel. Defaults to
                ``0.1``.

        Returns:
            Smeared distinct correlation function data.
        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.gdrt, sigma=sigma_n_bins)


def shell_volumes(intervals: NDArray[np.floating[Any]]) -> NDArray[np.floating[Any]]:
    """Compute the volumes of concentric spherical shells.

    Args:
        intervals: N radial boundaries used to define the set of N-1 shells.

    Returns:
        Volumes of each shell.
    """
    return 4.0 / 3.0 * np.pi * (intervals[1:] ** 3 - intervals[:-1] ** 3)
