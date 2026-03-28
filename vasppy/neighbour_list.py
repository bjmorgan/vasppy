"""Neighbour-list utilities following the scheme in Rabani et al.

Reference:
    Rabani et al. J. Chem. Phys. 1997
    doi: https://doi.org/10.1063/1.474927
"""

from __future__ import annotations

from pymatgen.core import Structure
from vasppy.utils import dr_ij
from typing import Any, cast
import numpy as np
from numpy.typing import NDArray


class NeighbourList:
    """A set of neighbour lists for a collection of central atoms.

    Attributes:
        vectors: An M x N boolean-integer array of neighbour lists. Each row
            corresponds to one central atom (from ``indices_i``) and each
            column to one potential neighbour (from ``indices_j``). An entry
            of 1 indicates a neighbour within the cutoff radius. Constructed
            following Rabani et al. J. Chem. Phys. 1997
            doi: https://doi.org/10.1063/1.474927
    """

    vectors: NDArray[np.integer[Any]]

    def __init__(
        self,
        structure: Structure,
        indices_i: list[int],
        indices_j: list[int],
        r_cut: float,
    ) -> None:
        """Initialise a NeighbourList instance.

        Args:
            structure: Pymatgen Structure object to parse.
            indices_i: List of indices of central atoms.
            indices_j: List of indices of potential neighbour atoms.
            r_cut: Neighbour cutoff distance.
        """
        all_dr_ij = dr_ij(
            structure=structure,
            indices_i=indices_i,
            indices_j=indices_j,
            self_reference=False,
        )
        self.vectors = (all_dr_ij <= r_cut).astype(int)

    @property
    def coordination_numbers(self) -> NDArray[np.integer[Any]]:
        """Return the coordination number of each central atom.

        Returns:
            A 1D array of coordination numbers, one per site in ``indices_i``.
        """
        return cast(NDArray[np.integer[Any]], np.sum(self.vectors, axis=1))

    def __eq__(self, other: object) -> bool:
        """Test whether two NeighbourList objects have equal vectors.

        Args:
            other: Object to compare against.

        Returns:
            ``True`` if the neighbour-list vectors are element-wise equal.
        """
        if not isinstance(other, NeighbourList):
            return NotImplemented
        return bool((self.vectors == other.vectors).all())

    @classmethod
    def from_species_strings(
        cls: type[NeighbourList],
        structure: Structure,
        species_i: str,
        species_j: str,
        r_cut: float,
    ) -> NeighbourList:
        """Initialise a NeighbourList by specifying species strings.

        Args:
            structure: A pymatgen Structure.
            species_i: String for species i, e.g. ``"Na"``.
            species_j: String for species j, e.g. ``"Cl"``.
            r_cut: Neighbour cutoff radius.

        Returns:
            A NeighbourList instance.
        """
        indices_i = [
            i for i, site in enumerate(structure) if site.species_string == species_i
        ]
        indices_j = [
            j for j, site in enumerate(structure) if site.species_string == species_j
        ]
        return cls(
            structure=structure, indices_i=indices_i, indices_j=indices_j, r_cut=r_cut
        )


def neighbour_list_correlation(
    nlist_i: NeighbourList, nlist_j: NeighbourList
) -> NDArray[np.floating[Any]]:
    """Compute the normalised correlation between two NeighbourList objects.

    For each neighbour-list vector, computes (l_i · l_j) / (l_i · l_i).
    See Rabani et al. J. Chem. Phys. 1997 doi:https://doi.org/10.1063/1.474927
    Eqn. 7 for details.

    Args:
        nlist_i: A NeighbourList object.
        nlist_j: A NeighbourList object.

    Returns:
        A 1D array of normalised correlation terms.

    Raises:
        ValueError: If the two NeighbourList objects have different vector
            shapes.
    """
    if nlist_i.vectors.shape != nlist_j.vectors.shape:
        raise ValueError(
            f"NeighbourList vector shapes are not equal: "
            f"{nlist_i.vectors.shape} != {nlist_j.vectors.shape}"
        )
    return cast(
        NDArray[np.floating[Any]],
        np.einsum("ij,ij->i", nlist_i.vectors, nlist_j.vectors)
        / np.einsum("ij,ij->i", nlist_i.vectors, nlist_i.vectors),
    )


def neighbour_list_n_out(
    nlist_i: NeighbourList, nlist_j: NeighbourList
) -> NDArray[np.floating[Any]]:
    """Compute n :sup:`out` between two NeighbourList objects.

    For each neighbour-list vector, computes (l_i · l_i) - (l_i · l_j).
    See Rabani et al. J. Chem. Phys. 1997 doi:https://doi.org/10.1063/1.474927
    Eqn. 8 for details.

    Args:
        nlist_i: A NeighbourList object for neighbour lists at time 0.
        nlist_j: A NeighbourList object for neighbour lists at time t.

    Returns:
        A 1D array of n :sup:`out` values.

    Raises:
        ValueError: If the two NeighbourList objects have different vector
            shapes.
    """
    if nlist_i.vectors.shape != nlist_j.vectors.shape:
        raise ValueError(
            f"NeighbourList vector shapes are not equal: "
            f"{nlist_i.vectors.shape} != {nlist_j.vectors.shape}"
        )
    return cast(
        NDArray[np.floating[Any]],
        np.einsum("ij,ij->i", nlist_i.vectors, nlist_i.vectors)
        - np.einsum("ij,ij->i", nlist_i.vectors, nlist_j.vectors),
    )


def neighbour_list_n_in(
    nlist_i: NeighbourList, nlist_j: NeighbourList
) -> NDArray[np.floating[Any]]:
    """Compute n :sup:`in` between two NeighbourList objects.

    For each neighbour-list vector, computes (l_j · l_j) - (l_i · l_j).
    See Rabani et al. J. Chem. Phys. 1997 doi:https://doi.org/10.1063/1.474927
    Eqn. 9 for details.

    Args:
        nlist_i: A NeighbourList object for neighbour lists at time 0.
        nlist_j: A NeighbourList object for neighbour lists at time t.

    Returns:
        A 1D array of n :sup:`in` values.

    Raises:
        ValueError: If the two NeighbourList objects have different vector
            shapes.
    """
    if nlist_i.vectors.shape != nlist_j.vectors.shape:
        raise ValueError(
            f"NeighbourList vector shapes are not equal: "
            f"{nlist_i.vectors.shape} != {nlist_j.vectors.shape}"
        )
    return cast(
        NDArray[np.floating[Any]],
        np.einsum("ij,ij->i", nlist_j.vectors, nlist_j.vectors)
        - np.einsum("ij,ij->i", nlist_i.vectors, nlist_j.vectors),
    )
