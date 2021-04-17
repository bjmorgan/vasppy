from __future__ import annotations

import numpy as np  # type: ignore
from scipy.ndimage.filters import gaussian_filter1d  # type: ignore
from pymatgen.core import Structure  
from typing import List, Optional, TypeVar, Type
from numpy.typing import ArrayLike


"""
This module provides classes for calculating radial distribution functions
and Van Hove correlation functions.
"""

class NeighbourList(object):
    """
    Class for computing a set of neighbour lists.
    
    Attributes:
        vectors (np.ndarray): An array of M neighbour lists, with each list a 1D
            vector over N potential neighbours, stored as a M x N numpy array.
            Each neighbour list is constructed following the scheme in 
            Rabani et al. J. Chem. Phys. 1997
            doi: https://doi.org/10.1063/1.474927
        
    """

    def __init__(self,
                 structure: Structure,
                 indices_i: List[int],
                 indices_j: List[int],
                 r_cut: float) -> None:
        """
        Initialise a NeighbourList instance.
    
        Args:
            structure (pymatgen.Structure): Pymatgen Structure object to parse.
            indices_i (list(int)): List of indices of central atoms.
            indices_j (list(int)): List of indices of potential neighbour atoms.
            r_cut (float): Neighbour cutoff distance. 

        """
        all_dr_ij = dr_ij(structure=structure,
                          indices_i=indices_i,
                          indices_j=indices_j,
                          self_reference=False)
        self.vectors = (all_dr_ij <= r_cut).astype(int)
        
    @property
    def coordination_numbers(self) -> np.ndarray:
        """
        Return the coordination number of each site i.
        
        Args:
            None
            
        Returns:
            None
            
        """
        return np.sum(self.vectors, axis=1)
        
    @classmethod
    def from_species_strings(cls: Type[NeighbourList],
                             structure: Structure,
                             species_i: str,
                             species_j: str,
                             r_cut: float) -> NeighbourList:
        """
        Initialise a NeighbourList instance by specifying species strings.
        
        Args:
            structure (pymatgen.Structure): A pymatgen Structure.
            species_i (str): String for species i, e.g., ``"Na"``.
            species_j (str): String for species j, e.g., ``"Cl"``.
            r_cut (float): Neighbour cutoff radius.
         
        Returns:
            (NeighbourList)
            
        """                   
        indices_i = [i for i, site in 
                     enumerate(structure)
                     if site.species_string is species_i]
        indices_j = [j for j, site in 
                     enumerate(structure)
                     if site.species_string is species_j]
        return cls(structure=structure,
                   indices_i=indices_i,
                   indices_j=indices_j,
                   r_cut=r_cut)

def neighbour_list_correlation(nlist_i: NeighbourList,
                               nlist_j: NeighbourList) -> np.ndarray:
    """
    Compute the normalised correlation between two NeighbourList object.
    
    Args:
        nlist_i (NeighbourList): A NeighbourList object.
        nlist_j (NeighbourList): A NeighbourList object.
        
    Returns:
        (np.ndarray(float)): A 1D array of normalised correlation terms.
        
    Raises:
        ValueError: If the two NeighbourList objects have different numbers
            of lengths of neighbour list vectors.
        
    Note:
        For each neighbour list vector, computes (l_i.l_j)/(l_i.l_i).
        See Rabani et al. J. Chem. Phys. 1997 doi:https://doi.org/10.1063/1.474927
        for details.

    """
    if nlist_i.vectors.shape != nlist_j.vectors.shape:
        raise ValueError(f'NeighbourList vector shapes are not equal: {nlist_i.vectors.shape} != {nlist_j.vectors.shape}')
    return (np.einsum('ij,ij->i', nlist_i.vectors, nlist_j.vectors) / 
            np.einsum('ij,ij->i', nlist_i.vectors, nlist_i.vectors))


class RadialDistributionFunction(object):
    """
    Class for computing radial distribution functions.

    Attributes:
        nbins (int): Number of bins.
        range ((float, float)): Minimum and maximum values of r.
        intervals (np.array(float)): r values of the bin edges.
        dr (float): bin width.
        r (float): mid-points of each bin.
        rdf (np.array(float)): RDF values.
        coordination_number (np.array(float)): Volume integral of the RDF.

    """

    def __init__(self,
                 structures: List[Structure],
                 indices_i: List[int],
                 indices_j: Optional[List[int]] = None,
                 nbins: int = 500,
                 r_min: float = 0.0,
                 r_max: float = 10.0,
                 weights: Optional[List[float]] = None) -> None:
        """
        Initialise a RadialDistributionFunction instance.

        Args:
            structures (list(pymatgen.Structure)): List of pymatgen Structure objects.
            indices_i (list(int)): List of indices for species i.
            indices_j (:obj:`list(int)`, optional): List of indices for species j. Optional,
                default is `None`.
            nbins (:obj:`int`, optional): Number of bins used for the RDF. Optional, default is 500.
            rmin (:obj:`float`, optional): Minimum r value. Optional, default is 0.0.
            rmax (:obj:`float`, optional): Maximum r value. Optional, default is 10.0.
            weights (:obj:`list(float)`, optional): List of weights for each structure.
                Optional, default is `None`.

        Returns:
             None

        """
        if weights:
            if len(weights) != len(structures):
                raise ValueError('List of structure weights needs to be the same length'
                                 ' as the list of structures.')
        else:
            weights = [1.0] * len(structures)
        self.self_reference = (not indices_j) or (indices_j == indices_i)
        if not indices_j:
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
        self.rdf = np.zeros((nbins), dtype=np.double)
        for structure, weight in zip(structures, weights):
            all_dr_ij = dr_ij(structure=structure,
                              indices_i=self.indices_i,
                              indices_j=self.indices_j,
                              self_reference=False).flatten()
            hist = np.histogram(all_dr_ij,
                                bins=nbins,
                                range=(r_min, r_max),
                                density=False)[0]
            rho = float(len(self.indices_i)) / structure.lattice.volume
            self.rdf += hist * weight / rho
            self.coordination_number += np.cumsum(hist)
        self.rdf = self.rdf / ff / sum(weights) / float(len(indices_j))
        self.coordination_number = self.coordination_number / \
            sum(weights) / float(len(self.indices_j))

    def smeared_rdf(self, 
                    sigma: float = 0.1) -> np.ndarray:
        """
        Smear the RDF with a Gaussian kernel.

        Args:
            sigma (:obj:`float`, optional): Standard deviation for Gaussian kernel. 
                Optional, default is 0.1.

        Returns:
            (np.array): Smeared RDF data.

        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.rdf, sigma=sigma_n_bins)

    @classmethod
    def from_species_strings(cls: Type[RadialDistributionFunction], 
                             structures: List[Structure], 
                             species_i: str, 
                             species_j: Optional[str] = None, 
                             **kwargs) -> RadialDistributionFunction:
        """
        Initialise a RadialDistributionFunction instance by specifying species strings.

        Args:
            structures (list(pymatgen.Structure)): List of pymatgen Structure objects.
            species_i (str): String for species i, e.g. ``"Na"``.
            species_j (:obj:`str`, optional): String for species j, e.g. ``"Cl"``. Optional,
                default is `None`. 
            **kwargs: Variable length keyword argument list. 
                See :func:`vasppy.rdf.RadialDistributionFunction`
                for the full list of accepted arguments.

        Returns:
            (RadialDistributionFunction)

        """
        indices_i: List[int]
        indices_j: Optional[List[int]]

        indices_i = [i for i, site in 
                     enumerate(structures[0]) if site.species_string is species_i]
        if species_j:
            indices_j = [j for j, site in 
                         enumerate(structures[0]) if site.species_string is species_j]
        else:
            indices_j = None
        return cls(structures=structures, 
                   indices_i=indices_i, 
                   indices_j=indices_j, 
                   **kwargs)


class VanHoveAnalysis(object):
    """
    Class for computing Van Hove correlation functions.

    Attributes:
        nbins (int): Number of bins.
        range ((float, float)): Minimum and maximum values of r.
        intervals (np.array(float)): r values of the bin edges.
        dr (float): bin width.
        r (float): mid-points of each bin.
        gsrt (np.array(float)): Self part of the Van Hove correlation function.
        gdrt (np.array(float)): Distinct part of the Van Hove correlation function.

    """

    def __init__(self,
                 structures: List[Structure],
                 indices: List[int],
                 d_steps: int,
                 nbins: int = 500,
                 r_min: float = 0.0,
                 r_max: float = 10.0):
        """
        Initialise a VanHoveCorrelationFunction instance.

        Args:
            structures (list(pymatgen.Structure)): List of pymatgen Structure objects.
            indices (list(int)): List of indices for species to consider.
            d_steps (int): number of steps between structures at dt=0 and dt=t.
            nbins (:obj:`int`, optional): Number of bins used for the RDF. Optional, default is 500.
            rmin (:obj:`float`, optional): Minimum r value. Optional, default is 0.0.
            rmax (:obj:`float`, optional): Maximum r value. Optional, default is 10.0.

        Returns:
             None

        """
        self.nbins = nbins
        self.range = (r_min, r_max)
        self.intervals = np.linspace(r_min, r_max, nbins + 1)
        self.dr = (r_max - r_min) / nbins
        self.r = self.intervals[:-1] + self.dr / 2.0
        self.gdrt = np.zeros((nbins), dtype=np.double)
        self.gsrt = np.zeros((nbins), dtype=np.double)
        rho = len(indices) / structures[0].lattice.volume
        lattice = structures[0].lattice
        ff = shell_volumes(self.intervals)
        rho = len(indices) / lattice.volume
        for struc_i, struc_j in zip(structures[:len(structures) - d_steps], structures[d_steps:]):
            i_frac_coords = struc_i.frac_coords[indices]
            j_frac_coords = struc_j.frac_coords[indices]
            dr_ij = lattice.get_all_distances(i_frac_coords, j_frac_coords)
            mask = np.ones(dr_ij.shape, dtype=bool)
            np.fill_diagonal(mask, 0)
            distinct_dr_ij = np.ndarray.flatten(dr_ij[mask])
            hist = np.histogram(distinct_dr_ij, bins=nbins,
                                range=(0.0, r_max), density=False)[0]
            self.gdrt += hist / rho
            self_dr_ij = np.ndarray.flatten(dr_ij[np.invert(mask)])
            hist = np.histogram(self_dr_ij, bins=nbins,
                                range=(0.0, r_max), density=False)[0]
            self.gsrt += hist / rho
        self.gdrt = self.gdrt / ff / \
            (len(structures) - d_steps) / float(len(indices))
        self.gsrt = self.gsrt / \
            (len(structures) - d_steps) / float(len(indices))

    def self(self, 
             sigma: Optional[float] = None) -> np.ndarray:
        if sigma:
            return self.smeared_gsrt(sigma=sigma)
        else:
            return self.gsrt

    def distinct(self,
                 sigma: Optional[float] = None) -> np.ndarray:
        if sigma:
            return self.smeared_gdrt(sigma=sigma)
        else:
            return self.gdrt

    def smeared_gsrt(self, 
                     sigma: float = 0.1) -> np.ndarray:
        """
        Smear the self part of the Van Hove correlation function with a Gaussian kernel.

        Args:
            sigma (:obj:`float`, optional): Standard deviation for Gaussian kernel. Optional, default is 0.1.

        Returns:
            (np.array): Smeared data.

        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.gsrt, sigma=sigma_n_bins)

    def smeared_gdrt(self,
                     sigma: float = 0.1) -> np.ndarray:
        """
        Smear the distinct part of the Van Hove correlation function with a Gaussian kernel.

        Args:
            sigma (:obj:`float`, optional): Standard deviation for Gaussian kernel. Optional, default is 0.1.

        Returns:
            (np.array): Smeared data.

        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.gdrt, sigma=sigma_n_bins)


def shell_volumes(intervals: np.ndarray) -> np.ndarray:
    """Volumes of concentric spherical shells.

    Args:
        intervals (np.array): N radial boundaries used to define the set of N-1 shells.

    Returns:
        np.array: Volumes of each shell.

    """
    return 4.0 / 3.0 * np.pi * (intervals[1:]**3 - intervals[:-1]**3)
    
def dr_ij(structure: Structure,
      indices_i: Optional[List[int]] = None,
      indices_j: Optional[List[int]] = None,
      self_reference: bool = False) -> np.ndarray:
    """
    Calculate all i-j interatomic distances for a single pymatgen Structure.

    Args:
        structure (:obj:`pymatgen.Structure`): A pymatgen Structure
        indices_i (:obj:`list(int)`, optional): List of indices for species i.
            Optional, default is `None`.
            If `indices_i` is not specified, then distances will be calculated between 
            all pairs of atoms.
        indices_j (:obj:`list(int)`, optional): List of indices for species j.
            Optional, default is `None`.
            If `indices_j` is not specified, then `indices_j` will be set equal
            to `indices_i`.
        self_reference (bool, optional): If computing distances for i==j, whether
            to include the i==j dr=0 terms. Optional, default is `False`.
            
    Returns:
        np.array: N_i x N_j numpy array of i-j minimum image distances.

    """
    if indices_i is None:
        indices_i = list(range(len(structure)))
    if indices_j is None:
        indices_j = indices_i
    lattice = structure.lattice
    i_frac_coords = structure.frac_coords[indices_i]
    j_frac_coords = structure.frac_coords[indices_j]
    dr_ij = lattice.get_all_distances(i_frac_coords, j_frac_coords)
    # If indices_i and indices_j contain common elements AND self_reference == False
    # then we want to mask dr_ij to remove the i==j dr=0 terms
    if (np.intersect1d(indices_i, indices_j).size > 0) and not self_reference:
        mask = np.ones_like(dr_ij, dtype=bool)
        for i_loc, i in enumerate(indices_i):
            for j_loc, j in enumerate(indices_j):
                if i == j:
                    mask[i_loc, j_loc] = 0
        to_return = dr_ij[mask].reshape(len(indices_i), -1)
    else:
        to_return = dr_ij
    return to_return

