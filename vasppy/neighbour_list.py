from __future__ import annotations

from pymatgen.core import Structure
from vasppy.utils import dr_ij
from typing import List, Type
import numpy as np

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
	
    def __eq__(self,
               other: NeighbourList) -> bool:
        """Test whether two NeighbourList objects have equal vectors."""
        return (self.vectors == other.vectors).all()

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
        Eqn. 7 for details.

    """
    if nlist_i.vectors.shape != nlist_j.vectors.shape:
        raise ValueError(f'NeighbourList vector shapes are not equal: {nlist_i.vectors.shape} != {nlist_j.vectors.shape}')
    return (np.einsum('ij,ij->i', nlist_i.vectors, nlist_j.vectors) / 
            np.einsum('ij,ij->i', nlist_i.vectors, nlist_i.vectors))

def neighbour_list_n_out(nlist_i: NeighbourList,
                         nlist_j: NeighbourList) -> np.ndarray:
    """
    Compute n^out between two NeighbourList object.
    
    Args:
        nlist_i (NeighbourList): A NeighbourList object for neighbour lists at time 0.
        nlist_j (NeighbourList): A NeighbourList object for neighbour lists at time t.
    
    Returns:
        (np.ndarray(float)): A 1D array of normalised correlation terms.
    
    Raises:
        ValueError: If the two NeighbourList objects have different numbers
            of lengths of neighbour list vectors.    
    
    Note:
        For each neighbour list vector, computes (l_i.l_i) - (l_i.l_j).
        See Rabani et al. J. Chem. Phys. 1997 doi:https://doi.org/10.1063/1.474927
        Eqn. 8 for details.
    
    """
    if nlist_i.vectors.shape != nlist_j.vectors.shape:
        raise ValueError(f'NeighbourList vector shapes are not equal: {nlist_i.vectors.shape} != {nlist_j.vectors.shape}')
    return (np.einsum('ij,ij->i', nlist_i.vectors, nlist_i.vectors) - 
            np.einsum('ij,ij->i', nlist_i.vectors, nlist_j.vectors))
            
def neighbour_list_n_in(nlist_i: NeighbourList,
                        nlist_j: NeighbourList) -> np.ndarray:
    """
    Compute n^in between two NeighbourList object.
    
    Args:
        nlist_i (NeighbourList): A NeighbourList object for neighbour lists at time 0.
        nlist_j (NeighbourList): A NeighbourList object for neighbour lists at time t.
    
    Returns:
        (np.ndarray(float)): A 1D array of normalised correlation terms.
    
    Raises:
        ValueError: If the two NeighbourList objects have different numbers
            of lengths of neighbour list vectors.    
    
    Note:
        For each neighbour list vector, computes (l_i.l_i) - (l_i.l_j).
        See Rabani et al. J. Chem. Phys. 1997 doi:https://doi.org/10.1063/1.474927
        Eqn. 9 for details.
    
    """
    if nlist_i.vectors.shape != nlist_j.vectors.shape:
        raise ValueError(f'NeighbourList vector shapes are not equal: {nlist_i.vectors.shape} != {nlist_j.vectors.shape}')
    return (np.einsum('ij,ij->i', nlist_j.vectors, nlist_j.vectors) - 
            np.einsum('ij,ij->i', nlist_i.vectors, nlist_j.vectors)) 
