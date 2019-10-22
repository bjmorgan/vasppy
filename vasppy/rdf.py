import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

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
    
    def __init__(self, structures, indices_i, indices_j=None, nbins=500, r_min=0.0, r_max=10.0):
        """
        Initialise a RadialDistributionFunction instance.

        Args:
            structures (list(pymatgen.Structure)): List of pymatgen Structure objects.
            indices_i (list(int)): List of indices for species i.
            indices_j (:obj:list(int), optional): List of indices for species j. Optional,
                default is `None`.
            nbins (:obj:`int`, optional): Number of bins used for the RDF. Optional, default is 500.
            rmin (:obj:`float`, optional): Minimum r value. Optional, default is 0.0.
            rmax (:obj:`float`, optional): Maximum r value. Optional, default is 10.0.

        Returns:
             None

        """
        self_reference = (indices_j is None) or (indices_j == indices_i)
        self.nbins = nbins
        self.range = (r_min, r_max)
        self.intervals = np.linspace(r_min, r_max, nbins+1)
        self.dr = (r_max - r_min)/nbins
        self.r = self.intervals[:-1]+self.dr/2.0
        ff = 4.0 / 3.0 * np.pi * (self.intervals[1:]**3 - self.intervals[:-1]**3)
        self.coordination_number = np.zeros(nbins)
        self.rdf = np.zeros((nbins), dtype=np.double)
        for structure in structures:
            lattice = structure.lattice
            rho = float(len(indices_i)) / lattice.volume
            i_frac_coords = structure.frac_coords[indices_i]
            j_frac_coords = structure.frac_coords[indices_j]
            dr_ij = lattice.get_all_distances(i_frac_coords, j_frac_coords)
            mask = np.ones(dr_ij.shape, dtype=bool)
            if self_reference:
                np.fill_diagonal(mask, 0)
            dr_ij = np.ndarray.flatten(dr_ij[mask])
            hist = np.histogram(dr_ij, bins=nbins, range=(r_min, r_max), density=False)[0]
            self.rdf += hist / rho
            self.coordination_number += np.cumsum(hist)
        self.rdf = self.rdf / ff / len(structures) / float(len(indices_j))
        self.coordination_number = self.coordination_number / len(structures) / float(len(indices_j))
        
    def smeared_rdf(self,sigma=0.1):
        """
        Smear the RDF with a Gaussian kernel.

        Args:
            sigma (:obj:`float`, optional): Standard deviation for Gaussian kernel. Optional, default is 0.1.

        Returns:
            (np.array): Smeared RDF data.

        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.rdf, sigma=sigma_n_bins)

    
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
    
    def __init__(self, structures, indices, d_steps, nbins=500, r_min=0.0, r_max=10.0):
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
        self.intervals = np.linspace(r_min, r_max, nbins+1)
        self.dr = (r_max - r_min)/nbins
        self.r = self.intervals[:-1]+self.dr/2.0
        self.gdrt = np.zeros((nbins), dtype=np.double)
        self.gsrt = np.zeros((nbins), dtype=np.double)
        rho = len(indices) / structures[0].lattice.volume
        lattice = structures[0].lattice
        ff = 4.0 / 3.0 * np.pi * (self.intervals[1:]**3 - self.intervals[:-1]**3)
        rho = len(indices) / lattice.volume
        for struc_i, struc_j in zip( structures[:len(structures)-d_steps], structures[d_steps:]):
            i_frac_coords = struc_i.frac_coords[indices]
            j_frac_coords = struc_j.frac_coords[indices]
            dr_ij = lattice.get_all_distances(i_frac_coords, j_frac_coords)
            mask = np.ones(dr_ij.shape, dtype=bool)
            np.fill_diagonal(mask, 0)
            distinct_dr_ij = np.ndarray.flatten(dr_ij[mask])
            hist = np.histogram(distinct_dr_ij, bins=nbins, range=(0.0, r_max), density=False)[0]
            self.gdrt += hist / rho
            self_dr_ij = np.ndarray.flatten(dr_ij[np.invert(mask)])
            hist = np.histogram(self_dr_ij, bins=nbins, range=(0.0, r_max), density=False)[0]
            self.gsrt += hist / rho
        self.gdrt = self.gdrt / ff / (len(structures)-d_steps) / float(len(indices))
        self.gsrt = self.gsrt / (len(structures)-d_steps) / float(len(indices))        
       
    def self(self, sigma=None):
        if sigma:
            return self.smeared_gsrt(sigma=sigma)
        else:
            return self.gsrt

    def distinct(self, sigma=None):
        if sigma:
            return self.smeared_gdrt(sigma=sigma)
        else:
            return self.gdrt
 
    def smeared_gsrt(self,sigma=0.1):
        """
        Smear the self part of the Van Hove correlation function with a Gaussian kernel.

        Args:
            sigma (:obj:`float`, optional): Standard deviation for Gaussian kernel. Optional, default is 0.1.

        Returns:
            (np.array): Smeared data.

        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.gsrt, sigma=sigma_n_bins)
    
    def smeared_gdrt(self,sigma=0.1):
        """
        Smear the distinct part of the Van Hove correlation function with a Gaussian kernel.

        Args:
            sigma (:obj:`float`, optional): Standard deviation for Gaussian kernel. Optional, default is 0.1.
        
        Returns:
            (np.array): Smeared data.

        """
        sigma_n_bins = sigma / self.dr
        return gaussian_filter1d(self.gdrt, sigma=sigma_n_bins)
