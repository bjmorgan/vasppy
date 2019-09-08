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
    
    def __init__(self, structures, indices_i, indices_j, nbins=500, r_min=0.0, r_max=10.0):
        """
        Initialise a RadialDistributionFunction instance.

        Args:
            structures (list(pymatgen.Structure)): List of pymatgen Structure objects.
            indices_i (list(int)): List of indices for species i.
            indices_j (list(int)): List of indices for species j.
            nbins (:obj:`int`, optional): Number of bins used for the RDF, Optional, default is 500.
            rmin (:obj:`floet`, optional): Minimum r value. Optional, default is 0.0.
            rmax (:obj:`floet`, optional): Maximum r value. Optional, default is 10.0.

        Returns:
             None

        """
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
            dr_ij = np.ndarray.flatten(lattice.get_all_distances(i_frac_coords, j_frac_coords))
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
