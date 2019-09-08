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

class Rdf:
    """
    class for radial distribution functions
    """

    def __init__( self, max_r, number_of_bins ):
        """
        Initialise a Rdf object for manipulating radial distribution functions.
        
        Args:
            max_r (Float): the maximum r value stored for g(r).
            number_of_bins (Int): number of bins for storing data about g(r).

        Returns:
            None
        """
        self.max_r = max_r
        self.number_of_bins = number_of_bins
        self.data = np.zeros( number_of_bins )
        self.dr = max_r / number_of_bins
 
    def add_dr( self, dr ):
        """
        Add an observed interatomic distance to the g(r) data at dr.

        Args:
            dr (Float): the interatomic distance, dr.

        Returns:
            None
        """ 
        this_bin = int( dr / self.dr ) 
        if this_bin > self.number_of_bins:
            raise IndexError( 'dr is larger than rdf max_r' )
        self.data[ this_bin ] += 1

    def normalised_data( self ):
        return( np.array( [ [ dr, g_of_r / self.volume_factor( dr ) ] for dr, g_of_r in zip( np.array( [ self.dr * ( r + 0.5 ) for r in range( 0, self.number_of_bins ) ] ), self.data / np.sum( self.data ) ) ] ) )  

    def __add__( self, other_rdf ):
        assert isinstance( other_rdf, Rdf )
        assert self.max_r == other_rdf.max_r
        assert self.number_of_bins == other_rdf.number_of_bins
        summed_rdf = Rdf( self.max_r, self.number_of_bins )
        summed_rdf.data = self.data + other_rdf.data
        return summed_rdf 

    def volume_factor( self, dr ):
        return ( ( dr + self.dr )**3 - dr**3 )
