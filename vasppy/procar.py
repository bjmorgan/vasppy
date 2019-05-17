import numpy as np
import re
import math
import warnings
from .units import angstrom_to_bohr, ev_to_hartree
from copy import deepcopy
import fortranformat as ff
from functools import reduce

class KPoint():

    def __init__( self, index, frac_coords, weight ):
        self.index = index
        self.frac_coords = frac_coords
        self.weight = weight

    def cart_coords( self, reciprocal_lattice ):
        """Convert the reciprocal fractional coordinates for this k-point to \
        reciprocal Cartesian coordinates.

        Args:
            reciprocal_lattice (np.array(float)): 3x3 numpy array containing the \
                Cartesian reciprocal lattice.

        Returns:
            (np.array): The reciprocal Cartesian coordinates of this k-point.

        """
        return np.dot( self.frac_coords, reciprocal_lattice )

    def __eq__( self, other ):
        return ( ( self.index == other.index ) &
                 ( self.frac_coords == other.frac_coords ).all() &
                 ( self.weight == other.weight ) )

    def __repr__( self ):
        return "k-point {}: {} weight = {}".format( self.index, ' '.join( [ str(c) for c in self.frac_coords ] ), self.weight )

def handle_occupancy( occupancy, negative_occupancies='warn' ):
    valid_negative_occupancies = [ 'warn', 'raise', 'ignore', 'zero' ]
    if negative_occupancies not in valid_negative_occupancies:
        raise ValueError( "valid options for negative_occupancies are {}".format( valid_negative_occupancies ) )
    if occupancy < 0:
       if negative_occupancies == 'warn':
           warnings.warn( "One or more occupancies in your PROCAR file are negative." )
       elif negative_occupancies == 'raise':
           raise ValueError( "One or more occupancies in your PROCAR file are negative." )
       elif negative_occupancies == 'ignore':
           pass
       elif negative_occupancies == 'zero':
           occupancy = 0.0
    return occupancy

class Band():

    def __init__( self, index, energy, occupancy, negative_occupancies='warn' ):
        self.index = index
        self.energy = energy
        self.occupancy = handle_occupancy( occupancy, negative_occupancies=negative_occupancies )

    def __eq__( self, other ):
        return ( ( self.index == other.index ) & 
                 ( self.energy == other.energy ) & 
                 ( self.occupancy == other.occupancy ) )

def get_numbers_from_string( string ):
    p = re.compile('-?\d+[.\d]*')
    return [ float( s ) for s in p.findall( string ) ]

def k_point_parser( string ):
    """Parse k-point data from a PROCAR string.

    Finds all lines of the form::

         k-point    1 :    0.50000000 0.25000000 0.75000000     weight = 0.00806452

    and extracts the k-point index, reciprocal fractional coordinates, and weight
    into a :obj:`procar.KPoint` object.

    Args:
        string (str): String containing a full PROCAR file.

    Returns:
        (list(:obj:`procar.KPoint`)): A list of :obj:`procar.KPoint` objects.

    """
    regex = re.compile( 'k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)' )
    captured = regex.findall( string ) 
    k_points = []
    for kp in captured:
        index = int(kp[0])
        frac_coords = np.array( [ float(s) for s in kp[1:4] ] )
        weight = float(kp[4])
        k_points.append( KPoint( index=index, frac_coords=frac_coords, weight=weight ) )
    return k_points

def projections_parser( string ):
    regex = re.compile( '([-.\d\se]+tot.+)\n' )
    data = regex.findall( string )
    data = [ x.replace( 'tot', '0' ) for x in data ]
    data = np.array( [ x.split() for x in data ], dtype = float )
    return data

def area_of_a_triangle_in_cartesian_space( a, b, c ):
    """
    Returns the area of a triangle defined by three points in Cartesian space.

    Args:
        a (np.array): Cartesian coordinates of point A.
        b (np.array): Cartesian coordinates of point B.
        c (np.array): Cartesian coordinates of point C.

    Returns:
        (float): the area of the triangle.
    """
    return 0.5 * np.linalg.norm( np.cross( b-a, c-a ) )

def points_are_in_a_straight_line( points, tolerance=1e-7 ):
    """
    Check whether a set of points fall on a straight line.
    Calculates the areas of triangles formed by triplets of the points.
    Returns False is any of these areas are larger than the tolerance.

    Args:
        points (list(np.array)): list of Cartesian coordinates for each point.
        tolerance (optional:float): the maximum triangle size for these points to be considered colinear. Default is 1e-7.

    Returns:
        (bool): True if all points fall on a straight line (within the allowed tolerance).
    """
    a = points[0]
    b = points[1]
    for c in points[2:]:
        if area_of_a_triangle_in_cartesian_space( a, b, c ) > tolerance:
            return False
    return True

def two_point_effective_mass( cartesian_k_points, eigenvalues ):
    """
    Calculate the effective mass given eigenvalues at two k-points.
    Reimplemented from Aron Walsh's original effective mass Fortran code.

    Args:
        cartesian_k_points (np.array): 2D numpy array containing the k-points in (reciprocal) Cartesian coordinates.
        eigenvalues (np.array):        numpy array containing the eigenvalues at each k-point.

    Returns:
        (float): The effective mass
    """
    assert( cartesian_k_points.shape[0] == 2 )
    assert( eigenvalues.size == 2 )
    dk = cartesian_k_points[ 1 ] - cartesian_k_points[ 0 ]
    mod_dk = np.sqrt( np.dot( dk, dk ) )
    delta_e = ( eigenvalues[ 1 ] - eigenvalues[ 0 ] ) * ev_to_hartree * 2.0
    effective_mass = mod_dk * mod_dk / delta_e
    return effective_mass

def least_squares_effective_mass( cartesian_k_points, eigenvalues ):
    """
    Calculate the effective mass using a least squares quadratic fit.

    Args:
        cartesian_k_points (np.array): Cartesian reciprocal coordinates for the k-points
        eigenvalues (np.array):        Energy eigenvalues at each k-point to be used in the fit.

    Returns:
        (float): The fitted effective mass

    Notes:
        If the k-points do not sit on a straight line a ValueError will be raised.
    """
    if not points_are_in_a_straight_line( cartesian_k_points ):
        raise ValueError( 'k-points are not collinear' )
    dk = cartesian_k_points - cartesian_k_points[0]
    mod_dk = np.linalg.norm( dk, axis = 1 )
    delta_e = eigenvalues - eigenvalues[0]
    effective_mass = 1.0 / ( np.polyfit( mod_dk, eigenvalues, 2 )[0] * ev_to_hartree * 2.0 )
    return effective_mass

class Procar:
    """
    Object for working with PROCAR data.

    Attributes:
        data (numpy.array(float)): A 5D numpy array that stores the projection data.

                    Axes are k-points, bands, spin-channels, ions and sum over ions, lm-projections.

        bands (numpy.array(:obj:`Band`)): A numpy array of ``Band`` objects, that contain band index, energy, and occupancy data.
        k_points (numpy.array(:obj:`KPoint`)): A numpy array of ``KPoint`` objects, that contain fractional coordinates and weights for each k-point.
        number_of_k_points (int): The number of k-points.
        number_of_bands (int): The number of bands.
        spin_channels (int): Number of spin channels in the PROCAR data:
 
                    -  1 for non-spin-polarised calculations.
                    -  2 for spin-polarised calculations.
                    -  4 for non-collinear calculations.
        number_of_ions (int): The number of ions.
        number_of_projections (int): The number of projections, e.g. TODO
        calculation (dict): Dictionary of True | False values describing the calculation type.

            Dictionary keys are 'non_spin_polarised', 'non_collinear', and 'spin_polarised'.
      
    """

    def __init__( self, spin=1, negative_occupancies='warn' ):
        self._spin_channels = spin # should be determined from PROCAR
        self._number_of_k_points = None
        self._number_of_ions = None
        self._number_of_bands = None
        self._number_of_projections = None
        self._k_point_blocks = None
        self._data = None
        self._bands = None
        self._k_points = None
        self.calculation = { 'non_spin_polarised': False, 'non_collinear': False, 'spin_polarised': False }
        if negative_occupancies not in [ 'warn', 'raise', 'zero' ]:
            raise ValueError( "negative_occupancies can be one of [ 'warn', 'raise', 'zero' ]" )
        self.negative_occupancies = negative_occupancies
        #self.non_spin_polarised = None

    @property
    def occupancy( self ):
        return np.array( [ [ band.index, band.occupancy ] for band in self._bands ] )

    def __add__( self, other ):
        if self.spin_channels != other.spin_channels:
            raise ValueError( 'Can only concatenate Procars with equal spin_channels: {}, {}'.format( self.spin_channels, other.spin_channels ) )
        if self.number_of_ions != other.number_of_ions:
            raise ValueError( 'Can only concatenate Procars with equal number_of_ions: {}, {}'.format( self.number_of_ions, other.number_of_ions ) )
        if self.number_of_bands != other.number_of_bands:
            raise ValueError( 'Can only concatenate Procars with equal number_of_bands: {}, {}'.format( self.number_of_bands, other.number_of_bands ) )
        if self.number_of_projections != other.number_of_projections:
            raise ValueError( 'Can only concatenate Procars with equal number_of_projections: {}, {}'.format( self.number_of_projections, other.number_of_projections ) )
        if self._k_point_blocks != other._k_point_blocks:
            raise ValueError( 'Can only concatenate Procars with equal k_point_blocks: {}, {}'.format( self._k_point_blocks, other._k_point_blocks ) )
        if self.calculation != other.calculation:
            raise ValueError( 'Can only concatenate Procars from equal calculations: {}, {}'.format( self.calculation, other.calculation ) )
        new_procar = deepcopy( self )
        new_procar._data = np.concatenate( ( self._data, other._data ), axis=0 )
        new_procar._number_of_k_points = self.number_of_k_points + other.number_of_k_points
        new_procar._bands = []
        new_procar._bands = np.ravel( np.concatenate( [ self.bands, other.bands ], axis=1 ) )
        new_procar._k_points = self._k_points + other._k_points
        for i, kp in enumerate( new_procar._k_points, 1 ):
            kp.index = i    
        new_procar.sanity_check()
        return new_procar
 
    def parse_projections( self ):
        self.projection_data = projections_parser( self.read_in )
        try:
            assert( self._number_of_bands * self._number_of_k_points == len( self.projection_data ) )
            self._spin_channels = 1 # non-magnetic, non-spin-polarised
            self._k_point_blocks = 1
            self.calculation[ 'non_spin_polarised' ] = True
        except:
            if self._number_of_bands * self._number_of_k_points * 4 == len( self.projection_data ):
                self._spin_channels = 4 # non-collinear (spin-orbit coupling)
                self._k_point_blocks = 1
                self.calculation[ 'non_collinear' ] = True
                pass
            elif self._number_of_bands * self._number_of_k_points * 2 == len( self.projection_data ):
                self._spin_channels = 2 # spin-polarised
                self._k_point_blocks = 2
                self.calculation[ 'spin_polarised' ] = True
                pass
            else:
                raise
        self._number_of_projections = int( self.projection_data.shape[1] / ( self._number_of_ions + 1 ) ) - 1

    def parse_k_points( self ):
        self._k_points = k_point_parser( self.read_in )[:self._number_of_k_points]

    def parse_bands( self ):
        band_data = re.findall( r"band\s*(\d+)\s*#\s*energy\s*([-.\d]+)\s?\s*#\s"r"*occ.\s*([-.\d]+)", self.read_in )
        self._bands = np.array( [ Band( float(i), float(e), float(o), negative_occupancies=self.negative_occupancies ) for i, e, o in band_data ] )

    def sanity_check( self ):
        assert( self._number_of_k_points == len( self._k_points ) ), "k-point number mismatch: {} in header; {} in file".format( self._number_of_k_points, len( self._k_points ) )
        read_bands = len( self._bands ) / self._number_of_k_points / self._k_point_blocks
        assert( self._number_of_bands == read_bands ), "band mismatch: {} in header; {} in file".format( self._number_of_bands, read_bands )

    @classmethod
    def from_files( cls, filenames, **kwargs ):
        """Create a :obj:`Procar` object by reading the projected wavefunction character of each band
        from a series of VASP ``PROCAR`` files.

        Useful when e.g. a band-structure calculation has been split over multiple VASP calculations,
        for example, when using hybrid functionals.

        Args:
            filename (str): Filename of the ``PROCAR`` file.
            **kwargs: See the ``from_file()`` method for a description of keyword arguments.

        Returns:
            (:obj:`vasppy.Procar`)
        
        """
        pcars = [ cls.from_file( f, **kwargs ) for f in filenames ]
        return reduce( cls.__add__, pcars )

    @classmethod
    def from_file( cls, filename, negative_occupancies='warn',
                   select_zero_weighted_k_points=False ):
        """Create a :obj:`Procar` object by reading the projected wavefunction character of each band
        from a VASP ``PROCAR`` file.

        Args:
            filename (str): Filename of the ``PROCAR`` file.
            negative_occupancies (:obj:`Str`, optional): Select how negative occupancies are handled.
                Options are:

                    - ``warn`` (default): Warn that some partial occupancies are negative.
                    - ``raise``:          Raise an AttributeError.
                    - ``ignore``:         Do nothing.
                    - ``zero``:           Negative partial occupancies will be set to zero.
            select_zero_weighted_k_points (:obj:`bool`, optional): Set to ``True`` to only 
                read in zero-weighted k-points from the ``PROCAR`` file. Default is ``False``.

        Returns:
            (:obj:`vasppy.Procar`)
        
        """
        pcar = cls( negative_occupancies=negative_occupancies )
        pcar._read_from_file( filename=filename )
        if select_zero_weighted_k_points:
            k_point_indices = [ i for i, kp in enumerate( pcar.k_points ) if kp.weight == 0.0 ]
            pcar = pcar.select_k_points( k_point_indices )
        return pcar
       
    def read_from_file( self, filename ):
        warnings.warn( "read_from_file() is deprecated as a part of the public API.\nPlease use Procar.from_file() or Procar.from_files() instead" )
        return self._read_from_file( filename=filename )
 
    def _read_from_file( self, filename ):
        """Reads the projected wavefunction character of each band from a VASP PROCAR file.

        Args:
            filename (str): Filename of the PROCAR file.

        Returns:
            None
        
        """
        with open( filename, 'r' ) as file_in:
            file_in.readline()
            self._number_of_k_points, self._number_of_bands, self._number_of_ions = [ int( f ) for f in get_numbers_from_string( file_in.readline() ) ]
            self.read_in = file_in.read()
        self.parse_k_points()
        self.parse_bands()
        self.parse_projections()
        self.sanity_check()
        self.read_in = None # clear memory
        if self.calculation[ 'spin_polarised' ]:
            self._data = self.projection_data.reshape( self._spin_channels, self._number_of_k_points, self._number_of_bands, self._number_of_ions+1, self._number_of_projections+1 )[:,:,:,:,1:].swapaxes( 0, 1).swapaxes( 1, 2 )
        else:
            self._data = self.projection_data.reshape( self._number_of_k_points, self._number_of_bands, self._spin_channels, self._number_of_ions+1, self._number_of_projections+1 )[:,:,:,:,1:]

    @property
    def number_of_k_points( self ):
        """The number of k-points described by this :obj:`Procar` object."""
        assert( self._number_of_k_points == self._data.shape[0] ), "Number of k-points in metadata ({}) not equal to number in PROCAR data ({})".format( self._number_of_k_points, self._data.shape[0] )
        return self._number_of_k_points

    @property
    def number_of_bands( self ):
        """The number of bands described by this :obj:`Procar` object."""
        assert( self._number_of_bands == self._data.shape[1] ), "Number of bands in metadata ({}) not equal to number in PROCAR data ({})".format( self._number_of_bands, self._data.shape[1] )
        return self._number_of_bands

    @property
    def spin_channels( self ):
        """The number of spin-channels described by this :obj:`Procar` object."""
        assert( self._spin_channels == self._data.shape[2] ), "Number of spin channels in metadata ({}) not equal to number in PROCAR data ({})".format( self._spin_channels, self._data.shape[2] )
        return self._spin_channels

    @property
    def number_of_ions( self ):
        """The number of ions described by thie :obj:`Procar` object."""
        assert( self._number_of_ions == self._data.shape[3]-1 ), "Number of ions in metadata ({}) not equal to number in PROCAR data ({})".format( self._number_of_ions, self._data.shape[3]-1 )
        return self._number_of_ions

    @property
    def number_of_projections( self ):
        """The number of lm-projections described by this :obj:`Procar` object."""
        assert (self._number_of_projections == self._data.shape[4]), "Number of projections in metadata ({}) not equal to number in PROCAR data ({})".format( self._number_of_projections, self._data.shape[4] ) 
        return self._number_of_projections

    def print_weighted_band_structure( self, spins=None, ions=None, orbitals=None, scaling=1.0, e_fermi=0.0, reciprocal_lattice=None ):
        band_structure_data = self.weighted_band_structure( spins=spins, ions=ions, orbitals=orbitals, scaling=scaling, e_fermi=e_fermi, reciprocal_lattice=reciprocal_lattice ) 
        for i, band_data in enumerate( band_structure_data, 1 ):
            print( '# band: {}'.format( i ) )
            for k_point_data in band_data:
                print( ' '.join( [ str(f) for f in k_point_data ] ) )
            print()

    def weighted_band_structure( self, spins=None, ions=None, orbitals=None, scaling=1.0, e_fermi=0.0, reciprocal_lattice=None ):
        if spins:
            spins = [ s - 1 for s in spins ]
        else:
            spins = list( range( self.spin_channels ) )
        if not ions:
            ions = [ self.number_of_ions ] # nions+1 is the `tot` index
        if not orbitals:
            orbitals = list( range( self.number_of_projections ) )
        if self.calculation[ 'spin_polarised' ]:
            band_energies = np.array( [ band.energy for band in self._bands ] ).reshape( self.spin_channels, self.number_of_k_points, self.number_of_bands )[ spins[0] ].T
        else:
            band_energies = np.array( [ band.energy for band in self._bands ] ).reshape( self.number_of_k_points, self.number_of_bands ).T
        orbital_projection = np.sum( self._data[ :, :, :, :, orbitals ], axis = 4 )
        ion_projection = np.sum( orbital_projection[ :, :, :, ions ], axis = 3 )
        spin_projection = np.sum( ion_projection[ :, :, spins ], axis = 2 )
        x_axis = self.x_axis( reciprocal_lattice )
        to_return = []
        for i in range( self.number_of_bands ):
            for k, ( e, p ) in enumerate( zip( band_energies[i], spin_projection.T[i] ) ):
                to_return.append( [ x_axis[ k ], e - e_fermi, p * scaling ] )
        to_return = np.array( to_return ).reshape( self.number_of_bands, -1, 3 )
        return to_return

    def effective_mass_calc( self, k_point_indices, band_index, reciprocal_lattice, spin=1, printing=False ):
        assert( spin <= self.k_point_blocks )
        assert( len( k_point_indices ) > 1 ) # we need at least 2 k-points
        band_energies = self._bands[:,1:].reshape( self.k_point_blocks, self.number_of_k_points, self.number_of_bands )
        frac_k_point_coords = np.array( [ self._k_points[ k - 1 ].frac_coords for k in k_point_indices ] )
        eigenvalues = np.array( [ band_energies[ spin - 1 ][ k - 1 ][ band_index - 1 ] for k in k_point_indices ] )
        if printing:
            print( '# h k l e' )
            [ print( ' '.join( [ str( f ) for f in row ] ) ) for row in np.concatenate( ( frac_k_point_coords, np.array( [ eigenvalues ] ).T ), axis = 1 ) ]
        reciprocal_lattice = reciprocal_lattice * 2 * math.pi * angstrom_to_bohr
        cart_k_point_coords = np.array( [ k.cart_coords( reciprocal_lattice ) for k in k_points ] ) # convert k-points to cartesian
        if len( k_point_indices ) == 2:
            effective_mass_function = two_point_effective_mass
        else:
            effective_mass_function = least_squares_effective_mass
        return effective_mass_function( cart_k_point_coords, eigenvalues )

    def x_axis( self, reciprocal_lattice=None ):
        """Generate the x-axis values for a band-structure plot.

        Returns an array of cumulative distances in reciprocal space between sequential k-points.

        Args:
            reciprocal_lattice (:obj:`np.array`, optional): 3x3 Cartesian reciprocal lattice.    
                Default is ``None``. If no reciprocal lattice is provided, the returned x-axis
                values will be sequential integers, giving even spacings between sequential
                k-points.
        
        Returns:
            (np.array): An array of x-axis values.
 
        """
        if reciprocal_lattice is not None:
            cartesian_k_points = np.array( [ k.cart_coords( reciprocal_lattice ) for k in k_points ] )
            x_axis = [ 0.0 ]
            for i in range( 1, len( cartesian_k_points ) ):
                dk = cartesian_k_points[ i - 1 ] - cartesian_k_points[ i ]
                mod_dk = np.sqrt( np.dot( dk, dk ) )
                x_axis.append( mod_dk + x_axis[-1] )
            x_axis = np.array( x_axis )
        else:
            x_axis = np.arange( len( self._k_points ) )
        return x_axis

    @property
    def bands( self ):
        return self._bands.reshape( self._k_point_blocks, 
                                    self._number_of_k_points, 
                                    self.number_of_bands )
       
    @property
    def k_points( self ):
        return self._k_points
 
    def select_bands_by_kpoint( self, band_indices ):
        return np.ravel( self.bands[:,band_indices,:] )

    def select_k_points( self, band_indices ):
        new_procar = deepcopy( self )
        new_procar._bands = np.ravel( new_procar.bands[:,band_indices,:] )
        new_procar._data = np.array( [ kp for i, kp in enumerate( new_procar._data ) if i in band_indices ] )
        new_procar._number_of_k_points = len( band_indices )
        new_procar._k_points = [ kp for i, kp in enumerate( new_procar._k_points ) if i in band_indices ]
        for i, kp in enumerate( new_procar._k_points, 1 ):
            kp.index = i
        new_procar.sanity_check()
        return new_procar

# TODO Need complete set of example PROCAR files before we start to write
# something that can write these all back out in the appropriate format
#    def write_file( self, filename ):
#        if self.calculation['non_collinear']:
#            raise NotImplementedError
#        with open( filename, 'w' ) as f:
#            f.write( 'PROCAR lm decomposed' )
#            for s in range( self.spin_channels ): # not sure what happens for non-collinear calculations
#                line = ff.FortranRecordWriter("/'# of k-points:',I5,9X,'# of bands:',I5,9X,'# of ions:',I5")
#                f.write( line.write( [ self.number_of_k_points, self.number_of_bands, self.number_of_ions ] )+'\n' )
#                for k_point in self.k_points:
#                    line = ff.FortranRecordWriter("/' k-point ',I5,' :',3X,3F11.8,'     weight = ',F10.8/")
#                    f.write( line.write( [ k_point.index, *k_point.frac_coords, k_point.weight ] ) )
#                    for band in self.bands:
#                        line = ff.FortranRecordWriter("/'band ',I5,' # energy',F14.8,' # occ.',F12.8/")
#                        f.write( line.write( [ band.index, band.energy, band.occupancy ] ) )
#                        f.write( '\nion      s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n' )
