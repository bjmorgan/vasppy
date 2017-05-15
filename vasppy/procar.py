import numpy as np
import re
import math

angstrom_to_bohr = 0.52918
ev_to_hartree = 0.036749309

def get_numbers_from_string( string ):
    p = re.compile('-?\d+[.\d]*')
    return [ float( s ) for s in p.findall( string ) ] 

def k_point_parser( string ):
    regex = re.compile( 'k-point\s+\d+\s*:\s+((?:[- ][01].\d{8}){3})' )
    return [ [ float(s) for s in [ x[0:11], x[11:22], x[22:33] ] ] for x in regex.findall( string ) ] 

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

    def __init__( self, spin = 1 ):
        self.spin_channels = spin # should be determined from PROCAR
        self.number_of_k_points = None
        self.number_of_ions = None
        self.number_of_bands = None
        self.data = None
        self.bands = None
        self.number_of_projections = None
        self.k_point_blocks = None
        self.calculation = { 'non_spin_polarised': False, 'non_collinear': False, 'spin_polarised': False }
        self.non_spin_polarised = None

    def parse_projections( self ):
        self.projection_data = projections_parser( self.read_in )
        try:
            assert( self.number_of_bands * self.number_of_k_points == len( self.projection_data ) )
            self.spin_channels = 1 # non-magnetic, non-spin-polarised
            self.k_point_blocks = 1
            self.calculation[ 'non_spin_polarised' ] = True
        except:
            if self.number_of_bands * self.number_of_k_points * 4 == len( self.projection_data ):
                self.spin_channels = 4 # non-collinear (spin-orbit coupling)
                self.k_point_blocks = 1
                self.calculation[ 'non_collinear' ] = True
                pass
            elif self.number_of_bands * self.number_of_k_points * 2 == len( self.projection_data ):
                self.spin_channels = 2 # spin-polarised
                self.k_point_blocks = 2
                self.calculation[ 'spin_polarised' ] = True
                pass
            else:
                raise
        self.number_of_projections = int( self.projection_data.shape[1] / ( self.number_of_ions + 1 ) )

    def parse_k_points( self ):
        k_points = k_point_parser( self.read_in )
        self.k_points = np.array( k_points, dtype = float )

    def parse_bands( self ):
        bands = re.findall( r"band\s*(\d+)\s*#\s*energy\s*([-.\d\s]+)", self.read_in )
        self.bands = np.array( bands, dtype = float )

    def sanity_check( self ):
        expected_k_points = self.number_of_k_points
        read_k_points = len( self.k_points ) / self.k_point_blocks
        assert( expected_k_points == read_k_points ), "k-point number mismatch: {} in header; {} in file".format( expected_k_points, read_k_points )
        expected_bands = self.number_of_bands
        read_bands = len( self.bands ) / self.number_of_k_points / self.k_point_blocks
        assert( expected_bands == read_bands ), "band mismatch: {} in header; {} in file".format( expected_bands, read_bands )
 
    def read_from_file( self, filename, bands_in_range = None ):
        with open( filename, 'r' ) as file_in:
            file_in.readline()
            self.number_of_k_points, self.number_of_bands, self.number_of_ions = [ int( f ) for f in get_numbers_from_string( file_in.readline() ) ]
            self.read_in = file_in.read()
        self.parse_k_points()
        self.parse_bands()
        self.parse_projections()
        self.sanity_check()
        self.read_in = None
        if self.calculation[ 'spin_polarised' ]:
            self.data = self.projection_data.reshape( self.spin_channels, self.number_of_k_points, self.number_of_bands, self.number_of_ions + 1, self.number_of_projections )[:,:,:,:,1:].swapaxes( 0, 1).swapaxes( 1, 2 )
        else:
            self.data = self.projection_data.reshape( self.number_of_k_points, self.number_of_bands, self.spin_channels, self.number_of_ions + 1, self.number_of_projections )[:,:,:,:,1:]

    def total_band_structure( self, spin ):
        # note: currently gives k-points linear spacing
        # if we know the k-vectors for each k-point can instead use their geometric separations to give the correct k-point density
        # note: correct k-spacing is already implemented in weighted_band_structure
        assert( self.bands.shape == ( self.number_of_bands * self.number_of_k_points, 2 ) )
        to_return = np.insert( band_energies, 0, range( 1, self.number_of_k_points + 1 ), axis = 1 )
        return to_return 

    def print_weighted_band_structure( self, spins = None, ions = None, orbitals = None, scaling = 1.0, e_fermi = 0.0, reciprocal_lattice = None ):
        if spins:
            spins = [ s - 1 for s in spins ]
        else:
            spins = list( range( self.spin_channels ) )
        if not ions:
            ions = [ self.number_of_ions ] 
        if not orbitals:
            orbitals = [ self.data.shape[-1]-1 ] # !! NOT TESTED YET FOR f STATES !!
        if self.calculation[ 'spin_polarised' ]:
            band_energies = self.bands[:,1:].reshape( self.spin_channels, self.number_of_k_points, self.number_of_bands )[ spins[0] ].T
        else:
            band_energies = self.bands[:,1:].reshape( self.number_of_k_points, self.number_of_bands ).T
        orbital_projection = np.sum( self.data[ :, :, :, :, orbitals ], axis = 4 )
        ion_projection = np.sum( orbital_projection[ :, :, :, ions ], axis = 3 ) 
        spin_projection = np.sum( ion_projection[ :, :, spins ], axis = 2 )
        x_axis = self.x_axis( reciprocal_lattice )
        for i in range( self.number_of_bands ):
            print( '# band: {}'.format( i + 1 ) )
            for k, ( e, p ) in enumerate( zip( band_energies[i], spin_projection.T[i] ) ):
                print( x_axis[ k ], e - e_fermi, p * scaling ) # k is the k_point index: currently gives linear k-point spacing
            print()

    def effective_mass_calc( self, k_point_indices, band_index, reciprocal_lattice, spin = 1, printing = False ):
        assert( spin <= self.k_point_blocks )
        assert( len( k_point_indices ) > 1 ) # we need at least 2 k-points
        band_energies = self.bands[:,1:].reshape( self.k_point_blocks, self.number_of_k_points, self.number_of_bands )
        k_points = np.array( [ self.k_points[ k - 1 ] for k in k_point_indices ] )
        eigenvalues = np.array( [ band_energies[ spin - 1 ][ k - 1 ][ band_index - 1 ] for k in k_point_indices ] )
        if printing:
            print( '# h k l e' )
            [ print( ' '.join( [ str( f ) for f in row ] ) ) for row in np.concatenate( ( k_points, np.array( [ eigenvalues ] ).T ), axis = 1 ) ]
        reciprocal_lattice = reciprocal_lattice * 2 * math.pi * angstrom_to_bohr
        cartesian_k_points = np.array( [ np.dot( k, reciprocal_lattice ) for k in k_points ] ) # convert k-points to cartesian
        if len( k_point_indices ) == 2:
            effective_mass_function = two_point_effective_mass
        else:
            effective_mass_function = least_squares_effective_mass
        return effective_mass_function( cartesian_k_points, eigenvalues )

    def x_axis( self, reciprocal_lattice ):
        if reciprocal_lattice is not None:
            cartesian_k_points = np.dot( self.k_points, reciprocal_lattice )
            x_axis = [ 0.0 ]
            for i in range( 1, len( cartesian_k_points ) ):
                dk = cartesian_k_points[ i - 1 ] - cartesian_k_points[ i ]
                mod_dk = np.sqrt( np.dot( dk, dk ) )
                x_axis.append( mod_dk + x_axis[-1] )
            x_axis = np.array( x_axis )
        else:
            x_axis = np.arange( len( self.k_points ) )
        return x_axis
