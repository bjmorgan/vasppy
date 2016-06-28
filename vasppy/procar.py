import numpy as np
import re
import math

angstrom_to_bohr = 0.52918
ev_to_hartree = 0.036749309

def get_numbers_from_string( string ):
    p = re.compile('-?\d+[.\d]*')
    return( [ float( s ) for s in p.findall( string ) ] )

class Procar:

    def __init__( self, spin = 1 ):
        self.spin_channels = spin # should be determined from PROCAR
        self.number_of_k_points = None
        self.number_of_ions = None
        self.number_of_bands = None
        self.data = None
        self.number_of_projections = None

    def parse_projections( self ):
        projection_data = re.findall( r"([-.\d\se]+tot.+)\n", self.read_in )
        projection_data = [ x.replace( 'tot', '0' ) for x in projection_data ]
        projection_data = [ x.split() for x in projection_data ]
        try:
            assert( self.number_of_bands * self.number_of_k_points == len( projection_data ) )
            self.spin_channels = 1 # non-magnetic, non-spin-polarised
        except:
            if self.number_of_bands * self.number_of_k_points * 4 == len( projection_data ):
                self.spin_channels = 4 # non-collinear (spin-orbit coupling)
                pass
            elif self.number_of_bands * self.number_of_k_points * 2 == len( projection_data ):
                self.spin_channels = 2 # spin-polarised
                pass
            else:
                raise
        self.projection_data = np.array( projection_data, dtype = float )
        self.number_of_projections = int( self.projection_data.shape[1] / ( self.number_of_ions + 1 ) )
        return( projection_data )

    def parse_k_points( self ):
        k_points = re.findall( r"k-point\s+\d+\s*:\s+([-.\d\s]+)", self.read_in )
        k_points = [ x.split() for x in k_points ]
        assert( self.number_of_k_points == len( k_points ) )
        self.k_points = np.array( k_points, dtype = float )
        return( k_points )

    def parse_bands( self ):
        bands = re.findall( r"band\s*(\d+)\s*#\s*energy\s*([-.\d\s]+)", self.read_in )
        assert( self.number_of_bands == len( bands ) / self.number_of_k_points )
        self.bands = np.array( bands, dtype = float )
        return( bands )
 
    def read_from_file( self, filename, bands_in_range = None ):
        with open( filename, 'r' ) as file_in:
            file_in.readline()
            self.number_of_k_points, self.number_of_bands, self.number_of_ions = [ int( f ) for f in get_numbers_from_string( file_in.readline() ) ]
            self.read_in = file_in.read()
        self.parse_k_points()
        self.parse_bands()
        self.parse_projections()
        self.read_in = None
        self.data = self.projection_data.reshape( self.number_of_k_points, self.number_of_bands, self.spin_channels, self.number_of_ions + 1, self.number_of_projections )[:,:,:,:,1:]

    def total_band_structure( self, spin ):
        # note: currently gives k-points linear spacing
        # if we know the k-vectors for each k-point can instead use their geometric separations to give the correct k-point density
        # note: correct k-spacing is already implemented in weighted_band_structure
        assert( self.bands.shape == ( self.number_of_bands * self.number_of_k_points, 2 ) )
        band_energies = self.bands[:,1:].reshape( self.number_of_k_points, self.number_of_bands )
        to_return = np.insert( band_energies, 0, range( 1, self.number_of_k_points + 1 ), axis = 1 )
        return to_return 

    def print_weighted_band_structure( self, spins = None, ions = None, orbitals = None, scaling = 1.0, e_fermi = 0.0, reciprocal_lattice = None ):
        if not spins:
            spins = list( range( self.spin_channels ) )
        if not ions:
            ions = [ self.number_of_ions ] 
        if not orbitals:
            orbitals = [ self.data.shape[-1]-1 ] # !! NOT TESTED YET FOR f STATES !!
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

    def effective_mass_calc( self, k_point_indices, band_index, reciprocal_lattice, printing = False ):
        assert( len( k_point_indices ) > 1 ) # we need at least 2 k-points
        band_energies = self.bands[:,1:].reshape( self.number_of_k_points, self.number_of_bands )
        k_points = np.array( [ self.k_points[ k - 1 ] for k in k_point_indices ] )
        eigenvalues = np.array( [ band_energies[ k - 1 ][ band_index - 1 ] for k in k_point_indices ] )
        if printing:
            print( '# h k l e' )
            [ print( ' '.join( [ str( f ) for f in row ] ) ) for row in np.concatenate( ( k_points, np.array( [ eigenvalues ] ).T ), axis = 1 ) ]
        reciprocal_lattice = reciprocal_lattice * 2 * math.pi * angstrom_to_bohr
        cartesian_k_points = np.array( [ np.dot( k, reciprocal_lattice ) for k in k_points ] ) # convert k-points to cartesian
        # TODO Check that the cartesian_k_points fall on a straight line, e.g. http://stackoverflow.com/questions/3813681/checking-to-see-if-3-points-are-on-the-same-line
        if len( k_point_indices ) == 2:
            # reimplemented from Aron's fortran version
            dk = cartesian_k_points[ 1 ] - cartesian_k_points[ 0 ]
            mod_dk = np.sqrt( np.dot( dk, dk ) )
            delta_e = ( eigenvalues[ 1 ] - eigenvalues[ 0 ] ) * ev_to_hartree * 2.0
            effective_mass = mod_dk * mod_dk / delta_e
        else:
            dk = cartesian_k_points - cartesian_k_points[0] 
            mod_dk = np.linalg.norm( dk, axis = 1 )
            delta_e = eigenvalues - eigenvalues[0]
            effective_mass = 1.0 / ( np.polyfit( mod_dk, eigenvalues, 2 )[0] * ev_to_hartree * 2.0 )
        return effective_mass

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
