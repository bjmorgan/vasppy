import numpy as np
import re
import sys

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
        self.data = self.projection_data.reshape( self.number_of_k_points, self.number_of_bands, self.spin_channels, self.number_of_ions + 1, 11 )[:,:,:,:,1:]

    def total_band_structure( self, spin ):
        # note: currently gives k-points linear spacing
        # if we know the k-vectors for each k-point can instead use their geometric separations to give the correct k-point density
        assert( self.bands.shape == ( self.number_of_bands * self.number_of_k_points, 2 ) )
        band_energies = self.bands[:,1:].reshape( self.number_of_k_points, self.number_of_bands )
        to_return = np.insert( band_energies, 0, range( 1, self.number_of_k_points + 1 ), axis = 1 )
        return to_return 

    def print_weighted_band_structure( self, spins, ions, orbitals, scaling = 1.0, e_fermi = 0.0 ):
        band_energies = self.bands[:,1:].reshape( self.number_of_k_points, self.number_of_bands ).T
        orbital_projection = np.sum( self.data[ :, :, :, :, orbitals ], axis = 4 )
        ion_projection = np.sum( orbital_projection[ :, :, :, ions ], axis = 3 ) 
        spin_projection = np.sum( ion_projection[ :, :, spins ], axis = 2 )
        for i in range( self.number_of_bands ):
            print( '# band: {}'.format( i + 1 ) )
            for k, ( e, p ) in enumerate( zip( band_energies[i], spin_projection.T[i] ) ):
                print( k, e - e_fermi, p * scaling ) # k is the k_point index: currently gives linear k-point spacing
            print()
