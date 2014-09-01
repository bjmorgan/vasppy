import numpy as np
import re

def get_numbers_from_string( string ):
    p = re.compile('-?\d+[.\d]*')
    return( [ float( s ) for s in p.findall( string ) ] )

class Occupation:

    def __init__( self, s, p_x, p_y, p_z, d_xy, d_yz, d_z2, d_xz, d_x2, total ):
        self.s = s
        self.p = [ p_x, p_y, p_z ]
        self.d = [ d_xy, d_yz, d_z2, d_xz, d_x2 ]
        self.total = total

class Ion:

    def __init__( self, number, occupations ):
        self.number      = number
        self.occupations = Occupation( *occupations )

class Band:

    def __init__( self, number, energy, occupation, ions, total ):
        self.number     = number
        self.energy     = energy
        self.occupation = occupation
        self.ions       = ions
        self.total      = Occupation( *total )

    def sorted_ions_by_occupation( self ):
        return( sorted( self.ions, key = lambda x: x.occupations.total ) )

class KPoint:

    def __init__( self, number, k_vector, weight, bands ):
        self.number   = number
        self.k_vector = k_vector
        self.weight   = weight
        self.bands    = bands

    def bands_by_number( self, numbers ):
        return( [ band for band in self.bands if band.number in numbers ] )

class Spin:

    def __init__( self, k_points ):
        self.k_points = k_points

class Procar:

    number_of_file_header_lines = 1
    number_of_spin_header_lines = 1
    number_of_kpoint_header_lines = 2
    number_of_band_header_lines = 3

    def __init__( self, spin = 1 ):
        self.spin_channels = spin

    def read_from_file( self, filename, bands_in_range = None ):
        with open( filename, 'r' ) as file_in:
            self.lines = [ line.strip() for line in file_in.readlines() ]
            self.lines.pop( Procar.number_of_file_header_lines - 1 )
            ( self.number_of_k_points, self.number_of_bands, self.number_of_ions ) = [ int( f ) for f in get_numbers_from_string( self.lines[0] ) ] 
            if bands_in_range == None:
                self.bands_in_range = [ 1, self.number_of_bands ]
            else:
                self.bands_in_range = bands_in_range # bands in range should be a list of the form ( start_number, end_number )
            print( self.number_of_k_points )
            self.spin = [ self.read_spin() for spin in range( self.spin_channels ) ]

    def next_line( self ):
        return( self.lines.pop( 0 ) )

    def pop_lines( self, i ):
        self.lines = self.lines[i:]
        # del( self.lines[ 0:i ] )

    def read_spin( self ):
        spin_data = self.next_line()
        print( "spin_data", spin_data )
        self.pop_lines( 1 )
        k_points = [ self.read_k_point() for k_point in range( self.number_of_k_points ) ]
        return( Spin( k_points ) )

    def read_k_point( self ):
        k_point_data = get_numbers_from_string( self.next_line() )
        number   = int( k_point_data[0] )
        k_vector = k_point_data[1:4]
        weight   = k_point_data[4]
        print( 'k_point data', number, k_vector, weight )
        self.pop_lines( 1 )
        bands = [ self.read_band_data() for band in range( self.number_of_bands ) ]
        bands = [ band for band in bands if band is not None ]
        self.pop_lines( 1 )
        return( KPoint( number, k_vector, weight, bands ) )

    def read_band_data( self ):
        ( number, energy, occupation ) = get_numbers_from_string( self.next_line() )
        number = int( number )
        if not self.bands_in_range or number >= self.bands_in_range[0] and number <= self.bands_in_range[1]:
            print( 'band_data', number, energy, occupation )
            self.pop_lines( 2 )
            ions  = [ self.read_ion_data() for ion in range( self.number_of_ions ) ] # read data for ion occupations
            total = [ float( s ) for s in self.next_line().split()[1:] ] # read band total data
            self.pop_lines( 1 )
            return( Band( number, energy, occupation, ions, total ) )
        else:
            self.pop_lines( self.number_of_ions + 4 )
            return( None )

    def read_ion_data( self ):
        ion_data = self.next_line().split()
        number = int( ion_data[0] )
        occupations = [ float( s ) for s in ion_data[1:] ]
        return( Ion( number, occupations ) )


