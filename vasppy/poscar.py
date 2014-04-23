import numpy as np
import sys
import re
import copy

def angle( x, y ):
    dot = np.dot( x, y )
    x_mod = np.linalg.norm( x )
    y_mod = np.linalg.norm( y )
    cos_angle = dot / ( x_mod * y_mod )
    return np.degrees( np.arccos( cos_angle ) )

def parity( list ):
    return( sum( list )%2 )

class Poscar:

    lines_offset = 9

    def __init__( self ):
        self.title = "Title"
        self.scaling = 1.0
        self.lattice = np.identity( 3 )
        self.atoms = [ 'A' ]
        self.atom_numbers = [ 1 ]
        self.coordinate_type = 'Direct'
        self.coordinates = np.array( [ [ 0.0, 0.0, 0.0 ] ] )
  
    def read_from( self, filename ):
        try:
            with open( filename ) as f:
                lines = f.readlines()
        except FileNotFoundError:
            print( "\"" + filename + "\" not found", file=sys.stderr ) 
            sys.exit( -2 )
        self.title = lines.pop(0).strip()
        self.scaling = float( lines.pop(0).strip() )
        self.lattice = np.array( [ [ float( e ) for e in lines.pop(0).split() ] for i in range( 3 ) ] )
        self.atoms = lines.pop(0).split()
        self.atom_numbers = [ int(element) for element in lines.pop(0).split() ]
        self.coordinate_type = lines.pop(0)
        self.coordinates = np.array( [ [ float( e ) for e in lines.pop(0).split()[0:3] ] for i in range( sum( self.atom_numbers ) ) ] )
        if self.coords_are_cartesian(): # Convert to direct coordinates
            self.coordinates = self.fractional_coordinates()
            self.coordinate_type = 'Direct'

    def in_bohr( self ):
        new_poscar = copy.deepcopy( self )
        bohr_to_angstrom = 0.529177211
        new_poscar.scaling *= bohr_to_angstrom
        new_poscar.lattice /= bohr_to_angstrom
        if new_poscar.coords_are_cartesian():
            new_poscar.coordinates /= bohr_to_angstrom
        return( new_poscar )

    def coords_are_fractional( self ):
        return re.match( r'\A[Dd]', self.coordinate_type )

    def coords_are_cartesian( self ):
        return not self.coords_are_fractional()

    def fractional_coordinates( self ):
        return ( self.coordinates if self.coords_are_fractional() else self.coordinates.dot( np.linalg.inv( self.lattice ) ) )

    def cartesian_coordinates( self ):
        return ( self.coordinates if self.coords_are_cartesian() else self.coordinates.dot( self.lattice ) )

    def coordinates_to_stdout( self, coordinate_type='Direct' ):
        coord_opts = { 'Direct'    : self.fractional_coordinates(), 
                       'Cartesian' : self.cartesian_coordinates() }
        try:
            [ print( ''.join( ['  {: .10f}'.format( element ) for element in row ] ) ) for row in coord_opts[ coordinate_type ] ]
        except KeyError: 
            raise Exception( 'Passed coordinate_type: ' + coordinate_type + '\nAccepted values: [ Direct | Cartesian ] ' )

    def labelled_coordinates_to_stdout( self, coordinate_type='Direct', label_pos='1' ):
        coord_opts = { 'Direct'    : self.fractional_coordinates(), 
                       'Cartesian' : self.cartesian_coordinates() } 
        if label_pos == 1:
            for ( row, label ) in zip( coord_opts[ coordinate_type ], self.labels() ):
                print( label.ljust(6) + ''.join( [ '  {: .10f}'.format( element ) for element in row ] ) )
        if label_pos == 4:
            for ( row, label ) in zip( coord_opts[ coordinate_type ], self.labels() ):
                print( ''.join( [ '  {: .10f}'.format( element ) for element in row ] ) + '  {:.4s}'.format( label ) )

    def output_coordinates_only( self, coordinate_type='Direct', label=None ):
        if label:
            self.labelled_coordinates_to_stdout( coordinate_type, label )
        else:
            self.coordinates_to_stdout( coordinate_type )

    def output( self, coordinate_type='Direct', label=None ):
        print( self.title )
        print( self.scaling )
        [ print( ''.join( ['   {: .10f}'.format( element ) for element in row ] ) ) for row in self.lattice ]
        print( ' '.join( self.atoms ) )
        print( ' '.join( [ str(n) for n in self.atom_numbers ] ) )
        print( coordinate_type )
        self.output_coordinates_only( coordinate_type, label )

    def write_to( self, filename, coordinate_type='Direct', label=None ):
        with open( filename, 'w' ) as sys.stdout:
            self.output( coordinate_type=coordinate_type, label=label )

    def output_as_xtl( self ):
        print( self.title )
        print( "CELL" )
        cell_lengths = self.cell_lengths()
        cell_angles  = self.cell_angles()
        cell_data = cell_lengths + cell_angles
        print( ''.join( ['   {: .8f}'.format( element ) for element in cell_data ] ) )
        print( " Symmetry label P1\n\nATOMS\nNAME      X       Y     Z" )
        self.output_coordinates_only( coordinate_type='Direct', label=True )

    def labels( self ):
        return( [ atom_name for ( atom_name, atom_number ) in zip( self.atoms, self.atom_numbers ) for __ in range( atom_number ) ] )

    def replicate( self, h, k, l, group=False ):
        lattice_scaling = np.array( [ h, k, l ], dtype=float )
        lattice_shift = np.reciprocal( lattice_scaling ) 
        new_poscar = Poscar()
        new_poscar.title = self.title
        new_poscar.scaling = self.scaling
        # print( lattice_scaling.dot( self.lattice ) )
        new_poscar.lattice = ( self.lattice.T * lattice_scaling ).T
        new_poscar.coordinate_type = self.coordinate_type
        new_coordinate_list = []
        cell_shift_indices = [ [a,b,c] for a in range(h) for b in range(k) for c in range(l) ]
        # generate grouped / ungrouped atom names and numbers for supercell
        if group:
            new_poscar.atoms = [ label + group for label in self.atoms for group in ('a','b') ]
            new_poscar.atom_numbers = [ int( num * h * k * l / 2 ) for num in self.atom_numbers for __ in ( 0, 1 )]
        else:
            new_poscar.atoms = self.atoms
            new_poscar.atom_numbers = [ num * h * k * l for num in self.atom_numbers ]
        # generate grouped / ungrouped atoms coordinates for supercell
        if group:
            for row in self.coordinates:
                pos_in_origin_cell = row / lattice_scaling
                for odd_even in ( 0, 1 ):
                    for ( a, b, c ) in [ cell_shift for cell_shift in cell_shift_indices if parity(cell_shift) == odd_even ]:
                        new_coordinate_list.append( [ pos_in_origin_cell + np.array( [ a, b, c ] * lattice_shift ) ][0].tolist() )
        else:
            for row in self.coordinates:
                pos_in_origin_cell = row / lattice_scaling
                for ( a, b, c ) in cell_shift_indices:
                    new_coordinate_list.append( [ pos_in_origin_cell + np.array( [ a, b, c ] * lattice_shift ) ][0].tolist() )
        new_poscar.coordinates = np.array( new_coordinate_list )
        return new_poscar    

    def cell_lengths( self ):
        return [ np.linalg.norm( row ) for row in self.lattice * self.scaling ]

    def cell_angles( self ):
        ( a, b, c ) = [ row for row in self.lattice ]
        return [ angle( b, c ), angle( a, c ), angle( a, b ) ]
    
