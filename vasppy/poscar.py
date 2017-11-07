import numpy as np
import sys
import re
import copy
from vasppy import configuration, atom, cell
from pymatgen import Lattice as pmg_Lattice
from pymatgen import Structure as pmg_Structure
from pymatgen.io.cif import CifWriter
from collections import Counter

# Ignore SIG_PIPE and don't throw exceptions on it... 
# http://newbebweb.blogspot.co.uk/2012/02/python-head-ioerror-errno-32-broken.html
from signal import signal, SIGPIPE, SIG_DFL
signal( SIGPIPE, SIG_DFL ) 

def parity( list ):
    return( sum( list )%2 )

def swap_axes( matrix, axes ):
    axes_index = { 'x': 0, 'y': 1, 'z': 2 }
    matrix[:, [ axes_index[ axes[ 0 ] ], axes_index[ axes[ 1 ] ] ] ] = matrix[:, [ axes_index[ axes[ 1 ] ], axes_index[ axes[ 0 ] ] ] ]
    return matrix

class Poscar:

    lines_offset = 9

    def __init__( self ):
        self.title = "Title"
        self.scaling = 1.0
        self.cell = cell.Cell( np.identity( 3 ) )
        self.atoms = [ 'A' ]
        self.atom_numbers = [ 1 ]
        self.coordinate_type = 'Direct'
        self.coordinates = np.array( [ [ 0.0, 0.0, 0.0 ] ] )
        self.selective_dynamics = False

    def coordinates_by_species( self, species ):
        return self.coordinates[ self.range_by_species( species ) ]
 
    def range_by_species( self, species ):
        i = 0
        atom_range = {}
        for a, n in zip( self.atoms, self.atom_numbers ):
            atom_range[ a ] = range(i, i+n)
            i += n
        return atom_range[ species ]

    def atom_number_by_species( self, species ):
        atom_numbers = { a : n for a, n in zip( self.atoms, self.atom_numbers ) }
        return atom_numbers[ species ]

    def sorted( self, species ):
        new_poscar = copy.deepcopy( self )
        new_poscar.atoms = species
        new_poscar.atom_numbers = [ self.atom_number_by_species( s ) for s in species ]
        new_poscar.coordinates = np.concatenate( [ self.coordinates_by_species( s ) for s in species ], axis = 0 )
        return new_poscar 
 
    def read_from( self, filename ):
        try:
            with open( filename ) as f:
                lines = f.readlines()
        except FileNotFoundError:
            print( "\"" + filename + "\" not found", file=sys.stderr ) 
            sys.exit( -2 )
        self.title = lines.pop(0).strip()
        self.scaling = float( lines.pop(0).strip() ) 
        self.cell.matrix = np.array( [ [ float( e ) for e in lines.pop(0).split() ] for i in range( 3 ) ] )
        self.cell.inv_matrix = np.linalg.inv( self.cell.matrix )
        self.atoms = lines.pop(0).split()
        self.atom_numbers = [ int(element) for element in lines.pop(0).split() ]
        self.coordinate_type = lines.pop(0).strip()
        if re.match( r'\A[Ss]', self.coordinate_type ): # test for 'Selective dynamics'
            self.selective_dynamics = True
            self.coordinate_type = lines.pop( 0 ).strip()
        self.coordinates = np.array( [ [ float( e ) for e in lines.pop(0).split()[0:3] ] for i in range( sum( self.atom_numbers ) ) ] )
        if self.coords_are_cartesian(): # Convert to direct coordinates
            self.coordinates = self.fractional_coordinates()
            self.coordinate_type = 'Direct'

    @classmethod
    def from_file( cls, filename ):
        poscar = cls()
        poscar.read_from( filename )
        return poscar

    def in_bohr( self ):
        new_poscar = copy.deepcopy( self )
        bohr_to_angstrom = 0.529177211
        new_poscar.scaling *= bohr_to_angstrom
        new_poscar.cell.matrix /= bohr_to_angstrom
        if new_poscar.coords_are_cartesian():
            new_poscar.coordinates /= bohr_to_angstrom
        return( new_poscar )

    def coords_are_fractional( self ): 
        return not self.coords_are_cartesian()

    def coords_are_cartesian( self ):
        return re.match( r'\A[CcKk]', self.coordinate_type )

    def fractional_coordinates( self ):
        return ( self.coordinates if self.coords_are_fractional() else self.coordinates.dot( np.linalg.inv( self.cell.matrix ) ) )
        # return ( self.coordinates if self.coords_are_fractional() else self.cell.cartesian_to_fractional_coordinates( coordinates )
            # need to check whether this alternative version works.

    def cartesian_coordinates( self ):
        return ( self.coordinates if self.coords_are_cartesian() else self.coordinates.dot( self.cell.matrix ) )

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

    def numbered_coordinates_to_stdout( self, coordinate_type='Direct' ):
        coord_opts = { 'Direct'    : self.fractional_coordinates(), 
                       'Cartesian' : self.cartesian_coordinates() } 
        try:
            for i, row in enumerate( coord_opts[ coordinate_type ] ):
                print( ''.join( ['  {: .10f}'.format( element ) for element in row ] ) + '  ' + str(i+1) )
        except KeyError: 
            raise Exception( 'Passed coordinate_type: ' + coordinate_type + '\nAccepted values: [ Direct | Cartesian ] ' )

    def output_coordinates_only( self, coordinate_type='Direct', opts=None ):
        if opts is None:
            opts = {}
        if opts.get( 'numbered' ):
            self.numbered_coordinates_to_stdout( coordinate_type )
        elif opts.get( 'label' ):
            self.labelled_coordinates_to_stdout( coordinate_type, opts[ 'label' ] )
        else:
            self.coordinates_to_stdout( coordinate_type )

    def output( self, coordinate_type='Direct', opts=None ):
        if opts is None:
            opts = {}
        if not opts.get( 'coordinates_only' ):
            print( self.title )
            print( self.scaling )
            [ print( ''.join( ['   {: .10f}'.format( element ) for element in row ] ) ) for row in self.cell.matrix ]
            print( ' '.join( self.atoms ) )
            print( ' '.join( [ str(n) for n in self.atom_numbers ] ) )
            print( coordinate_type )
        self.output_coordinates_only( coordinate_type=coordinate_type, opts=opts )

    def write_to( self, filename, coordinate_type='Direct', opts=None ):
        if opts is None:
            opts = {}
        with open( filename, 'w' ) as sys.stdout:
            self.output( coordinate_type=coordinate_type, opts=opts )
        sys.stdout = sys.__stdout__ # make sure sys.stdout is reset
    def output_as_xtl( self ):
        print( self.title )
        print( "CELL" )
        cell_lengths = self.cell_lengths()
        cell_angles  = self.cell_angles()
        cell_data = cell_lengths + cell_angles
        print( ''.join( ['   {: .8f}'.format( element ) for element in cell_data ] ) )
        print( " Symmetry label P1\n\nATOMS\nNAME      X       Y     Z" )
        output_opts = { 'label' : True }
        self.output_coordinates_only( coordinate_type='Direct', opts = output_opts )

    def output_as_cif( self, symprec = None ):
        print( CifWriter( self.to_pymatgen_structure(), symprec ) )

    def output_as_pimaim( self, to_bohr = True ):
        if to_bohr is True:
            unit_scaling = 0.52918
        else:
            unit_scaling = 1.0
        cell_lengths = self.cell.lengths() * self.scaling / unit_scaling
        print( "T\nF\nF\nF" )
        for row in self.cell_coordinates():
            print( ' '.join( [ str( x / unit_scaling ) for x in row ] ) )
        for row in self.cell.unit_vectors().transpose():
            print( ' '.join( [ str( x ) for x in row ] ) )
        for length in cell_lengths:
            print( length )

    def cell_coordinates( self ):
        return( self.coordinates * self.cell.lengths() * self.scaling )

    def labels( self ):
        return( [ atom_name for ( atom_name, atom_number ) in zip( self.atoms, self.atom_numbers ) for __ in range( atom_number ) ] )

    def replicate( self, h, k, l, group=False ):
        lattice_scaling = np.array( [ h, k, l ], dtype=float )
        lattice_shift = np.reciprocal( lattice_scaling ) 
        new_poscar = Poscar()
        new_poscar.title = self.title
        new_poscar.scaling = self.scaling
        # print( lattice_scaling.dot( self.lattice ) )
        new_poscar.cell.matrix = ( self.cell.matrix.T * lattice_scaling ).T
        new_poscar.coordinate_type = self.coordinate_type
        new_coordinate_list = []
        cell_shift_indices = [ [ a,b,c ] for a in range(h) for b in range(k) for c in range(l) ]
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
        return( ( self.cell.lengths() * self.scaling ).tolist() )
        # return [ np.linalg.norm( row ) for row in self.cell.matrix * self.scaling ]

    def cell_angles( self ):
        return( self.cell.angles() )
        # ( a, b, c ) = [ row for row in self.cell.matrix ]
        # return [ angle( b, c ), angle( a, c ), angle( a, b ) ]

    def to_configuration( self ):
        atoms = [ atom.Atom( label, coordinates ) for ( label, coordinates ) in zip( self.labels(), self.fractional_coordinates() ) ]
        config = configuration.Configuration( cell.Cell( matrix = self.cell.matrix * self.scaling ), atoms )
        return( config )

    def swap_axes( self, axes ):
        new_poscar = copy.deepcopy( self )
        new_poscar.cell = cell.Cell( swap_axes( self.cell.matrix, axes ) )
        return new_poscar

    def to_pymatgen_structure( self ):
        lattice = pmg_Lattice( self.cell.matrix * self.scaling )
        structure = pmg_Structure( lattice, self.labels(), self.coordinates )
        return structure 

    @property
    def stoichiometry( self ):
        """
        Stoichiometry for this POSCAR, as a Counter.
        e.g. AB_2O_4 -> Counter( { 'A': 1, 'B': 2, O: 4 } )
        
        Args:
            None

        Returns:
            None
        """
        return Counter( { label: number for label, number in zip( self.atoms, self.atom_numbers ) } )
