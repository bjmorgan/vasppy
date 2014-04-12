import numpy as np
import sys
from vasppy import poscar

class Grid:

    projections = { 'x' : 0, 'y' : 1, 'z' : 2 }

    def read_from_filename( self, filename ):
        self.filename = filename
        self.poscar = poscar.Poscar()
        self.poscar.read_from( self.filename )
        self.number_of_header_lines = sum( self.poscar.atom_numbers ) + poscar.Poscar.lines_offset
        self.read_dimensions()
        self.read_grid()

    def write_to_filename( self, filename ):
        with open( filename, 'w' ) as file_out:
            sys.stdout = file_out
            self.poscar.output()
            self.write_dimensions()
            self.write_grid()

    def read_dimensions( self ):
        with open( self.filename, 'r' ) as file_in:
            for i, line in enumerate( file_in ):
                if i == self.number_of_header_lines:
                    self.dimensions = [ int(i) for i in line.split() ]
                    break

    def write_dimensions( self ):
        print( "\n" + ' '.join( [ str(i) for i in self.dimensions ] ) ) 

    def read_grid( self ):
        self.grid = np.swapaxes( np.loadtxt( self.filename, skiprows = self.number_of_header_lines + 1 ).reshape( self.dimensions[::-1] ), 0, 2 )

    def write_grid( self ):
        for row in np.reshape( np.swapaxes( self.grid, 0, 2 ), ( -1, 5 ) ):
            print( ''.join( [ ' {:.11E}'.format( element ) for element in row ] ) )

    def average( self, normal_axis_label ):
        axes = [ 0, 1, 2 ]
        axes.remove( Grid.projections[ normal_axis_label ] )
        return( np.sum( np.sum( self.grid, axis=axes[1] ), axis=axes[0] ) / ( self.dimensions[0] * self.dimensions[1] ) )
