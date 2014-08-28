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
            sys.stdout.flush()
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
        grid_data = []
        grid_data_lines = ( self.dimensions[0] * self.dimensions[1] * self.dimensions[2] ) // 5
        with open( self.filename ) as file_in:
            for i, line in enumerate( file_in ):
                if ( i > self.number_of_header_lines ) and ( i <= self.number_of_header_lines + grid_data_lines + 1):
                    grid_data.append( line.strip() )
        grid_data = np.array( [ float( s ) for s in ' '.join( grid_data ).split() ] )
        print( grid_data.shape )
        self.grid = np.reshape( grid_data, tuple( self.dimensions ), order = 'F' )

    def write_grid( self ):
        np.savetxt( sys.stdout.buffer, np.swapaxes( self.grid, 0, 2 ).reshape( -1, 5 ), fmt='%.11E' )

    def average( self, normal_axis_label ):
        axes = [ 0, 1, 2 ]
        axes.remove( Grid.projections[ normal_axis_label ] )
        return( np.sum( np.sum( self.grid, axis=axes[1] ), axis=axes[0] ) / ( self.dimensions[0] * self.dimensions[1] ) )
