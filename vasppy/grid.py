import numpy as np
import sys
import math
from vasppy import poscar, cell

def interpolate( i, j, x ):
    return( ( i * ( 1.0 - x ) ) + ( j * x) )

def trilinear_interpolation( cube, r ):
    return( interpolate ( 
                interpolate( 
                    interpolate( cube[ 0, 0, 0 ], cube[ 1, 0, 0 ], r[ 0 ] ), # trilinear interpolation => http://en.wikipedia.org/wiki/Trilinear_interpolation
                    interpolate( cube[ 0, 1, 0 ], cube[ 1, 1, 0 ], r[ 0 ] ), 
                    r[ 1 ] ),
                interpolate( 
                    interpolate( cube[ 0, 0, 1 ], cube[ 1, 0, 1 ], r[ 0 ] ),
                    interpolate( cube[ 0, 1, 1 ], cube[ 1, 1, 1 ], r[ 0 ] ), 
                    r[ 1 ] ), 
                r[ 2 ] ) 
            )

class Grid:

    projections = { 'x' : 0, 'y' : 1, 'z' : 2 }

    def __init__( self, dimensions = [ 1, 1, 1 ] ):
        self.filename = None
        self.poscar = poscar.Poscar()
        self.number_of_header_lines = 0
        self.dimensions = dimensions
        self.spacing = np.array( [ 1.0 / number_of_points for number_of_points in self.dimensions ] )
        self.grid = np.zeros( self.dimensions )

    def read_from_filename( self, filename ):
        self.filename = filename
        self.poscar = poscar.Poscar()
        self.poscar.read_from( self.filename )
        self.number_of_header_lines = sum( self.poscar.atom_numbers ) + poscar.Poscar.lines_offset
        self.read_dimensions()
        self.read_grid()
        return self

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
        grid_data_lines = math.ceil( ( self.dimensions[0] * self.dimensions[1] * self.dimensions[2] ) / 5 )
        with open( self.filename ) as file_in:
            for i, line in enumerate( file_in ):
                if ( i > self.number_of_header_lines ) and ( i <= self.number_of_header_lines + grid_data_lines ):
                    grid_data.append( line.strip() )
        grid_data = np.array( [ float( s ) for s in ' '.join( grid_data ).split() ] )
        self.grid = np.reshape( grid_data, tuple( self.dimensions ), order = 'F' )

    def write_grid( self ):
        np.savetxt( sys.stdout.buffer, np.swapaxes( self.grid, 0, 2 ).reshape( -1, 5 ), fmt='%.11E' )

    def average( self, normal_axis_label ):
        axes = [ 0, 1, 2 ]
        axes.remove( Grid.projections[ normal_axis_label ] )
        return( np.sum( np.sum( self.grid, axis=axes[1] ), axis=axes[0] ) / ( self.dimensions[0] * self.dimensions[1] ) )

    def by_index( self, index ):
        return self.grid[ index[0], index[1], index[2] ]

    def fractional_coordinate_at_index( self, index ):      
        return( np.multiply( self.spacing, index ) )

    def cartesian_coordinate_at_index( self, index ):
        return( self.fractional_coordinate_at_index( index ).dot( self.poscar.cell.matrix ) )

    def cube_slice( self, x0, y0, z0 ):
        x1 = ( x0 + 1 ) % self.dimensions[0]
        y1 = ( y0 + 1 ) % self.dimensions[1]
        z1 = ( z0 + 1 ) % self.dimensions[2]
        cube = np.array( [ self.grid[ x0, y0, z0 ],
                           self.grid[ x0, y0, z1 ],
                           self.grid[ x0, y1, z0 ],
                           self.grid[ x0, y1, z1 ],
                           self.grid[ x1, y0, z0 ],
                           self.grid[ x1, y0, z1 ],
                           self.grid[ x1, y1, z0 ],
                           self.grid[ x1, y1, z1 ] ] ).reshape( 2, 2, 2 )
        return( cube )

    def interpolated_value_at_fractional_coordinate( self, coord ):
        point = np.multiply( np.array( self.dimensions ), coord ) # point contains the (fractional) index of the coordinate coord.
        origin = [ int( f ) for f in point ]                      # origin contains the 3D index of the lowest-index point in the cube surrounding point (i,j,k)
        delta = [ p - o for p, o in zip( point, origin ) ]        # delta contains the *fractional* offset of "point" from "origin"
        cube = self.cube_slice( *origin )                         # cube contains the data values at the 8 bounding grid points
        return( trilinear_interpolation( cube, delta ) )

    def interpolate_to_orthorhombic_grid( self, dimensions ): # warning. This may need a more robust minimim image function in Cell.py for highly non-orthorhombic cells
        old_grid = self
        old_cell = old_grid.poscar.cell
        new_grid = Grid( dimensions = dimensions )
        new_grid.poscar.cell = cell.Cell( np.diag( np.diag( self.poscar.cell.matrix ) ) )
        index_grid           = np.array( [ [ i, j, k ] for ( i, j, k ), value in np.ndenumerate( new_grid.grid ) ] )
        cart_coord_grid      = np.array( [ new_grid.cartesian_coordinate_at_index( index ) for index in index_grid ] )
        init_frac_coord_grid = np.array( [ old_cell.cartesian_to_fractional_coordinates( r ) for r in cart_coord_grid ] )
        frac_coord_grid      = np.array( [ old_cell.inside_cell( r ) for r in init_frac_coord_grid ] )
        new_grid_data        = np.array( [ old_grid.interpolated_value_at_fractional_coordinate( r ) for r in frac_coord_grid ] )
        new_grid.grid = new_grid_data.reshape( new_grid.dimensions )
        return( new_grid )
