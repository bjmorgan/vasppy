import math
import numpy as np

def angle( x, y ):
    dot = np.dot( x, y )
    x_mod = np.linalg.norm( x )
    y_mod = np.linalg.norm( y )
    cos_angle = dot / ( x_mod * y_mod )
    return np.degrees( np.arccos( cos_angle ) )

class Cell:

    def __init__( self, matrix ):
        assert type( matrix ) is np.ndarray
        assert matrix.shape == ( 3, 3 )
        self.matrix = matrix # 3 x 3 numpy Array
        self.inv_matrix = np.linalg.inv( matrix )

    def dr( self, r1, r2, cutoff = None ):
        delta_r_cartesian = ( r1 - r2 ).dot( self.matrix )
        delta_r_squared = sum( delta_r_cartesian**2 )
        if cutoff != None:
            cutoff_squared = cutoff ** 2
            if delta_r_squared > cutoff_squared:
                return None
        return( math.sqrt( delta_r_squared ) )

    def nearest_image( self, origin, point ):
        return( origin + self.minimum_image( origin, point ) )

    def minimum_image( self, r1, r2 ):
        delta_r = r2 - r1
        # print( delta_r )
        delta_r = np.array( [ x - math.copysign( 1.0, x ) if abs(x) > 0.5 else x for x in delta_r ] )
        # print( delta_r )
        return( delta_r )

    def minimum_image_dr( self, r1, r2, cutoff = None ):
        delta_r_vector = self.minimum_image( r1, r2 )
        return( self.dr( np.zeros( 3 ), delta_r_vector, cutoff ) )

    def lengths( self ):
        return( np.array( [ math.sqrt( sum( row**2 ) ) for row in self.matrix ] ) )

    def angles( self ):
        ( a, b, c ) = [ row for row in self.matrix ]
        return [ angle( b, c ), angle( a, c ), angle( a, b ) ]

    def centre( self, points ):
        print( points )

    def cartesian_to_fractional_coordinates( self, coordinates ):
        return( coordinates.dot( self.inv_matrix ) )

    def fractional_to_cartesian_coordinates( self, coordinates ):
        return( coordinates.dot( self.matrix ) )

    def inside_cell( self, r ):
        centre = np.array( [ 0.5, 0.5, 0.5 ] )
        # print( 'r=',r )
        new_r = self.nearest_image( centre, r )
        # print( 'new_r=',new_r )
        return new_r

    def volume( self ):
        return np.dot( self.matrix[0], np.cross( self.matrix[1], self.matrix[2] ) ) 

    def unit_vectors( self ):
        return( ( self.matrix.transpose() / self.lengths() ).transpose() )
