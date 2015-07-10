import math
import numpy as np

def angle( x, y ):
    dot = np.dot( x, y )
    x_mod = np.linalg.norm( x )
    y_mod = np.linalg.norm( y )
    cos_angle = dot / ( x_mod * y_mod )
    return np.degrees( np.arccos( cos_angle ) )

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray (axis )
    theta = np.asarray( theta )
    axis = axis / math.sqrt( np.dot( axis, axis ) )
    a = math.cos( theta / 2 )
    b, c, d = -axis * math.sin( theta / 2 )
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array( [ [ aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac) ],
                       [ 2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab) ],
                       [ 2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc ] ] )

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

    def rotate( self, axis, theta ):
        self.matrix = np.array( [ np.dot( rotation_matrix(axis,theta), v) for v in self.matrix ] )
