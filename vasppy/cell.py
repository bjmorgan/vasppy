import math
import numpy as np

def angle( x, y ):
    """
    Calculate the angle between two vectors, in degrees.

    Args:
        x (np.array): one vector.
        y (np.array): the other vector.

    Returns:
        (float):      the angle between x and y in degrees.
    """
    dot = np.dot( x, y )
    x_mod = np.linalg.norm( x )
    y_mod = np.linalg.norm( y )
    cos_angle = dot / ( x_mod * y_mod )
    return np.degrees( np.arccos( cos_angle ) )

def rotation_matrix(axis, theta):
    """
    Return the 3D rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    
    Args:
        axis (np.array): length 3 numpy array defining the axis of rotation.
        theta (float):   rotation angle in radians.

    Returns:
        (np.array):      the corredponding rotation matrix.
    """
    axis = np.asarray( axis )
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
        """
        Initialise a Cell object.

        Args:
            matrix (np.array): 3x3 numpy array containing the cell matrix.

        Returns:
            None
        """
        assert type( matrix ) is np.ndarray
        assert matrix.shape == ( 3, 3 )
        self.matrix = matrix # 3 x 3 numpy Array
        self.inv_matrix = np.linalg.inv( matrix )

    def dr( self, r1, r2, cutoff=None ):
        """
        Calculate the distance between two fractional coordinates in the cell.
        
        Args:
            r1 (np.array): fractional coordinates for position 1.
            r2 (np.array): fractional coordinates for position 2.
            cutoff (optional:Bool): If set, returns None for distances greater than the cutoff. Default None (unset).

        Returns:
            (float): the distance between r1 and r2.
        """
        delta_r_cartesian = ( r1 - r2 ).dot( self.matrix )
        delta_r_squared = sum( delta_r_cartesian**2 )
        if cutoff != None:
            cutoff_squared = cutoff ** 2
            if delta_r_squared > cutoff_squared:
                return None
        return( math.sqrt( delta_r_squared ) )

    def nearest_image( self, origin, point ):
        """
        Find the fractional_coordinates of the nearest periodic image to a point of origin.

        Args:
            origin (np.array): fractional coordinates of the point of origin.
            point  (np.array): fractional coordinates of the other point.

        Returns:
            (np.array): the fractional coordinates of the nearest image of `point` to `origin`.
        """
        return( origin + self.minimum_image( origin, point ) )

    def minimum_image( self, r1, r2 ):
        """
        Find the minimum image vector from point r1 to point r2.

        Args:
            r1 (np.array): fractional coordinates of point r1.
            r2 (np.array): fractional coordinates of point r2.

        Returns:
            (np.array): the fractional coordinate vector from r1 to the nearest image of r2.
        """
        delta_r = r2 - r1
        delta_r = np.array( [ x - math.copysign( 1.0, x ) if abs(x) > 0.5 else x for x in delta_r ] )
        return( delta_r )

    def minimum_image_dr( self, r1, r2, cutoff=None ):
        """
        Calculate the shortest distance between two points in the cell, 
        accounting for periodic boundary conditions.

        Args:
            r1 (np.array): fractional coordinates of point r1.
            r2 (np.array): fractional coordinates of point r2.
            cutoff (:obj: `float`, optional): if set, return zero if the minimum distance is greater than `cutoff`. Defaults to None.

        Returns:
            (float): The distance between r1 and r2.
        """
        delta_r_vector = self.minimum_image( r1, r2 )
        return( self.dr( np.zeros( 3 ), delta_r_vector, cutoff ) )

    def lengths( self ):
        """
        The cell lengths.

        Args:
            None

        Returns:
            (np.array(a,b,c)): The cell lengths.
        """
        return( np.array( [ math.sqrt( sum( row**2 ) ) for row in self.matrix ] ) )

    def angles( self ):
        """
        The cell angles (in degrees).

        Args:
            None

        Returns:
            (list(alpha,beta,gamma)): The cell angles.
        """
        ( a, b, c ) = [ row for row in self.matrix ]
        return [ angle( b, c ), angle( a, c ), angle( a, b ) ]

    def cartesian_to_fractional_coordinates( self, coordinates ):
        """
        Convert a set of Cartesian coordinates to fractional coordinates in the cell.

        Args:
            coordinates (np.array(dim(N,3))): The set of Cartesian coordinates.

        Returns:
            (np.array(dim(N,3))): The corresponding set of fractional coordinates.
        """
        return( coordinates.dot( self.inv_matrix ) )

    def fractional_to_cartesian_coordinates( self, coordinates ):
        """
        Convert a set of fractional coordinates in the cell to Cartesian coordinates.

        Args:
            coordinates (np.array(dim(N,3))): The set of fractional coordinates.

        Returns:
            (np.array(dim(N,3))): The corresponding set of Cartesian coordinates.
        """
        return( coordinates.dot( self.matrix ) )

    def inside_cell( self, r ):
        """
        Given a fractional-coordinate, if this lies outside the cell return the equivalent point inside the cell.

        Args:
            r (np.array): Fractional coordinates of a point (this may be outside the cell boundaries).

        Returns:
            (np.array): Fractional coordinates of an equivalent point, inside the cell boundaries.
        """
        centre = np.array( [ 0.5, 0.5, 0.5 ] )
        new_r = self.nearest_image( centre, r )
        return new_r

    def volume( self ):
        """
        The cell volume.

        Args:
            None

        Returns:
            (float): The cell volume.
        """
        return np.dot( self.matrix[0], np.cross( self.matrix[1], self.matrix[2] ) ) 

    def unit_vectors( self ):
        """
        The unit vectors for the cell vectors.

        Args:
            None

        Returns:
            (np.array): The unit vectors for the cell vectors.
        """
        return( ( self.matrix.transpose() / self.lengths() ).transpose() )

    def rotate( self, axis, theta ):
        self.matrix = np.array( [ np.dot( rotation_matrix(axis,theta), v) for v in self.matrix ] )
