import numpy as np
import copy

class Polyhedron:

    def __init__( self, vertices, cell, inside_point, cutoff ):
        self.vertices = vertices
        self.cell = cell
        self.inside_point = inside_point
        self.vertices = [ self.cell.nearest_image( inside_point, p ) for p in vertices ]
        self.inside_point = self.centre() # improved inside_point

    def centre( self ):
        return( sum( self.vertices ) / len( self.vertices ) )

    def print_points( self ):
        for point in self.vertices:
            print( point.dot( self.cell.matrix ) )