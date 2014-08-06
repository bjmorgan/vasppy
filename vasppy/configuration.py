from vasppy import atom, cell
import numpy as np

class Configuration:

    def __init__( self, cell, atoms ):
        self.cell  = cell 
        self.atoms = atoms

    def dr( self, atom1, atom2 ):
        return( self.cell.dr( atom1.r, atom2.r ) )

    def minimum_image_dr( self, atom1, atom2, cutoff = None ):
        return( self.cell.minimum_image_dr( atom1.r, atom2.r, cutoff = cutoff ) )

    def interatomic_distances( self, minimum_image_convention = True ):
        return( np.array( [ [ self.minimum_image_dr( atom_i, atom_j ) for atom_j in self.atoms ] for atom_i in self.atoms ] ) )

    def interatomic_distances_for_atom( self, atom1, minimum_image_convention = True ):
        return( np.array( [ self.minimum_image_dr( atom1, atom2 ) for atom2 in self.atoms ] ) )

        


