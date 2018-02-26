from vasppy import atom, cell, rdf
import numpy as np

class Configuration:
    """
    A Configuration object stores a single structure.
    """

    def __init__( self, cell, atoms ):
        self.cell  = cell 
        self.atoms = atoms

    def dr( self, atom1, atom2 ):
        """
        Calculate the distance between two atoms.

        Args:
            atom1 (vasppy.Atom): Atom 1.
            atom2 (vasppy.Atom): Atom 2.

        Returns:
            (float): The distance between Atom 1 and Atom 2.
        """
        return self.cell.dr( atom1.r, atom2.r )

    def minimum_image_dr( self, atom1, atom2, cutoff=None ):
        return self.cell.minimum_image_dr( atom1.r, atom2.r, cutoff=cutoff )

    def interatomic_distances( self, minimum_image_convention = True ):
        return np.array( [ [ self.minimum_image_dr( atom_i, atom_j ) for atom_j in self.atoms ] for atom_i in self.atoms ] )

    def interatomic_distances_for_atom( self, atom1, minimum_image_convention = True ):
        return np.array( [ self.minimum_image_dr( atom1, atom2 ) for atom2 in self.atoms ] )

    def atoms_with_label( self, label ):
        return filter( lambda atom: atom.label == label, self.atoms )        

    def partial_rdf( self, spec_i, spec_j, max_r, number_of_bins ):
        this_rdf = rdf.Rdf( max_r, number_of_bins )
        atoms_i = list( self.atoms_with_label( spec_i ) )
        atoms_j = list( self.atoms_with_label( spec_j ) )
        for atom_i in atoms_i:
            for atom_j in atoms_j:
                if atom_i is atom_j:
                    continue
                dr = self.minimum_image_dr( atom_i, atom_j )
                if dr <= max_r:
                    this_rdf.add_dr( dr )
        return this_rdf 

    def per_atom_rdf( self, spec_i, spec_j, max_r, number_of_bins ):
        rdfs = []
        atoms_i = list( self.atoms_with_label( spec_i ) )
        atoms_j = list( self.atoms_with_label( spec_j ) )
        for atom_i in atoms_i:
            this_rdf = rdf.Rdf( max_r, number_of_bins )
            for atom_j in atoms_j:
                if atom_i is atom_j:
                    continue
                dr = self.minimum_image_dr( atom_i, atom_j )
                try:
                    this_rdf.add_dr( dr )
                except IndexError:
                    pass
                except:
                    raise
            rdfs.append( this_rdf )
        return rdfs

