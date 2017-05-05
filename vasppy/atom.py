class Atom:
    """
    Class for individual atoms
    """

    def __init__( self, label, r ): # currently assume fractional coordinates
        """
        Initialise an Atom instance
    
        Args:
            label (Str): a label for this atom
            r (numpy.array): the atom coordinates
        
        Returns:
            None
        """
        self.label = label
        self.r = r
