import numpy as np

class Atom:
    """
    Class for individual atoms
    """

    def __init__(self,
                 label: str,
                 r: np.ndarray) -> None: # currently assume fractional coordinates
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
