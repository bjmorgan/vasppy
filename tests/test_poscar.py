import unittest
from unittest.mock import Mock, patch, mock_open

import vasppy.poscar 
import numpy as np
from collections import Counter

class PoscarTestCase( unittest.TestCase ):

    def test_stoichiometry( self ):
        poscar = vasppy.poscar.Poscar()
        poscar.atoms = [ 'A', 'B', 'C' ]
        poscar.atom_numbers = [ 1, 2, 3 ]
        self.assertEqual( poscar.stoichiometry, Counter( { 'A': 1, 'B': 2, 'C': 3 } ) )

if __name__ == '__main__':
    unittest.main()
