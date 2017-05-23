import unittest
from unittest.mock import Mock, patch
import numpy as np

from vasppy.configuration import Configuration
from vasppy.cell import Cell
from vasppy.atom import Atom

class TestConfiguration( unittest.TestCase ):

    def setUp( self ):
        self.cell = Mock( spec=Cell )
        self.atoms = [ Mock( spec=Atom ), Mock( spec=Atom ) ]
        self.configuration = Configuration( self.cell, self.atoms )

    def test_configuration_is_initialised( self ):
        self.assertEqual( self.configuration.cell, self.cell )
        self.assertEqual( self.configuration.atoms, self.atoms )

    def test_dr( self ):
        self.atoms[0].r = np.array( [ 0.0, 0.0, 0.0 ] )
        self.atoms[1].r = np.array( [ 0.1, 0.1, 0.1 ] )
        self.cell.dr = Mock( return_value=0.5477225575 )
        self.assertEqual( self.configuration.dr( self.atoms[0], self.atoms[1] ), 0.5477225575 )

if __name__ == '__main__':
    unittest.main()
