import unittest
from unittest.mock import Mock, patch
from io import StringIO

from vasppy.poscar import Poscar
from vasppy.cell import Cell
import numpy as np
from collections import Counter

class PoscarTestCase( unittest.TestCase ):

    def setUp( self ):
        self.poscar = Poscar()
        self.poscar.title = "Title"
        self.poscar.scaling = 1.0
        self.poscar.cell = Mock( spec=Cell )
        self.poscar.cell.matrix = np.identity( 3 )
        self.poscar.atoms = [ 'A' ]
        self.poscar.atom_numbers = [ 1 ]
        self.poscar.coordinate_type = 'Direct'
        self.poscar.coordinates = np.array( [ [ 0.0, 0.0, 0.0 ] ] )
        self.poscar.selective_dynamics = False

    def test_stoichiometry( self ):
        poscar = Poscar()
        poscar.atoms = [ 'A', 'B', 'C' ]
        poscar.atom_numbers = [ 1, 2, 3 ]
        self.assertEqual( poscar.stoichiometry, Counter( { 'A': 1, 'B': 2, 'C': 3 } ) )

    @patch('sys.stdout', new_callable=StringIO)
    def test_output_header( self, mock_stdout ):
        self.poscar.output_header()
        expected_header_string = ("Title\n"
                                  "1.0\n"
                                  "    1.0000000000    0.0000000000    0.0000000000\n"
                                  "    0.0000000000    1.0000000000    0.0000000000\n"
                                  "    0.0000000000    0.0000000000    1.0000000000\n"
                                  "A\n"
                                  "1\n"
                                  "Direct\n")
        self.assertEqual( mock_stdout.getvalue(), expected_header_string )

if __name__ == '__main__':
    unittest.main()
