import unittest
import os
from vasppy import procar

test_procar_filename = os.path.join( os.path.dirname( __file__ ), 'PROCAR_test' )

class ProcarTestCase( unittest.TestCase ):
    """Test for `procar.py`"""

    def setUp( self ):
        self.procar = procar.Procar()

    def test_procar_is_read_from_file( self ):
        """Checking that `PROCAR_test` is read"""
        pcar = procar.Procar()
        pcar.read_from_file( test_procar_filename )
        self.assertEqual( pcar.spin_channels, 4 )
        self.assertEqual( pcar.number_of_ions, 22 )
        self.assertEqual( pcar.number_of_bands, 4 )
        self.assertEqual( pcar.number_of_k_points, 2 )

class ParserTestCase( unittest.TestCase ):
    """Test for VASP output parsers"""

    def test_k_points_are_parsed( self ):
        """Checking that k-points are parsed from PROCAR format strings"""
        procar_string = " k-point    1 :    0.50000000 0.25000000 0.75000000     weight = 0.00806452\nk-point    2 :    0.50000000 0.25735294 0.74264706     weight = 0.00806452"
        self.assertEqual( procar.k_point_parser( procar_string ), 
            [ ['0.50000000', '0.25000000', '0.75000000'], 
              ['0.50000000', '0.25735294', '0.74264706'] ] )

    def test_negative_k_points_are_parsed( self ): 
        """Checking that negative k-points are parsed from PROCAR format strings"""
        procar_string = "  k-point  119 :   -0.01282051 0.00000000 0.00000000     weight = 0.00500000\n k-point  122 :    0.00000000-0.01282051 0.01282051     weight = 0.00500000\n k-point    1 :   -0.50000000 0.00000000-0.50000000     weight = 0.00500000"
        self.assertEqual( procar.k_point_parser( procar_string ), 
            [ ['-0.01282051',  '0.00000000',  '0.00000000'], 
              [ '0.00000000', '-0.01282051',  '0.01282051'],
              ['-0.50000000',  '0.00000000', '-0.50000000']  ] )

if __name__ == '__main__':
    unittest.main()
