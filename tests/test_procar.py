import unittest
import os
from vasppy import procar
import numpy as np

test_data_dir = 'test_data'
test_procar_filename = os.path.join( os.path.dirname( __file__ ), test_data_dir, 'PROCAR_test' )
test_procar_spin_polarised_filename = os.path.join( os.path.dirname( __file__ ), test_data_dir, 'PROCAR_spin_polarised_test' )

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

    def test_spin_polarised_procar_is_read_from_file( self ):
        """Checking that `PROCAR_spin_polarised_test` is read"""
        pcar = procar.Procar()
        pcar.read_from_file( test_procar_spin_polarised_filename )
        self.assertEqual( pcar.spin_channels, 2 )
        self.assertEqual( pcar.number_of_ions, 25 )
        self.assertEqual( pcar.number_of_bands, 112 )
        self.assertEqual( pcar.number_of_k_points, 8 )
         
class ParserTestCase( unittest.TestCase ):
    """Test for VASP output parsers"""

    def test_k_points_are_parsed( self ):
        """Checking that k-points are parsed from PROCAR format strings"""
        procar_string = " k-point    1 :    0.50000000 0.25000000 0.75000000     weight = 0.00806452\nk-point    2 :    0.50000000 0.25735294 0.74264706     weight = 0.00806452"
        self.assertEqual( procar.k_point_parser( procar_string ), 
            [ [ 0.50000000, 0.25000000, 0.75000000 ], 
              [ 0.50000000, 0.25735294, 0.74264706 ] ] )

    def test_negative_k_points_are_parsed( self ): 
        """Checking that negative k-points are parsed from PROCAR format strings"""
        procar_string = "  k-point  119 :   -0.01282051 0.00000000 0.00000000     weight = 0.00500000\n k-point  122 :    0.00000000-0.01282051 0.01282051     weight = 0.00500000\n k-point    1 :   -0.50000000 0.00000000-0.50000000     weight = 0.00500000"
        self.assertEqual( procar.k_point_parser( procar_string ), 
            [ [-0.01282051,  0.00000000,  0.00000000], 
              [ 0.00000000, -0.01282051,  0.01282051],
              [-0.50000000,  0.00000000, -0.50000000]  ] )

    def test_get_numbers_from_string( self ):
        """Checking function for extracting numbers from a string"""
        self.assertEqual( procar.get_numbers_from_string( 'asd834asd2.11 -23as' ), [ 834.0, 2.11, -23.0] )

    def test_projections_are_parsed( self ):
        """Checking that projections are parsed from PROCAR format strings"""
        procar_string = "ion      s      p      d    tot\n  1  0.006  0.000  0.000  0.006\n  2  0.009  0.000  0.000  0.009\ntot  0.835  0.021  0.012  0.868\nion      s      p      d    tot\n  1  0.006  0.000  0.000  0.006\n  2  0.009  0.000  0.000  0.009\ntot  0.835  0.021  0.012  0.868\n"
        self.assertEqual( procar.projections_parser( procar_string ).tolist(),
            np.array( [ [ 1.0, 0.006, 0.0, 0.0, 0.006, 2.0, 0.009, 0.0, 0.0, 0.009, 0.0, 0.835, 0.021, 0.012, 0.868 ], [ 1.0, 0.006, 0.0, 0.0, 0.006, 2.0, 0.009, 0.0, 0.0, 0.009, 0.0, 0.835, 0.021, 0.012, 0.868 ] ] ).tolist() )

class ProcarSupportFunctionsTestCase( unittest.TestCase ):
    """Test for the support functions in procar.py"""

    def test_area_of_a_triangle_in_cartesian_space( self ):
        a = np.array( [ 0.0, 0.0, 0.0 ] )
        b = np.array( [ 4.0, 0.0, 0.0 ] )
        c = np.array( [ 0.0, 3.0, 0.0 ] )
        self.assertEqual( procar.area_of_a_triangle_in_cartesian_space( a, b, c ), 6.0 )

if __name__ == '__main__':
    unittest.main()
