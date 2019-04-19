import unittest
import os
from vasppy import procar
import numpy as np
from unittest.mock import patch
import warnings

test_data_dir = 'test_data'
test_procar_filename = os.path.join( os.path.dirname( __file__ ), test_data_dir, 'PROCAR_test' )
test_procar_spin_polarised_filename = os.path.join( os.path.dirname( __file__ ), test_data_dir, 'PROCAR_spin_polarised_test' )

class ProcarTestCase( unittest.TestCase ):
    """Test for procar.py"""

    def setUp( self ):
        self.procar = procar.Procar()

    def test_procar_is_initialised( self ):
        pcar = procar.Procar()
        self.assertEqual( pcar.spin_channels, 1 )
        self.assertEqual( pcar.number_of_k_points, None )
        self.assertEqual( pcar.number_of_ions, None )
        self.assertEqual( pcar.number_of_bands, None )
        self.assertEqual( pcar.data, None )
        self.assertEqual( pcar.bands, None )
        self.assertEqual( pcar.occupancy, None )
        self.assertEqual( pcar.number_of_projections, None )
        self.assertEqual( pcar.k_point_blocks, None )
        self.assertEqual( pcar.calculation, { 'non_spin_polarised': False, 
                                              'non_collinear': False, 
                                              'spin_polarised': False } )
        
    def test_procar_is_read_from_file( self ):
        """Checking that `PROCAR_test` is read"""
        pcar = procar.Procar()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pcar.read_from_file( test_procar_filename )
        self.assertEqual( pcar.spin_channels, 4 )
        self.assertEqual( pcar.number_of_ions, 22 )
        self.assertEqual( pcar.number_of_bands, 4 )
        self.assertEqual( pcar.number_of_k_points, 2 )

    def test_procar_from_file_correctly_parses_bands( self ):
        pcar = procar.Procar()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pcar.read_from_file( test_procar_filename )
        np.testing.assert_equal( pcar.bands, 
                             np.array([[ 1., -13.17934476],
                                       [ 2., -13.17934476],
                                       [ 3., -13.16936722],
                                       [ 4., -13.16936722],
                                       [ 1., -13.1849117 ],
                                       [ 2., -13.1849117 ],
                                       [ 3., -13.16621473],
                                       [ 4., -13.16621472]] ) )

    def test_procar_from_file_correctly_parses_occupancies( self ):
        pcar = procar.Procar()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pcar.read_from_file( test_procar_filename )
        np.testing.assert_equal( pcar.occupancy,
                             np.array([[ 1., 1. ],
                                       [ 2., 1. ],
                                       [ 3., 1. ],
                                       [ 4., 1. ],
                                       [ 1., 1. ],
                                       [ 2., 1. ],
                                       [ 3., 1. ],
                                       [ 4., -0.03191968]] ) )

    def test_read_from_file_raises_valueerror_if_negative_occupancies_is_invalid( self ):
        pcar = procar.Procar()
        with self.assertRaises( ValueError ):
            pcar.read_from_file( test_procar_filename, negative_occupancies='foo' )

    def test_read_from_file_warns_about_negative_occupancies( self ):
        pcar = procar.Procar()
        with warnings.catch_warnings( record=True ) as w:
            pcar.read_from_file( test_procar_filename, negative_occupancies='warn' )
            self.assertEqual( len(w), 1 )
            self.assertTrue( "negative" in str(w[-1].message) )

    def test_read_from_file_raises_exception_if_negative_occupancies_is_raise( self ):
        pcar = procar.Procar()
        with self.assertRaises( ValueError ):
            pcar.read_from_file( test_procar_filename, negative_occupancies='raise' )

    def test_read_from_file_if_negative_occupancies_is_ignore( self ):
        pcar = procar.Procar()
        pcar.read_from_file( test_procar_filename, negative_occupancies='ignore' )

    def test_read_from_file_zeros_occupancies_if_negative_occupancies_is_zero( self ):
        pcar = procar.Procar()
        pcar.read_from_file( test_procar_filename, negative_occupancies='zero' )
        np.testing.assert_array_equal( pcar.occupancy[:,1],
            np.array( [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0 ] ) )

    def test_procar_occupation_values_are_read( self ):
        pcar = procar.Procar()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pcar.read_from_file( test_procar_filename )
            np.testing.assert_array_equal( pcar.occupancy[:,1], 
                np.array( [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -0.03191968 ] ) )

    def test_spin_polarised_procar_is_read_from_file( self ):
        """Checking that `PROCAR_spin_polarised_test` is read"""
        pcar = procar.Procar()
        pcar.read_from_file( test_procar_spin_polarised_filename )
        self.assertEqual( pcar.spin_channels, 2 )
        self.assertEqual( pcar.number_of_ions, 25 )
        self.assertEqual( pcar.number_of_bands, 112 )
        self.assertEqual( pcar.number_of_k_points, 8 )
        
    def test___add___( self ):
        pcar1 = procar.Procar()
        pcar2 = procar.Procar()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pcar1.read_from_file( test_procar_filename )
            pcar2.read_from_file( test_procar_filename )
        combined_pcar = pcar1 + pcar2
        self.assertEqual( combined_pcar.spin_channels, 4 )
        self.assertEqual( combined_pcar.number_of_ions, 22 )
        self.assertEqual( combined_pcar.number_of_bands, 4 )
        self.assertEqual( combined_pcar.number_of_k_points, 4 )
        np.testing.assert_equal( combined_pcar.occupancy, np.vstack( ( pcar1.occupancy, pcar2.occupancy ) ) )
        np.testing.assert_equal( combined_pcar.bands, np.vstack( ( pcar1.bands, pcar2.bands ) ) )
        np.testing.assert_equal( combined_pcar.k_points, np.vstack( ( pcar1.k_points, pcar2.k_points ) ) )

    def test___add___spin_polarised_procars( self ):
        pcar1 = procar.Procar()
        pcar2 = procar.Procar()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pcar1.read_from_file( test_procar_spin_polarised_filename )
            pcar2.read_from_file( test_procar_spin_polarised_filename )
        combined_pcar = pcar1 + pcar2
        self.assertEqual( combined_pcar.spin_channels, 2 )
        self.assertEqual( combined_pcar.number_of_ions, 25 )
        self.assertEqual( combined_pcar.number_of_bands, 112 )
        self.assertEqual( combined_pcar.number_of_k_points, 16 )
 
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

    def test_points_are_in_a_straight_line( self ):
        a = np.array( [ 0.0, 0.0, 0.0 ] )
        b = np.array( [ 2.0, 0.0, 0.0 ] )
        c = np.array( [ 3.0, 0.0, 0.0 ] )
        points = [ a, b, c ]
        tolerance = 1e-7
        self.assertEqual( procar.points_are_in_a_straight_line( points, tolerance ), True )

    def test_points_are_not_in_a_straight_line( self ):
        a = np.array( [ 0.0, 0.0, 0.0 ] )
        b = np.array( [ 2.0, 1.0, 0.0 ] )
        c = np.array( [ 3.0, 0.0, 0.0 ] )
        points = [ a, b, c ]
        tolerance = 1e-7
        self.assertEqual( procar.points_are_in_a_straight_line( points, tolerance ), False )

    def test_least_squares_effective_mass( self ):
        k_points = np.array( [ [ 0.0, 0.0, 0.0 ],
                               [ 1.0, 0.0, 0.0 ],
                               [ 2.0, 0.0, 0.0 ] ] )
        eigenvalues = np.array( [ 0.0, 1.0, 4.0 ] )
        with patch( 'vasppy.procar.points_are_in_a_straight_line' ) as mock_straight_line_test:
            mock_straight_line_test.return_value = True
            self.assertAlmostEqual( procar.least_squares_effective_mass( k_points, eigenvalues ), 13.6056930 )  
 
    def test_least_squares_effective_mass_raises_valueerror_if_points_are_not_collinear( self ):
        k_points = np.array( [ [ 0.0, 0.0, 0.0 ],
                               [ 1.0, 0.0, 1.0 ],
                               [ 2.0, 0.0, 0.0 ] ] )
        eigenvalues = np.array( [ 0.0, 1.0, 4.0 ] )
        with patch( 'vasppy.procar.points_are_in_a_straight_line' ) as mock_straight_line_test:
            mock_straight_line_test.return_value = False
            with self.assertRaises( ValueError ):
                procar.least_squares_effective_mass( k_points, eigenvalues )

    def test_two_point_effective_mass( self ):
        k_points = np.array( [ [ 0.0, 0.0, 0.0 ],
                               [ 1.0, 0.0, 0.0 ] ] )
        eigenvalues = np.array( [ 0.0, 1.0 ] )
        self.assertAlmostEqual( procar.two_point_effective_mass( k_points, eigenvalues ), 13.605693010 )

if __name__ == '__main__':
    unittest.main()
