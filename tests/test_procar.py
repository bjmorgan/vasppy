import unittest
import os
from vasppy import procar
import numpy as np
from unittest.mock import patch, call
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
        self.assertEqual( pcar._number_of_k_points, None )
        self.assertEqual( pcar._number_of_bands, None )
        self.assertEqual( pcar._spin_channels, 1 )
        self.assertEqual( pcar._number_of_ions, None )
        self.assertEqual( pcar._number_of_projections, None )
        self.assertEqual( pcar._k_point_blocks, None )
        self.assertEqual( pcar.data, None )
        self.assertEqual( pcar.bands, None )
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
        np.testing.assert_equal( [ b.index for b in pcar.bands ], 
                                 [ 1., 2., 3., 4., 1., 2., 3., 4. ] )
        np.testing.assert_equal( [ b.energy for b in pcar.bands ], 
                                 [ -13.17934476,
                                   -13.17934476,
                                   -13.16936722,
                                   -13.16936722,
                                   -13.1849117 ,
                                   -13.1849117 ,
                                   -13.16621473,
                                   -13.16621472] )
        np.testing.assert_equal( [ b.occupancy for b in pcar.bands ],
                                 [ 1., 1., 1., 1., 1., 1., 1., -0.03191968 ] )

    def test_band_init( self ):
        index = 2
        energy = 1.0
        occupancy = 0.5
        with patch( 'vasppy.procar.handle_occupancy' ) as mock_handle_occupancy:
            mock_handle_occupancy.return_value = 0.5
            band = procar.Band( index=index, energy=energy, occupancy=occupancy )
        self.assertEqual( index, band.index )
        self.assertEqual( energy, band.energy )
        self.assertEqual( occupancy, band.occupancy )
        mock_handle_occupancy.assert_has_calls( [call(0.5, negative_occupancies='warn')] )
 
    def test_handle_occupancy_raises_valuerror_if_negative_occupancies_is_invalid_keyword( self ):
        with self.assertRaises( ValueError ):
            procar.handle_occupancy( 0.5, negative_occupancies='foo' ) 

    def test_handle_occupancy_warns_about_negative_occupancies( self ):
        with warnings.catch_warnings( record=True ) as w:
            procar.handle_occupancy( -0.1, negative_occupancies='warn' )
            self.assertEqual( len(w), 1 )
            self.assertTrue( "negative" in str(w[-1].message) )

    def test_handle_occupancy_raises_exception_if_negative_occupancies_is_raise( self ):
        with self.assertRaises( ValueError ):
            procar.handle_occupancy( -0.1, negative_occupancies='raise' )

    def test_handle_occupancy_if_negative_occupancies_is_ignore( self ):
        self.assertEqual( procar.handle_occupancy( -0.1, negative_occupancies='ignore' ), -0.1 )

    def test_handle_occupancy_zeros_occupancies_if_negative_occupancies_is_zero( self ):
        self.assertEqual( procar.handle_occupancy( -0.1, negative_occupancies='zero' ), 0.0 )

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
        np.testing.assert_equal( combined_pcar.bands, pcar1.bands + pcar2.bands )
        for k1, k2 in zip( combined_pcar.k_points, pcar1.k_points + pcar2.k_points ):
            np.testing.assert_equal( k1.frac_coords, k2.frac_coords )
            self.assertEqual( k1.weight, k2.weight )
        self.assertEqual( [ k.index for k in combined_pcar.k_points ], [ 1, 2, 3, 4 ] )

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
    """Tests for VASP PROCAR parsers"""

    def test_k_points_are_parsed( self ):
        """Checking that k-points are parsed from PROCAR format strings"""
        procar_string = " k-point    1 :    0.50000000 0.25000000 0.75000000     weight = 0.00806452\nk-point    2 :    0.50000000 0.25735294 0.74264706     weight = 0.00806452"
        k_points = procar.k_point_parser( procar_string )
        np.testing.assert_array_equal( [ k.frac_coords for k in k_points ],
            [ [ 0.50000000, 0.25000000, 0.75000000 ], 
              [ 0.50000000, 0.25735294, 0.74264706 ] ] )
        self.assertEqual( [ k.weight for k in k_points ],
            [ 0.00806452, 0.00806452 ] )

    def test_negative_k_points_are_parsed( self ): 
        """Checking that negative k-points are parsed from PROCAR format strings"""
        procar_string = "  k-point  119 :   -0.01282051 0.00000000 0.00000000     weight = 0.00500000\n k-point  122 :    0.00000000-0.01282051 0.01282051     weight = 0.00500000\n k-point    1 :   -0.50000000 0.00000000-0.50000000     weight = 0.00500000"
        k_points = procar.k_point_parser( procar_string )
        np.testing.assert_array_equal( [ k.frac_coords for k in k_points ],
            [ [-0.01282051,  0.00000000,  0.00000000], 
              [ 0.00000000, -0.01282051,  0.01282051],
              [-0.50000000,  0.00000000, -0.50000000]  ] )
        self.assertEqual( [ k.weight for k in k_points ],
            [ 0.005, 0.005, 0.005 ] )

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
