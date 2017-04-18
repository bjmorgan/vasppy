import unittest
from unittest.mock import Mock, patch, mock_open

from vasppy.kpoints import AutoKPoints
import numpy as np

class AutoKPointsTestCase( unittest.TestCase ):

    def test_init_auto_kpoints( self ):
        title = 'title'
        subdivisions = np.array( [ 2, 2, 2 ] )
        auto_kpoints = AutoKPoints( title, subdivisions )
        self.assertEqual( auto_kpoints.title, title )
        np.testing.assert_array_equal( auto_kpoints.subdivisions, subdivisions )
        np.testing.assert_array_equal( auto_kpoints.grid_centering, 'G' )
        np.testing.assert_array_equal( auto_kpoints.shift, np.array( [ 0.0, 0.0, 0.0 ] ) )

    def test_init_auto_kpoints_mp( self ):
        title = 'title'
        subdivisions = np.array( [ 2, 2, 2 ] )
        grid_centering = 'MP'
        auto_kpoints = AutoKPoints( title, subdivisions, grid_centering=grid_centering )
        np.testing.assert_array_equal( auto_kpoints.grid_centering, 'MP' )

    def test_init_auto_kpoints_shift( self ):
        title = 'title'
        subdivisions = np.array( [ 2, 2, 2 ] )
        shift = np.array( [ 0.1, 0.2, 0.3 ] )
        auto_kpoints = AutoKPoints( title, subdivisions, shift=shift )
        np.testing.assert_array_equal( auto_kpoints.shift, shift )

    def test_init_auto_kpoints_invalid_grid_centering_raises_valueerror( self ):
        title = 'title'
        subdivisions = np.array( [ 2, 2, 2 ] )
        grid_centering = 'foo'
        with self.assertRaises( ValueError ):
            AutoKPoints( title, subdivisions, grid_centering=grid_centering )
        
if __name__ == '__main__':
    unittest.main()
