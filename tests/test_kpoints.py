import unittest

from vasppy.kpoints import AutoKPoints, get_convergence_testing_kspacing, get_subdivisions_from_kspacing
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

class KspacingTestCase(unittest.TestCase):

    def test_subdivisions_from_kspacing(self):
        kspacing = 0.2
        reciprocal_lattice_vectors = np.array([[0.16, 0, 0], [0, 0.19, 0], [0, 0, 0.20]])
        subdivisions = get_subdivisions_from_kspacing(kspacing, reciprocal_lattice_vectors)
        self.assertEqual(subdivisions, (6, 6, 7))

    def test_convergence_testing_kspacing(self):
        reciprocal_lattice_vectors = np.array([[0.16, 0, 0], [0, 0.19, 0], [0, 0, 0.20]])
        allowed_kspacing = get_convergence_testing_kspacing(reciprocal_lattice_vectors)
        self.assertEqual(allowed_kspacing, [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.3, 0.32, 0.34, 0.4, 0.42, 0.52, 0.6, 0.64])

if __name__ == '__main__':
    unittest.main()
