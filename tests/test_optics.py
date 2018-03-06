import unittest
from vasppy.optics import matrix_eigvals, to_matrix, parse_dielectric_data
import numpy as np
from unittest.mock import Mock, patch, call

class Test_Optics(unittest.TestCase):

    def test_matrix_eigvals(self):
        matrix = np.array( [ [ 2, 0, 3 ],
                             [ 0, 3, 0 ],
                             [ 0, 0, 3 ] ] )
        expected_eigenvalues = np.array( [ 2, 3, 3 ] )
        np.testing.assert_array_equal( matrix_eigvals( matrix ), expected_eigenvalues )

    def test_to_matrix(self):
        expected_matrix = np.array( [ [ 1, 2, 3 ],
                                      [ 2, 4, 5 ],
                                      [ 3, 5, 6 ] ] )
        np.testing.assert_array_equal( to_matrix( xx=1, yy=4, zz=6, xy=2, yz=5, xz=3 ),
                                       expected_matrix )

    @patch('vasppy.optics.to_matrix')
    @patch('vasppy.optics.matrix_eigvals')
    def test_parse_dielectric_data(self, mock_matrix_eigvals, mock_to_matrix):
        input_data = [ [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ],
                       [ 1.1, 1.2, 1.3, 1.4, 1.5, 1.6 ] ]
        matrix_a = np.array( [ [ 0.1, 0.4, 0.6 ], [ 0.4, 0.2, 0.5 ], [ 0.6, 0.5, 0.3 ] ] )
        matrix_b = np.array( [ [ 1.1, 1.4, 1.6 ], [ 1.4, 1.2, 1.5 ], [ 1.6, 1.5, 1.3 ] ] )
        mock_to_matrix.side_effect = [ matrix_a, matrix_b ]
        mock_matrix_eigvals.side_effect = [ np.array( [ 1, 2, 3 ] ), np.array( [ 4, 5, 6 ] ) ]
        expected_data = np.array( [ [ 1, 2, 3 ], [ 4, 5, 6 ] ] )
        np.testing.assert_array_equal( parse_dielectric_data( input_data ), expected_data )
        mock_to_matrix.assert_has_calls( [ call( *input_data[0] ), call( *input_data[1] ) ] )
        mock_matrix_eigvals.assert_has_calls( [ call( matrix_a ), call( matrix_b ) ] )
 
if __name__ == '__main__':
    unittest.main()
