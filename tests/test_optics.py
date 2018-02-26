import unittest
from vasppy.optics import matrix_eigvals
import numpy as np

class Test_Optics(unittest.TestCase):

    def test_matrix_eigvals(self):
        matrix = np.array( [ [ 2, 0, 3 ],
                             [ 0, 3, 0 ],
                             [ 0, 0, 3 ] ] )
        expected_eigenvalues = np.array( [ 2, 3, 3 ] )
        np.testing.assert_array_equal( matrix_eigvals( matrix ), expected_eigenvalues )

if __name__ == '__main__':
    unittest.main()
