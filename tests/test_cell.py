import unittest
from vasppy.cell import angle, rotation_matrix, Cell
from unittest.mock import patch, Mock
import numpy as np

class Test_Cell( unittest.TestCase ):

    def test_cell_init( self ):
        cell_matrix = np.array( [ [ 1.0, 1.0, 0.0 ],
                                  [ 1.0, 0.0, 1.0 ],
                                  [ 0.0, 1.0, 1.0 ] ] )
        with patch( 'numpy.linalg.inv' ) as mock_invert: 
            mock_invert.return_value = np.array( [ [ 0.0, 0.0, 1.0 ],
                                                   [ 0.0, 1.0, 0.0 ],
                                                   [ 1.0, 0.0, 0.0 ] ] )
            cell = Cell( cell_matrix )
            mock_invert.assert_called_with( cell_matrix )
            np.testing.assert_array_equal( cell.matrix, cell_matrix )
            np.testing.assert_array_equal( cell.inv_matrix, mock_invert.return_value )

    def test_angle( self ):
        test_data = [ [ np.array( [ 1.0, 0.0, 0.0 ] ), np.array( [ 0.0, 1.0, 0.0 ] ), 90.0 ],
                      [ np.array( [ 2.0, 2.0, 0.0 ] ), np.array( [ 0.5, 0.0, 0.0 ] ), 45.0 ] ]
        for data in test_data:
            self.assertAlmostEqual( angle( data[0], data[1] ), data[2] )

if __name__ == '__main__':
    unittest.main() 
