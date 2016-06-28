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

if __name__ == '__main__':
    unittest.main()
