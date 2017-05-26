import unittest
from unittest.mock import Mock, patch

from vasppy import xdatcar
from vasppy.poscar import Poscar

class TestXdatcar( unittest.TestCase ):

    def test_xdatcar_is_initialised( self ):
        xd = xdatcar.Xdatcar()
        self.assertEqual( type( xd.poscar[0] ), Poscar )

if __name__ == '__main__':
    unittest.main()
