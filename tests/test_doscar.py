import unittest
from unittest.mock import Mock, patch, mock_open

import vasppy.doscar as doscar
import numpy as np

class TestDoscar( unittest.TestCase ):

    def test_doscar_instance_is_initialised( self ):
        with patch( 'vasppy.doscar.Doscar.read_header' ) as mock_read_header:
            mock_read_header.return_value = 'header'
            d = doscar.Doscar( 'doscar filename' )
            self.assertEqual( d.filename, 'doscar filename' )
            self.assertEqual( mock_read_header.call_count, 1 )

if __name__ == '__main__':
    unittest.main()
