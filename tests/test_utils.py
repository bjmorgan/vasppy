import unittest
import hashlib
from vasppy.utils import md5sum, file_md5
from unittest.mock import patch, mock_open

class UtilsTestCase( unittest.TestCase ):

    def test_md5sum( self ):
        string = 'abcdefg'
        h = hashlib.new( 'md5' )
        h.update( string.encode( 'utf-8' ) )
        self.assertEqual( md5sum( string ), h.hexdigest() )

    def test_file_md5( self ):
        example_file = "abc\nabc\n"
        with patch( 'vasppy.utils.md5sum' ) as mock_md5sum:
            mock_md5sum.return_value = 'foo'
            with patch( 'vasppy.utils.zopen', mock_open( read_data=example_file ), create=True ) as m:
                self.assertEqual( file_md5( m ), 'foo' )
                mock_md5sum.assert_called_with( example_file )

if __name__ == '__main__':
    unittest.main()
