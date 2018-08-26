import unittest
import hashlib
from vasppy.utils import md5sum, file_md5, validate_checksum
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

    def test_validate_checksum( self ):
        with patch( 'vasppy.utils.match_filename' ) as mock_match_filename:
            with patch( 'vasppy.utils.file_md5' ) as mock_file_md5:
                mock_file_md5.return_value='abcdef'
                validate_checksum( filename='foo', md5sum='abcdef' )                

    def test_validate_checksum_raises_ValueError_if_not_matched( self ):
        with patch( 'vasppy.utils.match_filename' ) as mock_match_filename:
            mock_match_filename.return_value='bar'
            with patch( 'vasppy.utils.file_md5' ) as mock_file_md5:
                mock_file_md5.return_value='bcdefg'
                with self.assertRaises( ValueError ):
                    validate_checksum( filename='foo', md5sum='abcdef' )                

if __name__ == '__main__':
    unittest.main()
