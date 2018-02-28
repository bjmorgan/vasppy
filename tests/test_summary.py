import unittest
import numpy as np
import io
from unittest.mock import Mock, patch, call

from vasppy.summary import md5sum, potcar_spec

mock_potcar_string = """foo
End of Dataset
bar
End of Dataset
sds
End of Dataset
"""

mock_potcar_data = { 'PBE':    { 'A': '12',
                                 'B': '34' },
                     'PBE_52': { 'C': '01',
                                 'D': '23' },
                     'PBE_54': { 'E': '56',
                                 'F': '78' } }

class SummaryTestCase( unittest.TestCase ):

    pass

class SummaryHelperFunctionsTestCase( unittest.TestCase ):

    def test_md5sum( self ):
        self.assertEqual( md5sum('hello\n'), 'b1946ac92492d2347c6235b4d2611184' )

    def test_potcar_spec( self ):
        mock_potcar_filename = 'POTCAR'
        md5sum_return_values = ( '12', '56', '23' ) 
        with patch('builtins.open', return_value=io.StringIO(mock_potcar_string)) as mock_open:
            with patch('vasppy.summary.md5sum', side_effect=md5sum_return_values ) as mock_md5sum:
                with patch.dict('vasppy.data.potcar_md5sum_data.potcar_md5sum_data', mock_potcar_data, clear=True ):
                    p_spec = potcar_spec( mock_potcar_filename )
                    mock_open.assert_called_with( mock_potcar_filename, 'r' )
                    mock_md5sum.assert_has_calls( [ call('foo\nEnd of Dataset\n'), 
                                                    call('bar\nEnd of Dataset\n'),
                                                    call('sds\nEnd of Dataset\n') ] )
        self.assertEqual( p_spec, {'A': 'PBE', 'E': 'PBE_54', 'D': 'PBE_52'} )

    def test_potcar_spec_raises_valueerror_if_md5sum_not_matched( self ):
        mock_potcar_filename = 'POTCAR'
        md5sum_return_values = ( '12', '56', '90' )
        with patch('builtins.open', return_value=io.StringIO(mock_potcar_string)) as mock_open:
            with patch('vasppy.summary.md5sum', side_effect=md5sum_return_values ) as mock_md5sum:
                with patch.dict('vasppy.data.potcar_md5sum_data.potcar_md5sum_data', mock_potcar_data, clear=True ):
                    with self.assertRaises( ValueError ):
                        potcar_spec( mock_potcar_filename )
 
if __name__ == '__main__':
    unittest.main()
