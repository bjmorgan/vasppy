import unittest
import numpy as np
import io
import inspect
from unittest.mock import Mock, patch, call

from vasppy.summary import Summary, md5sum, potcar_spec, find_vasp_calculations

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

    @patch('vasppy.summary.VASPMeta')
    @patch('vasppy.summary.Summary.parse_vasprun')
    def test_summary_is_initialised( self, mock_parse_vasprun, MockVASPMeta ):
        MockVASPMeta.from_file = Mock( return_value='foo' )
        summary = Summary()
        mock_parse_vasprun.assert_called_once()
        expected_print_methods = [ 'title', 'type', 'status', 'stoichiometry',
                                   'potcar', 'eatom', 'energy', 'k-points',
                                   'functional', 'encut', 'plus_u', 'ediffg',
                                   'ibrion', 'converged', 'version', 'md5',
                                   'directory', 'lreal', 'vbm', 'cbm' ]
        for key in expected_print_methods:
            self.assertTrue(key in summary.print_methods)
            self.assertTrue( inspect.ismethod( summary.print_methods[ key ] ) )

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

    def test_find_vasp_calculations( self ):
        mock_glob_output = [ 'dir_A/vasprun.xml', 'dir_B/dir_C/vasprun.xml' ]
        with patch('glob.iglob', return_value=mock_glob_output) as mock_glob:
            v = find_vasp_calculations()
        self.assertEqual( v, [ './dir_A/', './dir_B/dir_C/' ] )
 
if __name__ == '__main__':
    unittest.main()
