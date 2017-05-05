import unittest
from unittest.mock import Mock, patch, mock_open

from vasppy.vaspmeta import VASPMeta

class VASPMetaTestCase( unittest.TestCase ):

    def test_init_vaspmeta( self ):
        title = 'title'
        description = 'description'
        notes = 'notes'
        valid_status = [ 'to-run', 'incomplete', 'finished', 'dropped' ]
        for s in valid_status:
            vaspmeta = VASPMeta( title, description, status=s, notes=notes )
            self.assertEqual( vaspmeta.title, title )
            self.assertEqual( vaspmeta.description, description )
            self.assertEqual( vaspmeta.status, s )
            self.assertEqual( vaspmeta.notes, notes )
    
    def test_init_vaspmeta_with_invalid_status_raises_valueerror( self ):
        title = 'title'
        description = 'description'
        invalid_status = 'foo'
        with self.assertRaises( ValueError ):
            VASPMeta( title, description, status=invalid_status )
    
    def test_init_vaspmeta_with_no_notes( self ):
        title = 'title'
        description = 'description'
        status = 'finished'
        vaspmeta = VASPMeta( title, description, status )
        self.assertEqual( vaspmeta.notes, None )

    def test_from_file( self ):
        example_file = """\
title: title
description: description
notes: notes
status: finished\
"""
        with patch( 'vasppy.vaspmeta.VASPMeta' ) as mock_VASPMeta:
            mock_VASPMeta.return_value = 'my VASP metadata'
            with patch( 'builtins.open', mock_open( read_data=example_file), create=True ) as m:
                vaspmeta = VASPMeta.from_file( example_file )
        mock_VASPMeta.assert_called_with( 'title', 'description', 'finished', notes='notes' )
        self.assertEqual( vaspmeta, mock_VASPMeta.return_value )

    def test_from_file_with_no_notes_entry( self ):
        example_file = """\
title: title
description: description
status: finished\
"""
        with patch( 'vasppy.vaspmeta.VASPMeta' ) as mock_VASPMeta:
            mock_VASPMeta.return_value = 'my VASP metadata'
            with patch( 'builtins.open', mock_open( read_data=example_file), create=True ) as m:
                vaspmeta = VASPMeta.from_file( example_file )
        mock_VASPMeta.assert_called_with( 'title', 'description', 'finished' )
        self.assertEqual( vaspmeta, mock_VASPMeta.return_value )

if __name__ == '__main__':
    unittest.main()
