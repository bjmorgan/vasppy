import unittest
import numpy as np

from vasppy.summary import md5sum

class SummaryTestCase( unittest.TestCase ):

    pass

class SummaryHelperFunctionsTestCase( unittest.TestCase ):

    def test_md5sum( self ):
        self.assertEqual( md5sum('hello\n'), 'b1946ac92492d2347c6235b4d2611184' )

if __name__ == '__main__':
    unittest.main()
