import unittest

from vasppy.rdf import Rdf

class RdfTestCase( unittest.TestCase ):

    def test_init_rdf( self ):
        max_r = 10.0
        number_of_bins = 10
        rdf = Rdf( max_r=max_r, number_of_bins=number_of_bins )
        self.assertEqual( rdf.number_of_bins, number_of_bins )
        self.assertEqual( rdf.max_r, max_r )
        self.assertEqual( rdf.dr, 1.0 )
        # check rdf.data is zeroed

    def test_add_dr( self ):
        max_r = 10.0
        number_of_bins = 10
        rdf = Rdf( max_r=max_r, number_of_bins=number_of_bins )
        rdf.add_dr( 5.3 )
        self.assertEqual( rdf.data[5], 1 )
        rdf.add_dr( 5.0 )
        self.assertEqual( rdf.data[5], 2 )

if __name__ == '__main__':
    unittest.main()
