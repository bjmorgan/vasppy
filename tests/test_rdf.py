import unittest
from vasppy.rdf import RadialDistributionFunction
import numpy as np
from unittest.mock import Mock, patch, call
from pymatgen import Structure

class TestRadialDistributionFunction(unittest.TestCase):

    def test_RDF_raises_ValueError_if_weights_doesnt_match_structures(self):
        mock_structures = [Mock(spec=Structure), Mock(spec=Structure)]
        weights = [1, 2, 3]
        indices_i = [0, 1]
        with self.assertRaises(ValueError):
            RadialDistributionFunction(structures=mock_structures,
                                       indices_i=indices_i,
                                       weights=weights)

if __name__ == '__main__':
    unittest.main()
