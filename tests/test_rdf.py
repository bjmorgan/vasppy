import unittest
from vasppy.rdf import RadialDistributionFunction
import numpy as np
from unittest.mock import Mock, patch, call
from pymatgen import Structure, Lattice

class TestRadialDistributionFunction(unittest.TestCase):

    def test_RadialDistributionFunction_raises_ValueError_if_weights_doesnt_match_structures(self):
        mock_structures = [Mock(spec=Structure), Mock(spec=Structure)]
        weights = [1, 2, 3]
        indices_i = [0, 1]
        with self.assertRaises(ValueError):
            RadialDistributionFunction(structures=mock_structures,
                                       indices_i=indices_i,
                                       weights=weights)

    def test_RadialDistributionFunction_init(self):
        mock_structures = [Mock(spec=Structure), Mock(spec=Structure)]
        for s in mock_structures:
            s.lattice = Mock(spec=Lattice)
            s.lattice.volume = 1.0
        indices_i = [0, 1]
        with patch('vasppy.rdf.RadialDistributionFunction.dr_ij') as mock_dr_ij:
            mock_dr_ij.side_effect = [np.array([5.0, 6.0]), np.array([6.0, 7.0])]
            rdf = RadialDistributionFunction(structures=mock_structures,
                                             indices_i=indices_i)
        self.assertEqual(rdf.indices_i, [0,1])
        self.assertEqual(rdf.indices_j, [0,1])
        self.assertEqual(rdf.nbins, 500)
        self.assertEqual(rdf.range, (0.0, 10.0))
        np.testing.assert_array_equal(rdf.intervals, np.linspace(0, 10, 501))
        self.assertEqual(rdf.dr, 0.02)
        np.testing.assert_array_equal(rdf.r, np.linspace(0.01, 9.99, 500))
        expected_rdf = np.zeros_like(rdf.r)
        expected_cn_basis = np.zeros_like(rdf.r)
        expected_rdf[250] = 0.019815
        expected_rdf[300] = 0.02753917
        expected_rdf[350] = 0.01012124
        expected_cn_basis[250] = 1.0
        expected_cn_basis[300] = 2.0
        expected_cn_basis[350] = 1.0
        expected_coordination_number = np.cumsum(expected_cn_basis / 2.0 / 2.0)
        np.testing.assert_array_almost_equal(rdf.rdf, expected_rdf)
        np.testing.assert_array_almost_equal(rdf.coordination_number, expected_coordination_number)
        mock_dr_ij.assert_has_calls([call(mock_structures[0]), call(mock_structures[1])])

if __name__ == '__main__':
    unittest.main()
