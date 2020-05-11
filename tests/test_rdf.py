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
        with patch('vasppy.rdf.RadialDistributionFunction._RadialDistributionFunction__dr_ij') as mock_dr_ij:
            mock_dr_ij.side_effect = [np.array([5.0, 6.0]), np.array([6.0, 7.0])]
            with patch('vasppy.rdf.shell_volumes') as mock_shell_volumes:
                mock_shell_volumes.return_value = np.ones(500)
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
        expected_rdf[250] = 0.125
        expected_rdf[300] = 0.25
        expected_rdf[350] = 0.125
        expected_coordination_number = np.cumsum(expected_rdf)*2.0
        np.testing.assert_array_almost_equal(rdf.rdf, expected_rdf)
        np.testing.assert_array_almost_equal(rdf.coordination_number, expected_coordination_number)
        mock_dr_ij.assert_has_calls([call(mock_structures[0]), call(mock_structures[1])])

    def test___dr_ij(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                                          alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords) 
        mock_structures = [Mock(spec=Structure)]
        for s in mock_structures:
            s.lattice = Mock(spec=Lattice)
            s.lattice.volume = 1.0
        indices_i = [0, 1]
        with patch('vasppy.rdf.RadialDistributionFunction._RadialDistributionFunction__dr_ij') as mock_dr_ij:
            rdf = RadialDistributionFunction(structures=mock_structures,
                                             indices_i=indices_i)
        rdf.self_reference = True
        np.testing.assert_array_almost_equal(rdf._RadialDistributionFunction__dr_ij(structure), 
                                             np.array([3.46410162, 3.46410162]))
        rdf.self_reference = False
        np.testing.assert_array_almost_equal(np.sort(rdf._RadialDistributionFunction__dr_ij(structure)), 
                                             np.array([0.0, 0.0, 3.46410162, 3.46410162]))
        rdf.indices_i = [0]
        rdf.indices_j = [1]
        np.testing.assert_array_almost_equal(rdf._RadialDistributionFunction__dr_ij(structure), 
                                             np.array([3.46410162]))

if __name__ == '__main__':
    unittest.main()
