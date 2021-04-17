import unittest
from vasppy.rdf import RadialDistributionFunction
from vasppy.rdf import NeighbourList
from vasppy.rdf import dr_ij
import numpy as np
from unittest.mock import Mock, patch, call, create_autospec
from pymatgen.core import Structure, Lattice

class TestRDFFunctions(unittest.TestCase):
    
    def test_dr_ij_default(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                                      alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords)
        dr = dr_ij(structure)
        np.testing.assert_array_almost_equal(dr,
            np.array([[3.46410162], 
                      [3.46410162]]))
            
    def test_dr_ij_all_distances(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                                  alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords) 
        dr = dr_ij(structure, self_reference=True)
        np.testing.assert_array_almost_equal(dr,
            np.array([[0.0, 3.46410162],
                      [3.46410162, 0.0]]))
                      
    def test_dr_ij_indices_i_set(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                                  alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords) 
        dr = dr_ij(structure, 
                   indices_i=[0, 1])
        np.testing.assert_array_almost_equal(dr,
            np.array([[3.46410162],
                      [3.46410162]]))
                      
    def test_dr_ij_indices_i_indices_j_set(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                                  alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords) 
        dr = dr_ij(structure, 
                   indices_i=[0],
                   indices_j=[1])
        np.testing.assert_array_almost_equal(dr,
            np.array([[3.46410162]]))
            
    def test_dr_ij_indices_i_indices_j_set_as_numpy_arrays(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                                  alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords) 
        dr = dr_ij(structure, 
                   indices_i=np.array([0]),
                   indices_j=np.array([1]))
        np.testing.assert_array_almost_equal(dr,
            np.array([[3.46410162]]))
            
    def test_dr_ij_masking_with_unsorted_indices(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                              alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords) 
        dr = dr_ij(structure, 
               indices_i=[1,0],
               indices_j=[0,1])
        np.testing.assert_array_almost_equal(dr,
        np.array([[3.46410162],
                  [3.46410162]]))
                  
    def test_dr_ij_masking_with_partially_overlapping_indices(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                              alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords) 
        dr = dr_ij(structure, 
               indices_i=[0],
               indices_j=[0,1])
        np.testing.assert_array_almost_equal(dr,
        np.array([[3.46410162]]))
        
    def test_dr_ij_masking_with_partially_overlapping_indices_with_self_reference(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                              alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords) 
        dr = dr_ij(structure, 
               indices_i=[0],
               indices_j=[0,1],
               self_reference=True)
        np.testing.assert_array_almost_equal(dr,
        np.array([[0.0, 3.46410162]]))

    
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
        with patch('vasppy.rdf.dr_ij') as mock_dr_ij:
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
        expected_calls = [call(structure=mock_structures[0],
                               indices_i=[0,1],
                               indices_j=[0,1],
                               self_reference=False),
                          call(structure=mock_structures[1],
                                 indices_i=[0,1],
                                 indices_j=[0,1],
                                 self_reference=False)]
        mock_dr_ij.assert_has_calls(expected_calls)
        

class TestNeighbourList(unittest.TestCase):

    def test_NeighbourList_one(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                                  alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords)
        indices_i = [0]
        indices_j = [1]
        r_cut = 3.0
        neighbour_list = NeighbourList(structure=structure,
                                       indices_i=indices_i,
                                       indices_j=indices_j,
                                       r_cut=r_cut)
        np.testing.assert_array_equal(neighbour_list.vectors, np.array([[0]]))
        
    def test_NeighbourList_two(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                                  alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords)
        indices_i = [0,1]
        indices_j = [0,1]
        r_cut = 3.5
        neighbour_list = NeighbourList(structure=structure,
                                       indices_i=indices_i,
                                       indices_j=indices_j,
                                       r_cut=r_cut)
        np.testing.assert_array_equal(neighbour_list.vectors, np.array([[1],[1]]))
        
    def test_coordination_numbers(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                              alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords)
        indices_i = [0,1]
        indices_j = [0,1]
        r_cut = 3.5
        neighbour_list = NeighbourList(structure=structure,
                                   indices_i=indices_i,
                                   indices_j=indices_j,
                                   r_cut=r_cut)
        np.testing.assert_array_equal(neighbour_list.coordination_numbers, np.array([1, 1]))
        
    def test_from_species_strings(self):
        coords = np.array([[0.5, 0.5, 0.5],
                           [0.0, 0.0, 0.0]])
        atom_list = ['Na', 'Cl']
        lattice = Lattice.from_parameters(a=4.0, b=4.0, c=4.0, 
                              alpha=90, beta=90, gamma=90)
        structure = Structure(lattice, atom_list, coords)
        r_cut = 3.5
        with patch("vasppy.rdf.NeighbourList.__init__") as mock_init:
            mock_init.return_value = None
            nlist = NeighbourList.from_species_strings(structure=structure,
                                               species_i="Na",
                                               species_j="Cl",
                                               r_cut=r_cut)
            mock_init.assert_called_with(structure=structure,
                                         indices_i=[0],
                                         indices_j=[1],
                                         r_cut=r_cut)
            self.assertTrue(isinstance(nlist, NeighbourList))
        
        
if __name__ == '__main__':
    unittest.main()
