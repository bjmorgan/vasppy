import unittest
from vasppy.neighbour_list import NeighbourList
from vasppy.neighbour_list import neighbour_list_correlation
from vasppy.neighbour_list import neighbour_list_n_out
from vasppy.neighbour_list import neighbour_list_n_in
import numpy as np
from unittest.mock import Mock, patch, call, create_autospec
from pymatgen.core import Structure, Lattice
from copy import deepcopy
        

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
        with patch("vasppy.neighbour_list.NeighbourList.__init__") as mock_init:
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
            
    def test_neighbour_list_correlation_single_lists(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 0, 1]])
        np.testing.assert_array_equal(neighbour_list_correlation(nlist_i, nlist_j),
                                      np.array([0.5]))
        np.testing.assert_array_equal(neighbour_list_correlation(nlist_i, nlist_i),
                                      np.array([1.0]))

            
    def test_neighbour_list_correlation_multiple_lists(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0],
                                    [1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 1, 0],
                                    [1, 1, 0, 0]])
        np.testing.assert_array_equal(neighbour_list_correlation(nlist_i, nlist_j),
                                      np.array([0.5, 1.0]))
        np.testing.assert_array_equal(neighbour_list_correlation(nlist_i, nlist_i),
                                      np.array([1.0, 1.0]))
                                      
    def test_neighbour_list_correlation_different_lengths_raises_valueerror(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0, 1]])
        nlist_j.vectors = np.array([[1, 0, 0, 1]])
        with self.assertRaises(ValueError):
            neighbour_list_correlation(nlist_i, nlist_j)                                      
                                      
    def test_neighbour_list_correlation_different_vector_numbers_raises_valueerror(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0],
                                    [1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 1, 0]])
        with self.assertRaises(ValueError):
            neighbour_list_correlation(nlist_i, nlist_j)
            
    def test_neighbour_list_n_out_single_lists(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 0, 1]])
        np.testing.assert_array_equal(neighbour_list_n_out(nlist_i, nlist_j),
                                      np.array([1]))
        np.testing.assert_array_equal(neighbour_list_n_out(nlist_i, nlist_i),
                                      np.array([0]))
                                      
    def test_neighbour_list_n_out_multiple_lists(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0],
                                    [1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 1, 0],
                                    [1, 1, 0, 0]])
        np.testing.assert_array_equal(neighbour_list_n_out(nlist_i, nlist_j),
                                      np.array([1, 0]))
        np.testing.assert_array_equal(neighbour_list_n_out(nlist_i, nlist_i),
                                      np.array([0, 0]))
                                      
    def test_neighbour_list_n_out_different_vector_numbers_raises_valueerror(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0],
                                    [1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 1, 0]])
        with self.assertRaises(ValueError):
            neighbour_list_n_out(nlist_i, nlist_j)
                
    def test_neighbour_list_n_out_multiple_lists_raises_valueerror(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0, 1]])
        nlist_j.vectors = np.array([[1, 0, 0, 1]])
        with self.assertRaises(ValueError):
            neighbour_list_n_out(nlist_i, nlist_j) 
            
    def test_neighbour_list_n_in_single_lists(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 0, 1]])
        np.testing.assert_array_equal(neighbour_list_n_in(nlist_i, nlist_j),
                                      np.array([1]))
        np.testing.assert_array_equal(neighbour_list_n_in(nlist_i, nlist_i),
                                      np.array([0]))
                                      
    def test_neighbour_list_n_in_multiple_lists(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0],
                                    [1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 1, 0],
                                    [1, 1, 1, 0]])
        np.testing.assert_array_equal(neighbour_list_n_in(nlist_i, nlist_j),
                                      np.array([1, 1]))
        np.testing.assert_array_equal(neighbour_list_n_in(nlist_i, nlist_i),
                                      np.array([0, 0]))
                                      
    def test_neighbour_list_n_in_different_vector_numbers_raises_valueerror(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0],
                                    [1, 1, 0, 0]])
        nlist_j.vectors = np.array([[1, 0, 1, 0]])
        with self.assertRaises(ValueError):
            neighbour_list_n_in(nlist_i, nlist_j)
                
    def test_neighbour_list_n_in_multiple_lists_raises_valueerror(self):
        nlist_i = Mock(spec=NeighbourList)
        nlist_j = Mock(spec=NeighbourList)
        nlist_i.vectors = np.array([[1, 1, 0, 0, 1]])
        nlist_j.vectors = np.array([[1, 0, 0, 1]])
        with self.assertRaises(ValueError):
            neighbour_list_n_in(nlist_i, nlist_j)

    def test_neighbour_lists_are_equal_returns_true(self):
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
        nlist_i = neighbour_list
        nlist_j = deepcopy(neighbour_list)
        nlist_i.vectors = np.array([[1, 1, 0, 0],
                                    [1, 0, 1, 1]])
        nlist_j.vectors = np.array([[1, 1, 0, 0],
                                    [1, 0, 1, 1]])
        self.assertTrue(nlist_i == nlist_j)
        
    def test_neighbour_lists_are_not_equal_returns_false(self):
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
        nlist_i = neighbour_list
        nlist_j = deepcopy(neighbour_list)
        nlist_i.vectors = np.array([[1, 1, 0, 0],
                                    [1, 0, 1, 1]])
        nlist_j.vectors = np.array([[1, 1, 1, 0],
                                    [1, 0, 1, 1]])   
        self.assertFalse(nlist_i == nlist_j)
                                      
if __name__ == '__main__':
    unittest.main()
