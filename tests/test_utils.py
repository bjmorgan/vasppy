import unittest
import hashlib
from vasppy.utils import md5sum, file_md5, validate_checksum, dr_ij
from unittest.mock import patch, mock_open
import numpy as np
from pymatgen.core import Lattice, Structure

class UtilsTestCase( unittest.TestCase ):

    def test_md5sum( self ):
        string = 'abcdefg'
        h = hashlib.new( 'md5' )
        h.update( string.encode( 'utf-8' ) )
        self.assertEqual( md5sum( string ), h.hexdigest() )

    def test_file_md5( self ):
        example_file = "abc\nabc\n"
        with patch( 'vasppy.utils.md5sum' ) as mock_md5sum:
            mock_md5sum.return_value = 'foo'
            with patch( 'vasppy.utils.zopen', mock_open( read_data=example_file ), create=True ) as m:
                self.assertEqual( file_md5( m ), 'foo' )
                mock_md5sum.assert_called_with( example_file )

    def test_validate_checksum( self ):
        with patch( 'vasppy.utils.match_filename' ):
            with patch( 'vasppy.utils.file_md5' ) as mock_file_md5:
                mock_file_md5.return_value='abcdef'
                validate_checksum( filename='foo', md5sum='abcdef' )                

    def test_validate_checksum_raises_ValueError_if_not_matched( self ):
        with patch( 'vasppy.utils.match_filename' ) as mock_match_filename:
            mock_match_filename.return_value='bar'
            with patch( 'vasppy.utils.file_md5' ) as mock_file_md5:
                mock_file_md5.return_value='bcdefg'
                with self.assertRaises( ValueError ):
                    validate_checksum( filename='foo', md5sum='abcdef' )
                    
class Test_drij(unittest.TestCase):
    
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
               

if __name__ == '__main__':
    unittest.main()
