import unittest
import io
from vasppy.doscar import Doscar
from vasppy.doscar import ispin_from_doscar
from unittest.mock import patch, mock_open

class DoscarFunctionsTestCase(unittest.TestCase):

    def test_ispin_from_doscar_returns_1_for_non_spin_polarised_input(self):
        doscar_filename = 'test_data/DOSCAR_non_spin_polarised'
        self.assertEqual(ispin_from_doscar(filename=doscar_filename), 1)

    def test_ispin_from_doscar_returns_2_for_spin_polarised_input(self):
        doscar_filename = 'test_data/DOSCAR_spin_polarised'
        self.assertEqual(ispin_from_doscar(filename=doscar_filename), 2)

    def test_ispin_from_doscar_raises_valueError_for_invalid_input(self):
        mock_file = io.StringIO('a\nb\nc\nd\ne\nf\ng')
        with patch('vasppy.doscar.open', return_value=mock_file, create=True):
            with self.assertRaises(ValueError):
                ispin_from_doscar(filename="bar")

if __name__ == '__main__':
    unittest.main()
