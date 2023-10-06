import warnings
import unittest
from unittest.mock import patch, call

from vasppy.band import Band, handle_occupancy


class BandTestCase(unittest.TestCase):
    """Tests for procar.Band class"""

    def test_band_is_initialised(self):
        """Test Band object is initialised"""
        index = 2
        energy = 1.0
        occupancy = 0.5
        with patch("vasppy.band.handle_occupancy") as mock_handle_occupancy:
            mock_handle_occupancy.return_value = 0.5
            band = Band(index=index, energy=energy, occupancy=occupancy)
        self.assertEqual(index, band.index)
        self.assertEqual(energy, band.energy)
        self.assertEqual(occupancy, band.occupancy)
        mock_handle_occupancy.assert_has_calls([call(0.5, negative_occupancies="warn")])

    def test_handle_occupancy_raises_valuerror_if_negative_occupancies_is_invalid_keyword(
        self,
    ):
        with self.assertRaises(ValueError):
            handle_occupancy(0.5, negative_occupancies="foo")

    def test_handle_occupancy_warns_about_negative_occupancies(self):
        with warnings.catch_warnings(record=True) as w:
            handle_occupancy(-0.1, negative_occupancies="warn")
            self.assertEqual(len(w), 1)
            self.assertTrue("negative" in str(w[-1].message))

    def test_handle_occupancy_raises_exception_if_negative_occupancies_is_raise(self):
        with self.assertRaises(ValueError):
            handle_occupancy(-0.1, negative_occupancies="raise")

    def test_handle_occupancy_if_negative_occupancies_is_ignore(self):
        self.assertEqual(handle_occupancy(-0.1, negative_occupancies="ignore"), -0.1)

    def test_handle_occupancy_zeros_occupancies_if_negative_occupancies_is_zero(self):
        self.assertEqual(handle_occupancy(-0.1, negative_occupancies="zero"), 0.0)


if __name__ == "__main__":
    unittest.main()
