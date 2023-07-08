import unittest
from vasppy.cell import angle, rotation_matrix, Cell
from unittest.mock import patch, Mock
import numpy as np
import math


class Test_Cell(unittest.TestCase):
    def test_cell_init(self):
        cell_matrix = np.array([[1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]])
        with patch("numpy.linalg.inv") as mock_invert:
            mock_invert.return_value = np.array(
                [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
            )
            cell = Cell(cell_matrix)
            mock_invert.assert_called_with(cell_matrix)
            np.testing.assert_array_equal(cell.matrix, cell_matrix)
            np.testing.assert_array_equal(cell.inv_matrix, mock_invert.return_value)

    def setUp(self):
        cell_matrix = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
        self.cell = Cell(cell_matrix)

    def test_dr(self):
        r1 = np.array([0.5, 0.1, 0.1])
        r2 = np.array([0.1, 0.4, 0.1])
        self.assertEqual(self.cell.dr(r1, r2), 5.0)

    def test_dr_cutoff(self):
        r1 = np.array([0.5, 0.1, 0.1])
        r2 = np.array([0.1, 0.4, 0.1])
        self.assertEqual(self.cell.dr(r1, r2, cutoff=1.0), None)

    def test_nearest_image(self):
        r1 = np.array([0.5, 0.1, 0.1])
        r2 = np.array([0.1, 0.4, 0.1])
        self.cell.minimum_image = Mock(return_value=np.array([-0.4, 0.3, 0.0]))
        np.testing.assert_array_almost_equal(self.cell.nearest_image(r1, r2), r2)

    def test_minimum_image(self):
        r1 = np.array([0.5, 0.1, 0.1])
        r2 = np.array([0.1, 0.4, 0.3])
        np.testing.assert_array_almost_equal(
            self.cell.minimum_image(r1, r2), np.array([-0.4, 0.3, 0.2])
        )

    def test_minimum_image_2(self):
        r1 = np.array([0.5, 0.1, 0.1])
        r2 = np.array([0.6, 0.8, 0.8])
        np.testing.assert_array_almost_equal(
            self.cell.minimum_image(r1, r2), np.array([0.1, -0.3, -0.3])
        )

    def test_minimum_image_dr(self):
        r1 = np.array([0.5, 0.1, 0.1])
        r2 = np.array([0.8, 0.7, 0.3])
        self.cell.minimum_image = Mock(return_value=np.array([0.3, -0.4, 0.2]))
        self.assertAlmostEqual(self.cell.minimum_image_dr(r1, r2), 5.385164807)

    def test_lengths(self):
        np.testing.assert_array_equal(self.cell.lengths(), np.array([10.0, 10.0, 10.0]))

    def test_angles(self):
        self.assertEqual(self.cell.angles(), [90.0, 90.0, 90.0])

    def test_cartesian_to_fractional_coordinates(self):
        coordinates = np.array([[1.0, 2.0, 3.0], [4.0, 6.0, 2.0]])
        expected_fractional_coordinates = np.array([[0.1, 0.2, 0.3], [0.4, 0.6, 0.2]])
        np.testing.assert_array_almost_equal(
            self.cell.cartesian_to_fractional_coordinates(coordinates),
            expected_fractional_coordinates,
        )

    def test_fractional_to_cartesian_coordinates(self):
        coordinates = np.array([[0.1, 0.2, 0.3], [0.4, 0.6, 0.2]])
        expected_cartesian_coordinates = np.array([[1.0, 2.0, 3.0], [4.0, 6.0, 2.0]])
        np.testing.assert_array_almost_equal(
            self.cell.fractional_to_cartesian_coordinates(coordinates),
            expected_cartesian_coordinates,
        )

    def test_inside_cell(self):
        c = np.array(
            [[[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]], [[0.3, 1.2, 1.3], [0.3, 0.2, 0.3]]]
        )
        for r1, r2 in c:
            np.testing.assert_array_almost_equal(self.cell.inside_cell(r1), r2)

    def test_volume(self):
        self.assertEqual(self.cell.volume(), 1000.0)

    def test_unit_vectors(self):
        np.testing.assert_array_equal(
            self.cell.unit_vectors(),
            np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        )


class Test_Cell_Support_Functions(unittest.TestCase):
    def test_angle(self):
        test_data = [
            [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), 90.0],
            [np.array([2.0, 2.0, 0.0]), np.array([0.5, 0.0, 0.0]), 45.0],
        ]
        for data in test_data:
            self.assertAlmostEqual(angle(data[0], data[1]), data[2])

    def test_rotation_matrix(self):
        test_matrix = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, math.sqrt(3) / 2, -0.5],
                [0.0, 0.5, math.sqrt(3) / 2],
            ]
        )
        angle = math.pi / 6
        axis = np.array([1.0, 0.0, 0.0])
        np.testing.assert_almost_equal(rotation_matrix(axis, angle), test_matrix)

    def test_rotation_matrix_2(self):
        test_matrix = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        angle = math.pi * 2 / 3
        axis = np.array([1.0, 1.0, 1.0])
        np.testing.assert_almost_equal(rotation_matrix(axis, angle), test_matrix)


if __name__ == "__main__":
    unittest.main()
