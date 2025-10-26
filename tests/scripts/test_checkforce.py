import unittest
from unittest.mock import patch, mock_open
import numpy as np

from vasppy.scripts.checkforce import get_forces_data, ForcesData


class GetForcesDataTestCase(unittest.TestCase):

    def test_get_forces_data_uses_forces_from_outcar(self):
        """Test that get_forces_data uses forces_from_outcar with last_one_only=True."""
        mock_forces = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
        outcar_content = "EDIFFG = -0.02\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            with patch("vasppy.scripts.checkforce.forces_from_outcar") as mock_forces_from_outcar:
                mock_forces_from_outcar.return_value = mock_forces

                forces_data = get_forces_data("test_outcar")

                # Verify forces_from_outcar was called with correct parameters
                mock_forces_from_outcar.assert_called_once_with("test_outcar", last_one_only=True)
                # Verify ForcesData was created correctly
                self.assertIsInstance(forces_data, ForcesData)
                self.assertEqual(forces_data.number_of_ions, 2)
                self.assertEqual(forces_data.convergence, 0.02)

    def test_get_forces_data_reads_ediffg_from_outcar(self):
        """Test that get_forces_data reads EDIFFG convergence criterion from OUTCAR."""
        mock_forces = np.array([[0.1, 0.2, 0.3]])
        outcar_content = "EDIFFG = -0.05\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            with patch("vasppy.scripts.checkforce.forces_from_outcar") as mock_forces_from_outcar:
                mock_forces_from_outcar.return_value = mock_forces

                forces_data = get_forces_data("test_outcar")

                self.assertEqual(forces_data.convergence, 0.05)

    def test_get_forces_data_uses_provided_convergence(self):
        """Test that get_forces_data uses provided convergence value instead of reading EDIFFG."""
        mock_forces = np.array([[0.1, 0.2, 0.3]])
        outcar_content = "EDIFFG = -0.05\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            with patch("vasppy.scripts.checkforce.forces_from_outcar") as mock_forces_from_outcar:
                mock_forces_from_outcar.return_value = mock_forces

                forces_data = get_forces_data("test_outcar", convergence=0.01)

                self.assertEqual(forces_data.convergence, 0.01)

    def test_get_forces_data_raises_error_if_ediffg_not_found(self):
        """Test that get_forces_data raises ValueError if EDIFFG cannot be read."""
        mock_forces = np.array([[0.1, 0.2, 0.3]])
        outcar_content = "No EDIFFG here\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            with patch("vasppy.scripts.checkforce.forces_from_outcar") as mock_forces_from_outcar:
                mock_forces_from_outcar.return_value = mock_forces

                with self.assertRaises(ValueError) as context:
                    get_forces_data("test_outcar")

                self.assertIn("Unable to read EDIFFG", str(context.exception))

    def test_get_forces_data_gets_number_of_ions_from_forces_shape(self):
        """Test that get_forces_data derives number_of_ions from forces array shape."""
        mock_forces = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
        outcar_content = "EDIFFG = -0.02\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            with patch("vasppy.scripts.checkforce.forces_from_outcar") as mock_forces_from_outcar:
                mock_forces_from_outcar.return_value = mock_forces

                forces_data = get_forces_data("test_outcar")

                self.assertEqual(forces_data.number_of_ions, 3)

class GetAllForcesDataTestCase(unittest.TestCase):
    """Tests for getting forces data from all ionic steps."""

    def test_get_all_forces_data_returns_generator_of_forces_data(self):
        """Test that we can get ForcesData for all ionic steps."""
        mock_all_forces = np.array([
            [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],  # Step 1
            [[0.7, 0.8, 0.9], [1.0, 1.1, 1.2]],  # Step 2
            [[1.3, 1.4, 1.5], [1.6, 1.7, 1.8]],  # Step 3
        ])
        outcar_content = "EDIFFG = -0.02\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            with patch("vasppy.scripts.checkforce.forces_from_outcar") as mock_forces_from_outcar:
                mock_forces_from_outcar.return_value = mock_all_forces

                from vasppy.scripts.checkforce import get_all_forces_data
                all_forces_data = list(get_all_forces_data("test_outcar"))

                # Should return 3 ForcesData objects
                self.assertEqual(len(all_forces_data), 3)
                for forces_data in all_forces_data:
                    self.assertIsInstance(forces_data, ForcesData)
                    self.assertEqual(forces_data.number_of_ions, 2)
                    self.assertEqual(forces_data.convergence, 0.02)

    def test_get_all_forces_data_with_provided_convergence(self):
        """Test that get_all_forces_data uses provided convergence."""
        mock_all_forces = np.array([
            [[0.1, 0.2, 0.3]],
            [[0.4, 0.5, 0.6]],
        ])
        outcar_content = "EDIFFG = -0.02\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            with patch("vasppy.scripts.checkforce.forces_from_outcar") as mock_forces_from_outcar:
                mock_forces_from_outcar.return_value = mock_all_forces

                from vasppy.scripts.checkforce import get_all_forces_data
                all_forces_data = list(get_all_forces_data("test_outcar", convergence=0.01))

                self.assertEqual(len(all_forces_data), 2)
                for forces_data in all_forces_data:
                    self.assertEqual(forces_data.convergence, 0.01)


class ParseArgsTestCase(unittest.TestCase):
    """Tests for command line argument parsing."""

    def test_parse_args_all_flag(self):
        """Test that --all flag is parsed correctly."""
        from vasppy.scripts.checkforce import parse_args
        with patch('sys.argv', ['checkforce', '--all']):
            args = parse_args()
            self.assertTrue(args.all)

    def test_parse_args_all_and_verbose_mutually_exclusive(self):
        """Test that --all and --verbose cannot be used together."""
        from vasppy.scripts.checkforce import parse_args
        with patch('sys.argv', ['checkforce', '--all', '--verbose']):
            with self.assertRaises(SystemExit):
                parse_args()

class ReadEdiffgTestCase(unittest.TestCase):
    """Tests for reading EDIFFG from OUTCAR."""

    def test_read_ediffg_from_outcar(self):
        """Test that read_ediffg_from_outcar reads EDIFFG correctly."""
        outcar_content = "EDIFFG = -0.03\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            from vasppy.scripts.checkforce import read_ediffg_from_outcar
            convergence = read_ediffg_from_outcar("test_outcar")
            self.assertEqual(convergence, 0.03)

    def test_read_ediffg_from_outcar_raises_error_if_not_found(self):
        """Test that read_ediffg_from_outcar raises ValueError if EDIFFG not found."""
        outcar_content = "No EDIFFG here\n"

        with patch("builtins.open", mock_open(read_data=outcar_content)):
            from vasppy.scripts.checkforce import read_ediffg_from_outcar
            with self.assertRaises(ValueError) as context:
                read_ediffg_from_outcar("test_outcar")
            self.assertIn("Unable to read EDIFFG", str(context.exception))


if __name__ == "__main__":
    unittest.main()
