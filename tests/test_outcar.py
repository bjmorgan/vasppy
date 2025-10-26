import unittest
import numpy as np
from unittest.mock import Mock, patch, mock_open, call

from vasppy.outcar import final_energy_from_outcar, potcar_eatom_list_from_outcar
from vasppy.outcar import forces_from_outcar


class OutcarTestCase(unittest.TestCase):

    def test_final_energy_from_outcar(self):
        example_file = """energy without entropy =    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy  without entropy=    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy without entropy =    -2997.63289805  energy(sigma->0) =    -2997.63289805\n"""
        with patch("builtins.open", mock_open(read_data=example_file), create=True):
            self.assertEqual(final_energy_from_outcar(), -2997.63289805)

    def test_final_energy_from_outcar_with_filename(self):
        example_file = """energy without entropy =    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy  without entropy=    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy without entropy =    -2997.63289805  energy(sigma->0) =    -2997.63289805\n"""
        with patch(
            "builtins.open", mock_open(read_data=example_file), create=True
        ) as m:
            final_energy_from_outcar(filename="foo")
            self.assertEqual(m.mock_calls[0], call("foo"))

    def test_potcar_eatom_list_from_potcar(self):
        example_file = """energy of atom  1       EATOM=-1042.3781\n
                           kinetic energy error for atom=    0.0024 (will be added to EATOM!!)\n
                           energy of atom  2       EATOM= -432.3788\n
                           kinetic energy error for atom=    0.0224 (will be added to EATOM!!)\n
                           energy of atom  3       EATOM= -659.6475\n
                           kinetic energy error for atom=    0.0354 (will be added to EATOM!!)\n"""
        with patch("builtins.open", mock_open(read_data=example_file), create=True):
            self.assertEqual(
                potcar_eatom_list_from_outcar(), [-1042.3781, -432.3788, -659.6475]
            )

    def test_forces_from_outcar_returns_all_steps_by_default(self):
        """Test that forces_from_outcar returns all ionic steps by default."""
        mock_forces = [
            [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],  # Step 1: 2 ions
            [[0.7, 0.8, 0.9], [1.0, 1.1, 1.2]],  # Step 2: 2 ions
            [[1.3, 1.4, 1.5], [1.6, 1.7, 1.8]],  # Step 3: 2 ions
        ]
        with patch("vasppy.outcar.Outcar") as mock_outcar_class:
            mock_outcar_instance = Mock()
            mock_outcar_instance.read_table_pattern.return_value = mock_forces
            mock_outcar_class.return_value = mock_outcar_instance
    
            forces = forces_from_outcar("test_outcar")
    
            # Should return all 3 steps
            self.assertEqual(forces.shape, (3, 2, 3))
            np.testing.assert_array_equal(forces, np.array(mock_forces))
            # Verify last_one_only=False was passed
            mock_outcar_instance.read_table_pattern.assert_called_once()
            call_kwargs = mock_outcar_instance.read_table_pattern.call_args[1]
            self.assertEqual(call_kwargs["last_one_only"], False)
   
    def test_forces_from_outcar_returns_only_last_step_when_requested(self):
        """Test that forces_from_outcar returns only the last step when last_one_only=True."""
        mock_last_step_forces = [[1.3, 1.4, 1.5], [1.6, 1.7, 1.8]]  # 2 ions
        
        with patch("vasppy.outcar.Outcar") as mock_outcar_class:
            mock_outcar_instance = Mock()
            # When last_one_only=True, pymatgen returns a single table, not [table]
            mock_outcar_instance.read_table_pattern.return_value = mock_last_step_forces
            mock_outcar_class.return_value = mock_outcar_instance
            
            forces = forces_from_outcar("test_outcar", last_one_only=True)
            
            # Should return only the last step (2 ions, 3 components)
            self.assertEqual(forces.shape, (2, 3))
            np.testing.assert_array_equal(forces, np.array(mock_last_step_forces))
            # Verify last_one_only=True was passed
            mock_outcar_instance.read_table_pattern.assert_called_once()
            call_kwargs = mock_outcar_instance.read_table_pattern.call_args[1]
            self.assertEqual(call_kwargs["last_one_only"], True)

    def test_forces_from_outcar_with_filename(self):
        """Test that forces_from_outcar uses the provided filename."""
        mock_forces = [[[0.1, 0.2, 0.3]]]
    
        with patch("vasppy.outcar.Outcar") as mock_outcar_class:
            mock_outcar_instance = Mock()
            mock_outcar_instance.read_table_pattern.return_value = mock_forces
            mock_outcar_class.return_value = mock_outcar_instance
    
            forces_from_outcar("custom_outcar_file")
    
            mock_outcar_class.assert_called_once_with("custom_outcar_file")


if __name__ == "__main__":
    unittest.main()
