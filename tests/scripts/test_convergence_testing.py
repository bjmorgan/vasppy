"""Unit tests for convergence_testing script."""
import unittest
from unittest.mock import patch, Mock, call
import numpy as np


class ParseArgsTestCase(unittest.TestCase):
    """Tests for argument parsing."""

    def test_parse_args_requires_incar_and_poscar(self):
        """Test that incar and poscar arguments are required."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', ['convergence_testing']):
            with self.assertRaises(SystemExit):
                parse_args()

    def test_parse_args_default_values(self):
        """Test that default values are set correctly."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', ['convergence_testing', '-i', 'INCAR', '-p', 'POSCAR']):
            args = parse_args()
            
            self.assertEqual(args.incar, 'INCAR')
            self.assertEqual(args.poscar, 'POSCAR')
            self.assertEqual(args.encut, (100, 700, 50))
            self.assertEqual(args.kspacing, (0.1, 0.8, 0.02))
            self.assertEqual(args.directory, './convergence_testing')
            self.assertEqual(args.base_encut, 400)
            self.assertEqual(args.base_kspacing, 0.3)
            self.assertIsNone(args.pseudopotentials)

    def test_parse_args_custom_values(self):
        """Test that custom argument values are parsed correctly."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'custom_INCAR',
            '-p', 'custom_POSCAR',
            '-e', '200', '500', '100',
            '-k', '0.2', '0.6', '0.05',
            '-d', './custom_dir',
            '--base-encut', '300',
            '--base-kspacing', '0.4',
            '--pseudopotentials', 'Li_sv', 'O'
        ]):
            args = parse_args()
            
            self.assertEqual(args.incar, 'custom_INCAR')
            self.assertEqual(args.poscar, 'custom_POSCAR')
            self.assertEqual(args.encut, ['200', '500', '100'])
            self.assertEqual(args.kspacing, ['0.2', '0.6', '0.05'])
            self.assertEqual(args.directory, './custom_dir')
            self.assertEqual(args.base_encut, '300')
            self.assertEqual(args.base_kspacing, '0.4')
            self.assertEqual(args.pseudopotentials, ['Li_sv', 'O'])


class CreateDirectoryStructureTestCase(unittest.TestCase):
    """Tests for create_directory_structure function."""

    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_creates_base_directory(self, mock_mkdir):
        """Test that base directory is created."""
        from vasppy.scripts.convergence_testing import create_directory_structure
        
        create_directory_structure('./test_dir')
        
        self.assertIn(call('./test_dir'), mock_mkdir.call_args_list)

    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_creates_encut_subdirectory(self, mock_mkdir):
        """Test that ENCUT subdirectory is created."""
        from vasppy.scripts.convergence_testing import create_directory_structure
        
        create_directory_structure('./test_dir')
        
        self.assertIn(call('./test_dir/ENCUT'), mock_mkdir.call_args_list)

    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_creates_kspacing_subdirectory(self, mock_mkdir):
        """Test that KSPACING subdirectory is created."""
        from vasppy.scripts.convergence_testing import create_directory_structure
        
        create_directory_structure('./test_dir')
        
        self.assertIn(call('./test_dir/KSPACING'), mock_mkdir.call_args_list)

    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_creates_all_directories_in_order(self, mock_mkdir):
        """Test that all three directories are created."""
        from vasppy.scripts.convergence_testing import create_directory_structure
        
        create_directory_structure('./test_dir')
        
        self.assertEqual(mock_mkdir.call_count, 3)
        expected_calls = [
            call('./test_dir'),
            call('./test_dir/ENCUT'),
            call('./test_dir/KSPACING'),
        ]
        mock_mkdir.assert_has_calls(expected_calls, any_order=False)


class WriteVaspInputFilesTestCase(unittest.TestCase):
    """Tests for write_vasp_input_files function."""

    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_writes_incar_file(self, mock_incar_class):
        """Test that INCAR file is written."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        
        incar_params = {'ENCUT': 400, 'KSPACING': 0.3}
        write_vasp_input_files('./test_dir', incar_params, mock_structure)
        
        mock_incar_class.from_dict.assert_called_once_with(incar_params)
        mock_incar_instance.write_file.assert_called_once_with('./test_dir/INCAR')

    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_writes_poscar_file(self, mock_incar_class):
        """Test that POSCAR file is written."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        
        incar_params = {'ENCUT': 400}
        write_vasp_input_files('./test_dir', incar_params, mock_structure)
        
        mock_structure.to.assert_called_once_with(fmt='poscar', filename='./test_dir/POSCAR')

    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_writes_potcar_when_provided(self, mock_incar_class):
        """Test that POTCAR is written when potcar object is provided."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        mock_potcar = Mock()
        
        incar_params = {'ENCUT': 400}
        write_vasp_input_files('./test_dir', incar_params, mock_structure, potcar=mock_potcar)
        
        mock_potcar.write_file.assert_called_once_with('./test_dir/POTCAR')

    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_does_not_write_potcar_when_not_provided(self, mock_incar_class):
        """Test that POTCAR is not written when potcar object is None."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        
        incar_params = {'ENCUT': 400}
        write_vasp_input_files('./test_dir', incar_params, mock_structure, potcar=None)
        
        # Just verify no exception is raised


class GenerateEncutInputsTestCase(unittest.TestCase):
    """Tests for generate_encut_inputs function."""

    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_creates_directory_for_each_encut_value(self, mock_mkdir, mock_write_files):
        """Test that a directory is created for each ENCUT value."""
        from vasppy.scripts.convergence_testing import generate_encut_inputs
        
        encut_values = np.array([100.0, 200.0, 300.0])
        mock_structure = Mock()
        incar_dict = {'ALGO': 'Normal'}
        
        generate_encut_inputs('./test_dir', encut_values, 0.3, incar_dict, mock_structure)
        
        expected_calls = [
            call('./test_dir/ENCUT/100.0'),
            call('./test_dir/ENCUT/200.0'),
            call('./test_dir/ENCUT/300.0'),
        ]
        mock_mkdir.assert_has_calls(expected_calls, any_order=False)

    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_writes_files_with_correct_encut_parameters(self, mock_mkdir, mock_write_files):
        """Test that files are written with correct ENCUT values."""
        from vasppy.scripts.convergence_testing import generate_encut_inputs
        
        encut_values = np.array([100.0, 200.0])
        mock_structure = Mock()
        incar_dict = {'ALGO': 'Normal'}
        
        generate_encut_inputs('./test_dir', encut_values, 0.3, incar_dict, mock_structure)
        
        # Check first call
        first_call_params = mock_write_files.call_args_list[0][0][1]
        self.assertEqual(first_call_params['ENCUT'], 100.0)
        self.assertEqual(first_call_params['KSPACING'], 0.3)
        
        # Check second call
        second_call_params = mock_write_files.call_args_list[1][0][1]
        self.assertEqual(second_call_params['ENCUT'], 200.0)
        self.assertEqual(second_call_params['KSPACING'], 0.3)

    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_preserves_base_incar_parameters(self, mock_mkdir, mock_write_files):
        """Test that base INCAR parameters are preserved."""
        from vasppy.scripts.convergence_testing import generate_encut_inputs
        
        encut_values = np.array([100.0])
        mock_structure = Mock()
        incar_dict = {'ALGO': 'Normal', 'EDIFF': 1e-6}
        
        generate_encut_inputs('./test_dir', encut_values, 0.3, incar_dict, mock_structure)
        
        call_params = mock_write_files.call_args_list[0][0][1]
        self.assertEqual(call_params['ALGO'], 'Normal')
        self.assertEqual(call_params['EDIFF'], 1e-6)

    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_passes_potcar_to_write_function(self, mock_mkdir, mock_write_files):
        """Test that potcar is passed to write function."""
        from vasppy.scripts.convergence_testing import generate_encut_inputs
        
        encut_values = np.array([100.0])
        mock_structure = Mock()
        mock_potcar = Mock()
        incar_dict = {'ALGO': 'Normal'}
        
        generate_encut_inputs('./test_dir', encut_values, 0.3, incar_dict, mock_structure, mock_potcar)
        
        # Check potcar was passed
        self.assertEqual(mock_write_files.call_args_list[0][0][3], mock_potcar)


class GenerateKspacingInputsTestCase(unittest.TestCase):
    """Tests for generate_kspacing_inputs function."""

    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_creates_directory_for_each_kspacing_value(self, mock_mkdir, mock_write_files):
        """Test that a directory is created for each KSPACING value."""
        from vasppy.scripts.convergence_testing import generate_kspacing_inputs
        
        kspacing_values = (0.1, 0.2, 0.3)
        mock_structure = Mock()
        incar_dict = {'ALGO': 'Normal'}
        
        generate_kspacing_inputs('./test_dir', kspacing_values, 400, incar_dict, mock_structure)
        
        expected_calls = [
            call('./test_dir/KSPACING/0.1'),
            call('./test_dir/KSPACING/0.2'),
            call('./test_dir/KSPACING/0.3'),
        ]
        mock_mkdir.assert_has_calls(expected_calls, any_order=False)

    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_writes_files_with_correct_kspacing_parameters(self, mock_mkdir, mock_write_files):
        """Test that files are written with correct KSPACING values."""
        from vasppy.scripts.convergence_testing import generate_kspacing_inputs
        
        kspacing_values = (0.1, 0.2)
        mock_structure = Mock()
        incar_dict = {'ALGO': 'Normal'}
        
        generate_kspacing_inputs('./test_dir', kspacing_values, 400, incar_dict, mock_structure)
        
        # Check first call
        first_call_params = mock_write_files.call_args_list[0][0][1]
        self.assertEqual(first_call_params['KSPACING'], 0.1)
        self.assertEqual(first_call_params['ENCUT'], 400)
        
        # Check second call
        second_call_params = mock_write_files.call_args_list[1][0][1]
        self.assertEqual(second_call_params['KSPACING'], 0.2)
        self.assertEqual(first_call_params['ENCUT'], 400)

    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_preserves_base_incar_parameters(self, mock_mkdir, mock_write_files):
        """Test that base INCAR parameters are preserved."""
        from vasppy.scripts.convergence_testing import generate_kspacing_inputs
        
        kspacing_values = (0.1,)
        mock_structure = Mock()
        incar_dict = {'ALGO': 'Normal', 'EDIFF': 1e-6}
        
        generate_kspacing_inputs('./test_dir', kspacing_values, 400, incar_dict, mock_structure)
        
        call_params = mock_write_files.call_args_list[0][0][1]
        self.assertEqual(call_params['ALGO'], 'Normal')
        self.assertEqual(call_params['EDIFF'], 1e-6)

    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_passes_potcar_to_write_function(self, mock_mkdir, mock_write_files):
        """Test that potcar is passed to write function."""
        from vasppy.scripts.convergence_testing import generate_kspacing_inputs
        
        kspacing_values = (0.1,)
        mock_structure = Mock()
        mock_potcar = Mock()
        incar_dict = {'ALGO': 'Normal'}
        
        generate_kspacing_inputs('./test_dir', kspacing_values, 400, incar_dict, mock_structure, mock_potcar)
        
        # Check potcar was passed
        self.assertEqual(mock_write_files.call_args_list[0][0][3], mock_potcar)


class MainIntegrationTestCase(unittest.TestCase):
    """Integration tests for main function."""

    @patch('vasppy.scripts.convergence_testing.generate_kspacing_inputs')
    @patch('vasppy.scripts.convergence_testing.generate_encut_inputs')
    @patch('vasppy.scripts.convergence_testing.create_directory_structure')
    @patch('vasppy.scripts.convergence_testing.Potcar')
    @patch('vasppy.scripts.convergence_testing.get_convergence_testing_kspacing')
    @patch('vasppy.scripts.convergence_testing.Incar')
    @patch('vasppy.scripts.convergence_testing.Structure')
    @patch('vasppy.scripts.convergence_testing.parse_args')
    def test_main_orchestrates_all_functions(
        self, mock_parse_args, mock_structure_class, mock_incar_class,
        mock_get_kspacing, mock_potcar_class, mock_create_dirs,
        mock_gen_encut, mock_gen_kspacing
    ):
        """Test that main calls all functions in correct order."""
        from vasppy.scripts.convergence_testing import main
        
        # Setup mocks
        mock_args = Mock(
            poscar='POSCAR',
            incar='INCAR',
            encut=(100, 200, 100),
            kspacing=(0.1, 0.2, 0.1),
            directory='./test_dir',
            base_encut=400,
            base_kspacing=0.3,
            pseudopotentials=None
        )
        mock_parse_args.return_value = mock_args
        
        mock_structure = Mock()
        mock_structure.lattice.reciprocal_lattice_crystallographic.matrix = \
            np.array([[0.16, 0, 0], [0, 0.19, 0], [0, 0, 0.20]])
        mock_structure_class.from_file.return_value = mock_structure
        
        mock_incar = Mock()
        mock_incar.as_dict.return_value = {'ALGO': 'Normal'}
        mock_incar_class.from_file.return_value = mock_incar
        
        mock_get_kspacing.return_value = (0.1, 0.2)
        
        # Execute
        main()
        
        # Verify all functions were called
        mock_create_dirs.assert_called_once_with('./test_dir')
        mock_gen_encut.assert_called_once()
        mock_gen_kspacing.assert_called_once()

    @patch('vasppy.scripts.convergence_testing.generate_kspacing_inputs')
    @patch('vasppy.scripts.convergence_testing.generate_encut_inputs')
    @patch('vasppy.scripts.convergence_testing.create_directory_structure')
    @patch('vasppy.scripts.convergence_testing.Potcar')
    @patch('vasppy.scripts.convergence_testing.get_convergence_testing_kspacing')
    @patch('vasppy.scripts.convergence_testing.Incar')
    @patch('vasppy.scripts.convergence_testing.Structure')
    @patch('vasppy.scripts.convergence_testing.parse_args')
    def test_main_creates_potcar_when_pseudopotentials_specified(
        self, mock_parse_args, mock_structure_class, mock_incar_class,
        mock_get_kspacing, mock_potcar_class, mock_create_dirs,
        mock_gen_encut, mock_gen_kspacing
    ):
        """Test that main creates Potcar when pseudopotentials are specified."""
        from vasppy.scripts.convergence_testing import main
        
        mock_args = Mock(
            poscar='POSCAR',
            incar='INCAR',
            encut=(100, 200, 100),
            kspacing=(0.1, 0.2, 0.1),
            directory='./test_dir',
            base_encut=400,
            base_kspacing=0.3,
            pseudopotentials=['Li_sv', 'O']
        )
        mock_parse_args.return_value = mock_args
        
        mock_structure = Mock()
        mock_structure.lattice.reciprocal_lattice_crystallographic.matrix = \
            np.array([[0.16, 0, 0], [0, 0.19, 0], [0, 0, 0.20]])
        mock_structure_class.from_file.return_value = mock_structure
        
        mock_incar = Mock()
        mock_incar.as_dict.return_value = {'ALGO': 'Normal'}
        mock_incar_class.from_file.return_value = mock_incar
        
        mock_get_kspacing.return_value = (0.1, 0.2)
        
        main()
        
        mock_potcar_class.assert_called_once_with(['Li_sv', 'O'])


if __name__ == '__main__':
    unittest.main()
