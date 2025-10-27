"""Unit tests for convergence_testing script."""
import unittest
from unittest.mock import patch, Mock, call
import numpy as np
from io import StringIO
from pathlib import Path
from vasppy.scripts.convergence_testing import create_directory_structure
from vasppy.scripts.convergence_testing import main
from vasppy.scripts.convergence_testing import ConvergenceTarget


class ExecuteTargetsTestCase(unittest.TestCase):
    """Tests for execute_targets function."""
    
    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_creates_directories_for_all_targets(self, mock_mkdir, mock_write_files):
        """Test that directories are created for each target."""
        from vasppy.scripts.convergence_testing import execute_targets, ConvergenceTarget
        
        targets = [
            ConvergenceTarget(
                path=Path('test_dir/ENCUT/100'),
                encut=100,
                kspacing=None,
                base_incar={'ALGO': 'Normal'}
            ),
            ConvergenceTarget(
                path=Path('test_dir/ENCUT/200'),
                encut=200,
                kspacing=None,
                base_incar={'ALGO': 'Normal'}
            ),
        ]
        mock_structure = Mock()
        mock_potcar = Mock()
        
        execute_targets(targets, mock_structure, mock_potcar)
        
        expected_calls = [
            call('test_dir/ENCUT/100'),
            call('test_dir/ENCUT/200'),
        ]
        mock_mkdir.assert_has_calls(expected_calls)
    
    @patch('vasppy.scripts.convergence_testing.write_vasp_input_files')
    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_writes_files_for_all_targets(self, mock_mkdir, mock_write_files):
        """Test that files are written for each target."""
        from vasppy.scripts.convergence_testing import execute_targets, ConvergenceTarget
        
        targets = [
            ConvergenceTarget(
                path=Path('test_dir/ENCUT/100'),
                encut=100,
                kspacing=None,
                base_incar={'ALGO': 'Normal'}
            ),
            ConvergenceTarget(
                path=Path('test_dir/ENCUT/200'),
                encut=200,
                kspacing=None,
                base_incar={'ALGO': 'Normal'}
            ),
        ]
        mock_structure = Mock()
        mock_potcar = Mock()
        
        execute_targets(targets, mock_structure, mock_potcar)
        
        self.assertEqual(mock_write_files.call_count, 2)
        # Check first call uses incar_params from target
        call_path = mock_write_files.call_args_list[0][0][0]
        call_params = mock_write_files.call_args_list[0][0][1]
        self.assertEqual(call_path, Path('test_dir/ENCUT/100'))
        self.assertEqual(call_params['ENCUT'], 100)
    
class ParseArgsTestCase(unittest.TestCase):
    """Tests for argument parsing."""
    
    def test_parse_args_default_values(self):
        """Test default argument values."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR',
            '--potcar-file', 'POTCAR'
        ]):
            args = parse_args()
            
            self.assertEqual(args.encut, (100, 700, 50))  # List, not tuple
            self.assertEqual(args.kspacing, (0.1, 0.8, 0.02))  # List, not tuple
            self.assertEqual(args.directory, './convergence_testing')
            self.assertEqual(args.base_encut, 400)
            self.assertEqual(args.base_kspacing, 0.3)
    
    def test_parse_args_custom_values(self):
        """Test parsing custom argument values."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'my_incar',
            '-p', 'my_poscar',
            '--pseudopotentials', 'Li_sv', 'O',
            '-e', '200', '800', '100',
            '-k', '0.2', '0.6', '0.05',
            '-d', './my_tests',
            '--base-encut', '500',
            '--base-kspacing', '0.4'
        ]):
            args = parse_args()
            
            self.assertEqual(args.incar, 'my_incar')
            self.assertEqual(args.poscar, 'my_poscar')
            self.assertEqual(args.encut, (200, 800, 100))
            self.assertEqual(args.kspacing, (0.2, 0.6, 0.05))
            self.assertEqual(args.directory, './my_tests')
            self.assertEqual(args.base_encut, 500)
            self.assertEqual(args.base_kspacing, 0.4)
    
    def test_parse_args_accepts_dry_run_flag(self):
        """Test that --dry-run flag is accepted."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR',
            '--potcar-file', 'POTCAR',
            '--dry-run'
        ]):
            args = parse_args()
            
            self.assertTrue(args.dry_run)
    
    def test_parse_args_dry_run_defaults_to_false(self):
        """Test that dry_run defaults to False."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR',
            '--potcar-file', 'POTCAR'
        ]):
            args = parse_args()
            
            self.assertFalse(args.dry_run)
    
    def test_parse_args_accepts_potcar_file(self):
        """Test that --potcar-file argument is accepted."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR',
            '--potcar-file', 'POTCAR'
        ]):
            args = parse_args()
            
            self.assertEqual(args.potcar_file, 'POTCAR')
            self.assertIsNone(args.pseudopotentials)
    
    def test_parse_args_accepts_pseudopotentials(self):
        """Test that --pseudopotentials argument is accepted."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR',
            '--pseudopotentials', 'Li_sv', 'O'
        ]):
            args = parse_args()
            
            self.assertEqual(args.pseudopotentials, ['Li_sv', 'O'])
            self.assertIsNone(args.potcar_file)
    
    def test_parse_args_rejects_both_pseudopotentials_and_potcar_file(self):
        """Test that --pseudopotentials and --potcar-file are mutually exclusive."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR',
            '--pseudopotentials', 'Li_sv', 'O',
            '--potcar-file', 'POTCAR'
        ]):
            with self.assertRaises(SystemExit):
                parse_args()
    
    def test_parse_args_rejects_missing_potcar_arguments(self):
        """Test that at least one POTCAR argument is required."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR'
        ]):
            with self.assertRaises(SystemExit):
                parse_args()
                
    def test_parse_args_accepts_job_script(self):
        """Test that --job-script argument is accepted."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR',
            '--potcar-file', 'POTCAR',
            '--job-script', 'submit.sh'
        ]):
            args = parse_args()
            
            self.assertEqual(args.job_script, 'submit.sh')
    
    def test_parse_args_job_script_defaults_to_none(self):
        """Test that job_script defaults to None."""
        from vasppy.scripts.convergence_testing import parse_args
        
        with patch('sys.argv', [
            'convergence_testing',
            '-i', 'INCAR',
            '-p', 'POSCAR',
            '--potcar-file', 'POTCAR'
        ]):
            args = parse_args()
            
            self.assertIsNone(args.job_script)

class CreateDirectoryStructureTestCase(unittest.TestCase):
    """Tests for create_directory_structure function."""

    @patch('vasppy.scripts.convergence_testing.os.mkdir')
    def test_creates_all_directories_in_order(self, mock_mkdir):
        """Test that all three directories are created."""
        
        create_directory_structure('./test_dir')
        
        self.assertEqual(mock_mkdir.call_count, 3)
        expected_calls = [
            call('test_dir'),
            call('test_dir/ENCUT'),
            call('test_dir/KSPACING'),
        ]
        mock_mkdir.assert_has_calls(expected_calls, any_order=False)


class WriteVaspInputFilesTestCase(unittest.TestCase):
    """Tests for write_vasp_input_files function."""

    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_writes_incar_file(self, mock_incar_class):
        """Test that INCAR file is written."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_potcar = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        
        incar_params = {'ENCUT': 400, 'KSPACING': 0.3}
        write_vasp_input_files('./test_dir', incar_params, mock_structure, mock_potcar)
        
        mock_incar_class.from_dict.assert_called_once_with(incar_params)
        mock_incar_instance.write_file.assert_called_once_with('test_dir/INCAR')
    
    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_writes_poscar_file(self, mock_incar_class):
        """Test that POSCAR file is written."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_potcar = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        
        incar_params = {'ENCUT': 400}
        write_vasp_input_files('./test_dir', incar_params, mock_structure, mock_potcar)
        
        mock_structure.to.assert_called_once_with(fmt='poscar', filename='test_dir/POSCAR')

    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_writes_potcar_file(self, mock_incar_class):
        """Test that POTCAR file is written."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        mock_potcar = Mock()
        
        incar_params = {'ENCUT': 400}
        write_vasp_input_files('./test_dir', incar_params, mock_structure, mock_potcar)
        
        mock_potcar.write_file.assert_called_once_with('test_dir/POTCAR')

    @patch('vasppy.scripts.convergence_testing.shutil')
    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_copies_job_script_when_provided(self, mock_incar_class, mock_shutil):
        """Test that job script is copied when provided."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_potcar = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        
        incar_params = {'ENCUT': 400}
        write_vasp_input_files('./test_dir', incar_params, mock_structure, mock_potcar, job_script='submit.sh')
        
        mock_shutil.copy2.assert_called_once_with('submit.sh', 'test_dir/submit.sh')
    
    @patch('vasppy.scripts.convergence_testing.shutil')
    @patch('vasppy.scripts.convergence_testing.Incar')
    def test_does_not_copy_job_script_when_none(self, mock_incar_class, mock_shutil):
        """Test that job script is not copied when None."""
        from vasppy.scripts.convergence_testing import write_vasp_input_files
        
        mock_structure = Mock()
        mock_potcar = Mock()
        mock_incar_instance = Mock()
        mock_incar_class.from_dict.return_value = mock_incar_instance
        
        incar_params = {'ENCUT': 400}
        write_vasp_input_files('./test_dir', incar_params, mock_structure, mock_potcar, job_script=None)
        
        mock_shutil.copy2.assert_not_called()


class MainTestCase(unittest.TestCase):
    """Unit tests for main function."""
    
    @patch('vasppy.scripts.convergence_testing.execute_convergence_tests')
    @patch('vasppy.scripts.convergence_testing.print_dry_run_summary')
    @patch('vasppy.scripts.convergence_testing.calculate_kspacing_targets')
    @patch('vasppy.scripts.convergence_testing.calculate_encut_targets')
    @patch('vasppy.scripts.convergence_testing.get_convergence_testing_kspacing')
    @patch('vasppy.scripts.convergence_testing.load_potcar')
    @patch('vasppy.scripts.convergence_testing.load_inputs')
    @patch('vasppy.scripts.convergence_testing.validate_inputs')
    @patch('vasppy.scripts.convergence_testing.parse_args')
    def test_main_orchestrates_all_functions(
        self, mock_parse_args, mock_validate_inputs, mock_load_inputs, mock_load_potcar,
        mock_get_kspacing, mock_calc_encut, mock_calc_kspacing, mock_print_dry_run, mock_execute
    ):
        """Test that main calls all functions in correct order."""
        from vasppy.scripts.convergence_testing import main
        
        mock_args = Mock(
            poscar='POSCAR',
            incar='INCAR',
            encut=(100, 200, 100),
            kspacing=(0.1, 0.2, 0.1),
            directory='./test_dir',
            base_encut=400,
            base_kspacing=0.3,
            pseudopotentials=None,
            potcar_file='POTCAR',
            job_script=None,
            dry_run=False
        )
        mock_parse_args.return_value = mock_args
        
        mock_structure = Mock()
        mock_structure.lattice.reciprocal_lattice_crystallographic.matrix = \
            np.array([[0.16, 0, 0], [0, 0.19, 0], [0, 0, 0.20]])
        mock_load_inputs.return_value = (mock_structure, {'ALGO': 'Normal'})
        mock_load_potcar.return_value = Mock()
        
        mock_get_kspacing.return_value = (0.1, 0.2)
        mock_calc_encut.return_value = [Mock(), Mock()]
        mock_calc_kspacing.return_value = [Mock(), Mock()]
        
        main()
        
        # Verify orchestration
        mock_validate_inputs.assert_called_once_with('POSCAR', 'INCAR', './test_dir', job_script=None)
        mock_load_inputs.assert_called_once_with('POSCAR', 'INCAR')
        mock_load_potcar.assert_called_once_with(None, 'POTCAR')
        mock_calc_encut.assert_called_once()
        mock_calc_kspacing.assert_called_once()
        mock_execute.assert_called_once()
    
    @patch('vasppy.scripts.convergence_testing.execute_convergence_tests')
    @patch('vasppy.scripts.convergence_testing.print_dry_run_summary')
    @patch('vasppy.scripts.convergence_testing.calculate_kspacing_targets')
    @patch('vasppy.scripts.convergence_testing.calculate_encut_targets')
    @patch('vasppy.scripts.convergence_testing.get_convergence_testing_kspacing')
    @patch('vasppy.scripts.convergence_testing.load_potcar')
    @patch('vasppy.scripts.convergence_testing.load_inputs')
    @patch('vasppy.scripts.convergence_testing.validate_inputs')
    @patch('vasppy.scripts.convergence_testing.parse_args')
    def test_main_passes_potcar_to_execute(
        self, mock_parse_args, mock_validate_inputs, mock_load_inputs, mock_load_potcar,
        mock_get_kspacing, mock_calc_encut, mock_calc_kspacing, mock_print_dry_run, mock_execute
    ):
        """Test that main passes potcar from load_potcar to execute_convergence_tests."""
        from vasppy.scripts.convergence_testing import main
        
        mock_args = Mock(
            poscar='POSCAR',
            incar='INCAR',
            encut=(100, 200, 100),
            kspacing=(0.1, 0.2, 0.1),
            directory='./test_dir',
            base_encut=400,
            base_kspacing=0.3,
            pseudopotentials=['Li_sv', 'O'],
            potcar_file=None,
            job_script=None,
            dry_run=False
        )
        mock_parse_args.return_value = mock_args
        
        mock_structure = Mock()
        mock_structure.lattice.reciprocal_lattice_crystallographic.matrix = \
            np.array([[0.16, 0, 0], [0, 0.19, 0], [0, 0, 0.20]])
        mock_load_inputs.return_value = (mock_structure, {'ALGO': 'Normal'})
        
        mock_potcar = Mock()
        mock_load_potcar.return_value = mock_potcar
        
        mock_get_kspacing.return_value = (0.1, 0.2)
        mock_encut_targets = [Mock()]
        mock_kspacing_targets = [Mock()]
        mock_calc_encut.return_value = mock_encut_targets
        mock_calc_kspacing.return_value = mock_kspacing_targets
        
        main()
        
        # Check load_potcar was called with correct args
        mock_load_potcar.assert_called_once_with(['Li_sv', 'O'], None)
        
        # Check execute_convergence_tests was called with potcar and job_script
        mock_execute.assert_called_once_with(
            './test_dir',
            mock_encut_targets,
            mock_kspacing_targets,
            mock_structure,
            mock_potcar,
            None  # job_script
        )
    
    @patch('vasppy.scripts.convergence_testing.execute_convergence_tests')
    @patch('vasppy.scripts.convergence_testing.print_dry_run_summary')
    @patch('vasppy.scripts.convergence_testing.calculate_kspacing_targets')
    @patch('vasppy.scripts.convergence_testing.calculate_encut_targets')
    @patch('vasppy.scripts.convergence_testing.get_convergence_testing_kspacing')
    @patch('vasppy.scripts.convergence_testing.load_potcar')
    @patch('vasppy.scripts.convergence_testing.load_inputs')
    @patch('vasppy.scripts.convergence_testing.validate_inputs')
    @patch('vasppy.scripts.convergence_testing.parse_args')
    def test_main_uses_dry_run_when_flag_set(
        self, mock_parse_args, mock_validate_inputs, mock_load_inputs, mock_load_potcar,
        mock_get_kspacing, mock_calc_encut, mock_calc_kspacing, mock_print_dry_run, mock_execute
    ):
        """Test that main uses print_dry_run_summary when dry_run is True."""
        from vasppy.scripts.convergence_testing import main
        
        mock_args = Mock(
            poscar='POSCAR',
            incar='INCAR',
            encut=(100, 200, 100),
            kspacing=(0.1, 0.2, 0.1),
            directory='./test_dir',
            base_encut=400,
            base_kspacing=0.3,
            pseudopotentials=None,
            potcar_file='POTCAR',
            job_script=None,
            dry_run=True
        )
        mock_parse_args.return_value = mock_args
        
        mock_structure = Mock()
        mock_structure.lattice.reciprocal_lattice_crystallographic.matrix = \
            np.array([[0.16, 0, 0], [0, 0.19, 0], [0, 0, 0.20]])
        mock_load_inputs.return_value = (mock_structure, {'ALGO': 'Normal'})
        mock_potcar = Mock()
        mock_load_potcar.return_value = mock_potcar
        
        mock_get_kspacing.return_value = (0.1, 0.2)
        mock_encut_targets = [Mock()]
        mock_kspacing_targets = [Mock()]
        mock_calc_encut.return_value = mock_encut_targets
        mock_calc_kspacing.return_value = mock_kspacing_targets
        
        main()
        
        # Should call print_dry_run_summary instead of execute
        mock_print_dry_run.assert_called_once_with(
            './test_dir',
            mock_encut_targets,
            mock_kspacing_targets,
            mock_potcar,
            None  # job_script
        )
        mock_execute.assert_not_called()


class InputValidationTestCase(unittest.TestCase):
    """Tests for input file validation."""
    
    def test_validates_poscar_exists_before_starting(self):
        """Test that POSCAR file existence is validated."""
        from vasppy.scripts.convergence_testing import validate_inputs
        
        with patch('os.path.isfile', return_value=False), \
             patch('os.path.exists', return_value=False):
            
            with self.assertRaises(FileNotFoundError):
                validate_inputs('missing.poscar', 'INCAR', './test_dir', job_script=None)
    
    def test_validates_incar_exists_before_starting(self):
        """Test that INCAR file existence is validated."""
        from vasppy.scripts.convergence_testing import validate_inputs
        
        with patch('os.path.isfile') as mock_isfile:
            mock_isfile.side_effect = lambda path: path == 'POSCAR'
            
            with patch('os.path.exists', return_value=False):
                with self.assertRaises(FileNotFoundError):
                    validate_inputs('POSCAR', 'missing.incar', './test_dir', job_script=None)
    
    def test_validates_output_directory_does_not_exist(self):
        """Test that output directory is checked for existence."""
        from vasppy.scripts.convergence_testing import validate_inputs
        
        with patch('os.path.isfile', return_value=True), \
             patch('os.path.exists', return_value=True):
            
            with self.assertRaises(FileExistsError):
                validate_inputs('POSCAR', 'INCAR', './existing_dir', job_script=None)
        
    def test_validate_inputs_checks_job_script_exists(self):
        """Test that job script file existence is validated when provided."""
        from vasppy.scripts.convergence_testing import validate_inputs
        
        with patch('os.path.isfile') as mock_isfile:
            mock_isfile.side_effect = lambda path: path in ['POSCAR', 'INCAR']
            
            with self.assertRaises(FileNotFoundError) as context:
                validate_inputs('POSCAR', 'INCAR', './new_dir', job_script='submit.sh')
            
            self.assertIn('submit.sh', str(context.exception))
            self.assertIn('job script', str(context.exception).lower())
    
    def test_validate_inputs_skips_job_script_check_when_none(self):
        """Test that job script validation is skipped when None."""
        from vasppy.scripts.convergence_testing import validate_inputs
        
        with patch('os.path.isfile', return_value=True), \
             patch('os.path.exists', return_value=False):
            # Should not raise - job_script=None means skip validation
            validate_inputs('POSCAR', 'INCAR', './new_dir', job_script=None)


class ErrorMessageQualityTestCase(unittest.TestCase):
    """Tests for error message quality and helpfulness."""

    @patch('vasppy.scripts.convergence_testing.os.path.isfile')
    @patch('vasppy.scripts.convergence_testing.parse_args')
    def test_poscar_error_message_is_helpful(self, mock_parse_args, mock_isfile):
        """Test that POSCAR error message suggests how to fix the problem."""
        
        mock_args = Mock(
            poscar='my_structure.poscar',
            incar='INCAR',
            directory='./test_dir'
        )
        mock_parse_args.return_value = mock_args
        mock_isfile.return_value = False
        
        with self.assertRaises(FileNotFoundError) as context:
            main()
        
        error_msg = str(context.exception)
        # Should mention the file that's missing
        self.assertIn('my_structure.poscar', error_msg)
        # Should be clear it's the POSCAR file
        self.assertIn('POSCAR', error_msg)

    @patch('vasppy.scripts.convergence_testing.os.path.exists')
    @patch('vasppy.scripts.convergence_testing.os.path.isfile')
    @patch('vasppy.scripts.convergence_testing.parse_args')
    def test_directory_exists_error_message_is_helpful(
        self, mock_parse_args, mock_isfile, mock_exists
    ):
        """Test that directory exists error suggests solutions."""
        
        mock_args = Mock(
            poscar='POSCAR',
            incar='INCAR',
            directory='./my_tests'
        )
        mock_parse_args.return_value = mock_args
        mock_isfile.return_value = True
        mock_exists.return_value = True
        
        with self.assertRaises(FileExistsError) as context:
            main()
        
        error_msg = str(context.exception)
        # Should mention the directory
        self.assertIn('my_tests', error_msg)
        # Should suggest a solution
        self.assertTrue(
            'different directory' in error_msg.lower() or 
            'choose another' in error_msg.lower() or
            'already exists' in error_msg.lower()
        )

        
class ConvergenceTargetTestCase(unittest.TestCase):
    """Tests for ConvergenceTarget dataclass."""
    
    def test_creates_target_with_encut_only(self):
        """Test creating a target that varies ENCUT."""
        
        
        target = ConvergenceTarget(
            path=Path('test_dir/ENCUT/100'),
            encut=100,
            kspacing=None,
            base_incar={'ALGO': 'Normal'}
        )
        
        self.assertEqual(target.encut, 100)
        self.assertIsNone(target.kspacing)
    
    def test_creates_target_with_kspacing_only(self):
        """Test creating a target that varies KSPACING."""
        
        
        target = ConvergenceTarget(
            path=Path('test_dir/KSPACING/0.3'),
            encut=None,
            kspacing=0.3,
            base_incar={'ALGO': 'Normal'}
        )
        
        self.assertIsNone(target.encut)
        self.assertEqual(target.kspacing, 0.3)
    
    def test_incar_params_includes_encut_when_set(self):
        """Test that incar_params includes ENCUT when specified."""
        
        target = ConvergenceTarget(
            path=Path('test_dir/ENCUT/100'),
            encut=100,
            kspacing=None,
            base_incar={'ALGO': 'Normal', 'EDIFF': 1e-6}
        )
        
        params = target.incar_params
        self.assertEqual(params['ENCUT'], 100)
        self.assertEqual(params['ALGO'], 'Normal')
        self.assertEqual(params['EDIFF'], 1e-6)
    
    def test_incar_params_includes_kspacing_when_set(self):
        """Test that incar_params includes KSPACING when specified."""
        
        target = ConvergenceTarget(
            path=Path('test_dir/KSPACING/0.3'),
            encut=None,
            kspacing=0.3,
            base_incar={'ALGO': 'Normal'}
        )
        
        params = target.incar_params
        self.assertEqual(params['KSPACING'], 0.3)
        self.assertEqual(params['ALGO'], 'Normal')
    
    def test_incar_params_includes_both_when_both_set(self):
        """Test that incar_params includes both ENCUT and KSPACING when specified."""
        
        target = ConvergenceTarget(
            path=Path('test_dir/matrix/100_0.3'),
            encut=100,
            kspacing=0.3,
            base_incar={'ALGO': 'Normal'}
        )
        
        params = target.incar_params
        self.assertEqual(params['ENCUT'], 100)
        self.assertEqual(params['KSPACING'], 0.3)
    
    def test_incar_params_does_not_modify_base_incar(self):
        """Test that accessing incar_params doesn't modify base_incar."""
        
        base = {'ALGO': 'Normal'}
        target = ConvergenceTarget(
            path=Path('test_dir/ENCUT/100'),
            encut=100,
            kspacing=None,
            base_incar=base
        )
        
        _ = target.incar_params
        
        # Original dict should be unchanged
        self.assertEqual(base, {'ALGO': 'Normal'})
        self.assertNotIn('ENCUT', base)
    
    def test_raises_error_when_neither_param_set(self):
        """Test that error is raised when both encut and kspacing are None."""
        
        with self.assertRaises(ValueError) as context:
            ConvergenceTarget(
                path=Path('test_dir'),
                encut=None,
                kspacing=None,
                base_incar={'ALGO': 'Normal'}
            )
        
        self.assertIn('at least one', str(context.exception).lower())
        
class CalculateEncutTargetsTestCase(unittest.TestCase):
    """Tests for calculate_encut_targets function."""

    def test_calculates_correct_number_of_targets(self):
        """Test that correct number of targets is returned."""
        from vasppy.scripts.convergence_testing import calculate_encut_targets
        
        encut_values = np.array([100, 200, 300])
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_encut_targets('./test_dir', encut_values, 0.3, incar_dict)
        
        self.assertEqual(len(targets), 3)

    def test_sets_correct_paths(self):
        """Test that paths are constructed correctly."""
        from vasppy.scripts.convergence_testing import calculate_encut_targets
        
        encut_values = np.array([100, 200])
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_encut_targets('./test_dir', encut_values, 0.3, incar_dict)
        
        self.assertEqual(str(targets[0].path), 'test_dir/ENCUT/100')
        self.assertEqual(str(targets[1].path), 'test_dir/ENCUT/200')

    def test_sets_encut_values_correctly(self):
        """Test that ENCUT values are set in targets."""
        from vasppy.scripts.convergence_testing import calculate_encut_targets
        
        encut_values = np.array([100, 200])
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_encut_targets('./test_dir', encut_values, 0.3, incar_dict)
        
        self.assertEqual(targets[0].encut, 100)
        self.assertEqual(targets[1].encut, 200)

    def test_sets_kspacing_to_none(self):
        """Test that KSPACING is None for ENCUT tests."""
        from vasppy.scripts.convergence_testing import calculate_encut_targets
        
        encut_values = np.array([100])
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_encut_targets('./test_dir', encut_values, 0.3, incar_dict)
        
        self.assertIsNone(targets[0].kspacing)

    def test_incar_params_includes_base_kspacing(self):
        """Test that base KSPACING is included in incar_params."""
        from vasppy.scripts.convergence_testing import calculate_encut_targets
        
        encut_values = np.array([100])
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_encut_targets('./test_dir', encut_values, 0.3, incar_dict)
        
        # Base KSPACING should be in base_incar, not as the kspacing attribute
        self.assertIn('KSPACING', targets[0].base_incar)
        self.assertEqual(targets[0].base_incar['KSPACING'], 0.3)


class CalculateKspacingTargetsTestCase(unittest.TestCase):
    """Tests for calculate_kspacing_targets function."""

    def test_calculates_correct_number_of_targets(self):
        """Test that correct number of targets is returned."""
        from vasppy.scripts.convergence_testing import calculate_kspacing_targets
        
        kspacing_values = (0.1, 0.2, 0.3)
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_kspacing_targets('./test_dir', kspacing_values, 400, incar_dict)
        
        self.assertEqual(len(targets), 3)

    def test_returns_convergence_target_objects(self):
        """Test that ConvergenceTarget objects are returned."""
        from vasppy.scripts.convergence_testing import calculate_kspacing_targets, ConvergenceTarget
        
        kspacing_values = (0.1, 0.2)
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_kspacing_targets('./test_dir', kspacing_values, 400, incar_dict)
        
        self.assertIsInstance(targets[0], ConvergenceTarget)
        self.assertIsInstance(targets[1], ConvergenceTarget)

    def test_sets_correct_paths(self):
        """Test that paths are constructed correctly."""
        from vasppy.scripts.convergence_testing import calculate_kspacing_targets
        
        kspacing_values = (0.1, 0.2)
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_kspacing_targets('./test_dir', kspacing_values, 400, incar_dict)
        
        self.assertEqual(str(targets[0].path), 'test_dir/KSPACING/0.1')
        self.assertEqual(str(targets[1].path), 'test_dir/KSPACING/0.2')

    def test_sets_kspacing_values_correctly(self):
        """Test that KSPACING values are set in targets."""
        from vasppy.scripts.convergence_testing import calculate_kspacing_targets
        
        kspacing_values = (0.1, 0.2)
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_kspacing_targets('./test_dir', kspacing_values, 400, incar_dict)
        
        self.assertEqual(targets[0].kspacing, 0.1)
        self.assertEqual(targets[1].kspacing, 0.2)

    def test_sets_encut_to_none(self):
        """Test that ENCUT is None for KSPACING tests."""
        from vasppy.scripts.convergence_testing import calculate_kspacing_targets
        
        kspacing_values = (0.1,)
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_kspacing_targets('./test_dir', kspacing_values, 400, incar_dict)
        
        self.assertIsNone(targets[0].encut)

    def test_incar_params_includes_base_encut(self):
        """Test that base ENCUT is included in incar_params."""
        from vasppy.scripts.convergence_testing import calculate_kspacing_targets
        
        kspacing_values = (0.1,)
        incar_dict = {'ALGO': 'Normal'}
        
        targets = calculate_kspacing_targets('./test_dir', kspacing_values, 400, incar_dict)
        
        # Base ENCUT should be in base_incar, not as the encut attribute
        self.assertIn('ENCUT', targets[0].base_incar)
        self.assertEqual(targets[0].base_incar['ENCUT'], 400)
        
class LoadInputsTestCase(unittest.TestCase):
    """Tests for load_inputs function."""

    @patch('vasppy.scripts.convergence_testing.Incar')
    @patch('vasppy.scripts.convergence_testing.Structure')
    def test_loads_structure_from_poscar(self, mock_structure_class, mock_incar_class):
        """Test that Structure is loaded from POSCAR file."""
        from vasppy.scripts.convergence_testing import load_inputs
        
        mock_structure = Mock()
        mock_structure_class.from_file.return_value = mock_structure
        mock_incar = Mock()
        mock_incar.as_dict.return_value = {'ALGO': 'Normal'}
        mock_incar_class.from_file.return_value = mock_incar
        
        structure, _ = load_inputs('POSCAR', 'INCAR')
        
        mock_structure_class.from_file.assert_called_once_with('POSCAR')
        self.assertEqual(structure, mock_structure)

    @patch('vasppy.scripts.convergence_testing.Incar')
    @patch('vasppy.scripts.convergence_testing.Structure')
    def test_loads_incar_from_file(self, mock_structure_class, mock_incar_class):
        """Test that INCAR is loaded from file."""
        from vasppy.scripts.convergence_testing import load_inputs
        
        mock_structure = Mock()
        mock_structure_class.from_file.return_value = mock_structure
        mock_incar = Mock()
        mock_incar.as_dict.return_value = {'ALGO': 'Normal', 'EDIFF': 1e-6}
        mock_incar_class.from_file.return_value = mock_incar
        
        _, incar_dict = load_inputs('POSCAR', 'INCAR')
        
        mock_incar_class.from_file.assert_called_once_with('INCAR')
        self.assertEqual(incar_dict, {'ALGO': 'Normal', 'EDIFF': 1e-6})


class PrintDryRunSummaryTestCase(unittest.TestCase):
    """Tests for print_dry_run_summary function."""

    @patch('sys.stdout', new_callable=StringIO)
    def test_prints_base_directory(self, mock_stdout):
        """Test that base directory is printed."""
        from vasppy.scripts.convergence_testing import print_dry_run_summary, ConvergenceTarget
        
        encut_targets = []
        kspacing_targets = []
        
        print_dry_run_summary('./test_dir', encut_targets, kspacing_targets, None)
        
        output = mock_stdout.getvalue()
        self.assertIn('test_dir', output)
        self.assertIn('base directory', output.lower())

    @patch('sys.stdout', new_callable=StringIO)
    def test_prints_encut_target_count(self, mock_stdout):
        """Test that ENCUT target count is printed."""
        from vasppy.scripts.convergence_testing import print_dry_run_summary, ConvergenceTarget
        
        encut_targets = [
            ConvergenceTarget(Path('test/ENCUT/100'), encut=100, kspacing=None, base_incar={}),
            ConvergenceTarget(Path('test/ENCUT/200'), encut=200, kspacing=None, base_incar={}),
        ]
        kspacing_targets = []
        
        print_dry_run_summary('./test_dir', encut_targets, kspacing_targets, None)
        
        output = mock_stdout.getvalue()
        self.assertIn('2', output)
        self.assertIn('ENCUT', output)

    @patch('sys.stdout', new_callable=StringIO)
    def test_prints_kspacing_target_count(self, mock_stdout):
        """Test that KSPACING target count is printed."""
        from vasppy.scripts.convergence_testing import print_dry_run_summary, ConvergenceTarget
        
        encut_targets = []
        kspacing_targets = [
            ConvergenceTarget(Path('test/KSPACING/0.1'), encut=None, kspacing=0.1, base_incar={}),
            ConvergenceTarget(Path('test/KSPACING/0.2'), encut=None, kspacing=0.2, base_incar={}),
            ConvergenceTarget(Path('test/KSPACING/0.3'), encut=None, kspacing=0.3, base_incar={}),
        ]
        
        print_dry_run_summary('./test_dir', encut_targets, kspacing_targets, None)
        
        output = mock_stdout.getvalue()
        self.assertIn('3', output)
        self.assertIn('KSPACING', output)

    @patch('sys.stdout', new_callable=StringIO)
    def test_indicates_dry_run_mode(self, mock_stdout):
        """Test that output indicates dry-run mode."""
        from vasppy.scripts.convergence_testing import print_dry_run_summary
        
        print_dry_run_summary('./test_dir', [], [], None)
        
        output = mock_stdout.getvalue()
        self.assertIn('DRY RUN', output)
        
    @patch('sys.stdout', new_callable=StringIO)
    def test_prints_potcar_info_when_provided(self, mock_stdout):
        """Test that POTCAR info is printed when potcar provided."""
        from vasppy.scripts.convergence_testing import print_dry_run_summary
        
        mock_potcar = Mock()
        
        print_dry_run_summary('./test_dir', [], [], mock_potcar)
        
        output = mock_stdout.getvalue()
        self.assertIn('POTCAR', output)
    
    @patch('sys.stdout', new_callable=StringIO)
    def test_does_not_print_potcar_when_none(self, mock_stdout):
        """Test that POTCAR info is not printed when None."""
        from vasppy.scripts.convergence_testing import print_dry_run_summary
        
        print_dry_run_summary('./test_dir', [], [], None)
        
        output = mock_stdout.getvalue()
        # Should complete without error, no POTCAR mention
        self.assertIn('DRY RUN', output)
        
    @patch('sys.stdout', new_callable=StringIO)
    def test_prints_job_script_info_when_provided(self, mock_stdout):
        """Test that job script info is printed when provided."""
        from vasppy.scripts.convergence_testing import print_dry_run_summary
        
        mock_potcar = Mock()
        
        print_dry_run_summary('./test_dir', [], [], mock_potcar, job_script='submit.sh')
        
        output = mock_stdout.getvalue()
        self.assertIn('submit.sh', output)
        self.assertIn('job script', output.lower())


class ExecuteConvergenceTestsTestCase(unittest.TestCase):
    """Tests for execute_convergence_tests function."""
    
    @patch('vasppy.scripts.convergence_testing.execute_targets')
    @patch('vasppy.scripts.convergence_testing.create_directory_structure')
    def test_creates_directory_structure(self, mock_create_dirs, mock_execute):
        """Test that directory structure is created."""
        from vasppy.scripts.convergence_testing import execute_convergence_tests
        
        mock_structure = Mock()
        
        execute_convergence_tests('./test_dir', [], [], mock_structure, None)
        
        mock_create_dirs.assert_called_once_with('./test_dir')
    
    @patch('vasppy.scripts.convergence_testing.execute_targets')
    @patch('vasppy.scripts.convergence_testing.create_directory_structure')
    def test_executes_encut_targets(self, mock_create_dirs, mock_execute):
        """Test that ENCUT targets are executed."""
        from vasppy.scripts.convergence_testing import execute_convergence_tests
        
        mock_structure = Mock()
        encut_targets = [Mock(), Mock()]
        kspacing_targets = []
        
        execute_convergence_tests('./test_dir', encut_targets, kspacing_targets, mock_structure, None)
        
        # First call should be for encut_targets
        first_call = mock_execute.call_args_list[0]
        self.assertEqual(first_call[0][0], encut_targets)
    
    @patch('vasppy.scripts.convergence_testing.execute_targets')
    @patch('vasppy.scripts.convergence_testing.create_directory_structure')
    def test_executes_kspacing_targets(self, mock_create_dirs, mock_execute):
        """Test that KSPACING targets are executed."""
        from vasppy.scripts.convergence_testing import execute_convergence_tests
        
        mock_structure = Mock()
        encut_targets = []
        kspacing_targets = [Mock(), Mock(), Mock()]
        
        execute_convergence_tests('./test_dir', encut_targets, kspacing_targets, mock_structure, None)
        
        # Second call should be for kspacing_targets
        second_call = mock_execute.call_args_list[1]
        self.assertEqual(second_call[0][0], kspacing_targets)
    
    @patch('vasppy.scripts.convergence_testing.execute_targets')
    @patch('vasppy.scripts.convergence_testing.create_directory_structure')
    def test_passes_potcar_to_execute_targets(self, mock_create_dirs, mock_execute):
        """Test that potcar is passed to execute_targets."""
        from vasppy.scripts.convergence_testing import execute_convergence_tests
        
        mock_structure = Mock()
        mock_potcar = Mock()
        
        execute_convergence_tests('./test_dir', [], [], mock_structure, mock_potcar)
        
        # Both execute_targets calls should receive the potcar
        self.assertEqual(mock_execute.call_args_list[0][0][2], mock_potcar)
        self.assertEqual(mock_execute.call_args_list[1][0][2], mock_potcar)
    
    @patch('vasppy.scripts.convergence_testing.execute_targets')
    @patch('vasppy.scripts.convergence_testing.create_directory_structure')
    def test_passes_none_when_no_potcar(self, mock_create_dirs, mock_execute):
        """Test that None is passed when no potcar provided."""
        from vasppy.scripts.convergence_testing import execute_convergence_tests
        
        mock_structure = Mock()
        
        execute_convergence_tests('./test_dir', [], [], mock_structure, None)
        
        # Both execute_targets calls should receive None for potcar
        self.assertIsNone(mock_execute.call_args_list[0][0][2])
        self.assertIsNone(mock_execute.call_args_list[1][0][2])
        
class LoadPotcarTestCase(unittest.TestCase):
    """Tests for load_potcar function."""
    
    @patch('vasppy.scripts.convergence_testing.Potcar')
    def test_creates_potcar_from_pseudopotentials(self, mock_potcar_class):
        """Test that Potcar is created from pseudopotentials list."""
        from vasppy.scripts.convergence_testing import load_potcar
        
        mock_potcar = Mock()
        mock_potcar_class.return_value = mock_potcar
        
        result = load_potcar(['Li_sv', 'O'], None)
        
        mock_potcar_class.assert_called_once_with(['Li_sv', 'O'])
        self.assertEqual(result, mock_potcar)
    
    @patch('vasppy.scripts.convergence_testing.Potcar')
    def test_loads_potcar_from_file(self, mock_potcar_class):
        """Test that Potcar is loaded from file."""
        from vasppy.scripts.convergence_testing import load_potcar
        
        mock_potcar = Mock()
        mock_potcar_class.from_file.return_value = mock_potcar
        
        result = load_potcar(None, 'POTCAR')
        
        mock_potcar_class.from_file.assert_called_once_with('POTCAR')
        self.assertEqual(result, mock_potcar)
    
    @patch('vasppy.scripts.convergence_testing.Potcar')
    def test_raises_error_when_no_potcar_specified(self, mock_potcar_class):
        """Test that ValueError is raised when no potcar specified."""
        from vasppy.scripts.convergence_testing import load_potcar
        
        with self.assertRaises(ValueError) as context:
            load_potcar(None, None)
        
        self.assertIn('must be provided', str(context.exception).lower())
        mock_potcar_class.assert_not_called()
        mock_potcar_class.from_file.assert_not_called()
    
    @patch('vasppy.scripts.convergence_testing.Potcar')
    def test_prefers_pseudopotentials_over_file(self, mock_potcar_class):
        """Test that pseudopotentials takes precedence if both provided."""
        from vasppy.scripts.convergence_testing import load_potcar
        
        mock_potcar = Mock()
        mock_potcar_class.return_value = mock_potcar
        
        result = load_potcar(['Li_sv'], 'POTCAR')
        
        # Should use pseudopotentials, not file
        mock_potcar_class.assert_called_once_with(['Li_sv'])
        mock_potcar_class.from_file.assert_not_called()


if __name__ == '__main__':
    unittest.main()
