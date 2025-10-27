"""Generate a series of VASP inputs for convergence testing."""
import argparse
import os
import numpy as np
import shutil
from pymatgen.core import Structure
from pymatgen.io.vasp import Incar, Potcar
from vasppy.kpoints import get_convergence_testing_kspacing
from pathlib import Path
from dataclasses import dataclass
from typing import Any


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.
    
    Returns:
        Parsed command line arguments.
        
    """
    parser = argparse.ArgumentParser(
        description='Generate a series of VASP inputs for convergence testing. '
                    'Requires a template INCAR, the geometry of the system (a POSCAR), '
                    'and pseudopotentials (via --pseudopotentials or --potcar-file).'
    )
    parser.add_argument('-i', '--incar', required=True, help='Specify the template INCAR.')
    parser.add_argument('-p', '--poscar', required=True, help='Specify the geometry of the system (POSCAR).')
    potcar_group = parser.add_mutually_exclusive_group(required=True)
    potcar_group.add_argument(
        '--pseudopotentials', nargs='+', default=None,
        help='Specify the pseudopotentials (POTCARs) to be used e.g. Li_sv Mn_pv O. '
             'Mutually exclusive with --potcar-file. '
             'Defaults to None i.e. no POTCARs will be written. '
             'This functionality requires Pymatgen to be set up for VASP POTCARs.'
    )
    potcar_group.add_argument(
        '--potcar-file', default=None,
        help='Path to existing POTCAR file to use for all calculations. '
             'Mutually exclusive with --pseudopotentials.'
    )
    parser.add_argument(
        '-e', '--encut', type=int, nargs=3, default=(100, 700, 50),
        help='Set the upper/lower bounds and step size for the basis set size (ENCUT). Defaults to 100 700 50.'
    )
    parser.add_argument(
        '-k', '--kspacing', type=float, nargs=3, default=(0.1, 0.8, 0.02),
        help='Set the upper/lower bounds and step size for the minimum allowed distance between k-points (KSPACING). '
             'Defaults to 0.1 0.8 0.02'
    )
    parser.add_argument(
        '-d', '--directory', default='./convergence_testing',
        help='Specify the directory in which to place the generated VASP inputs. Defaults to ./convergence_testing.'
    )
    parser.add_argument(
        '--base-encut', type=int, default=400,
        help='Set the value of ENCUT for the KSPACING convergence tests. Defaults to 400.'
    )
    parser.add_argument(
        '--base-kspacing', type=float, default=0.3,
        help='Set the value of KSPACING for the ENCUT convergence tests. Defaults to 0.3.'
    )
    parser.add_argument(
        '--job-script', default=None,
        help='Path to job submission script to copy into each convergence test directory.'
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Show what would be created without actually creating files or directories.'
    )
    args = parser.parse_args()
    args.encut = tuple(args.encut)
    args.kspacing = tuple(args.kspacing)
    return args

@dataclass
class ConvergenceTarget:
    """Represents a single convergence test calculation.
    
    Attributes:
        path: Directory path for this calculation.
        encut: ENCUT value to test (None if not varying).
        kspacing: KSPACING value to test (None if not varying).
        base_incar: Base INCAR parameters shared across tests.
        
    """
    path: Path
    encut: int | None
    kspacing: float | None
    base_incar: dict[str, Any]
    
    # Parameters that should be included in INCAR
    CONVERGENCE_PARAMS = ['encut', 'kspacing']
    
    def __post_init__(self):
        """Validate that at least one parameter is set."""
        if self.encut is None and self.kspacing is None:
            raise ValueError(
                "At least one of encut or kspacing must be set for a convergence test"
            )
    
    @property
    def incar_params(self) -> dict[str, Any]:
        """Get complete INCAR parameters for writing.
        
        Returns:
            Dictionary of INCAR parameters including base params and test values.
            
        """
        params = self.base_incar.copy()
        
        for param_name in self.CONVERGENCE_PARAMS:
            value = getattr(self, param_name)
            if value is not None:
                params[param_name.upper()] = value
        
        return params
            
def create_directory_structure(base_dir: str | Path) -> None:
    """Create the base directory structure for convergence testing.
    
    Args:
        base_dir: Base directory path for convergence tests.
        
    """
    base_path = Path(base_dir)
    os.mkdir(str(base_path))
    os.mkdir(str(base_path / 'ENCUT'))
    os.mkdir(str(base_path / 'KSPACING'))
    
    
def write_vasp_input_files(
    directory: str | Path,
    incar_params: dict,
    structure: Structure,
    potcar: Potcar,
    job_script: str | None = None
) -> None:
    """Write VASP input files to a directory.
    
    Args:
        directory: Directory path to write files to.
        incar_params: Dictionary of INCAR parameters.
        structure: Pymatgen Structure object.
        potcar: Potcar object.
        job_script: Optional path to job script file to copy.
        
    """
    dir_path = Path(directory)
    incar = Incar.from_dict(incar_params)
    incar.write_file(str(dir_path / 'INCAR'))
    structure.to(fmt='poscar', filename=str(dir_path / 'POSCAR'))
    potcar.write_file(str(dir_path / 'POTCAR'))
    
    if job_script is not None:
        shutil.copy2(job_script, str(dir_path / Path(job_script).name))
    
def calculate_encut_targets(
    base_dir: str | Path,
    encut_values: np.ndarray,
    base_kspacing: float,
    incar_dict: dict
) -> list[ConvergenceTarget]:
    """Calculate convergence test targets for ENCUT.
    
    Args:
        base_dir: Base directory for convergence tests.
        encut_values: Array of ENCUT values to test.
        base_kspacing: Fixed KSPACING value to use.
        incar_dict: Base INCAR parameters dictionary.
        
    Returns:
        List of ConvergenceTarget objects for ENCUT tests.
        
    """
    base_path = Path(base_dir)
    targets = []
    
    # Add base KSPACING to base_incar for these tests
    base_incar_with_kspacing = incar_dict.copy()
    base_incar_with_kspacing['KSPACING'] = base_kspacing
    
    for energy_cutoff in encut_values:
        path = base_path / 'ENCUT' / str(energy_cutoff)
        target = ConvergenceTarget(
            path=path,
            encut=int(energy_cutoff),
            kspacing=None,
            base_incar=base_incar_with_kspacing
        )
        targets.append(target)
    
    return targets


def calculate_kspacing_targets(
    base_dir: str | Path,
    kspacing_values: tuple[float, ...],
    base_encut: int,
    incar_dict: dict
) -> list[ConvergenceTarget]:
    """Calculate convergence test targets for KSPACING.
    
    Args:
        base_dir: Base directory for convergence tests.
        kspacing_values: Tuple of KSPACING values to test.
        base_encut: Fixed ENCUT value to use.
        incar_dict: Base INCAR parameters dictionary.
        
    Returns:
        List of ConvergenceTarget objects for KSPACING tests.
        
    """
    base_path = Path(base_dir)
    targets = []
    
    # Add base ENCUT to base_incar for these tests
    base_incar_with_encut = incar_dict.copy()
    base_incar_with_encut['ENCUT'] = base_encut
    
    for minimum_distance in kspacing_values:
        path = base_path / 'KSPACING' / str(minimum_distance)
        target = ConvergenceTarget(
            path=path,
            encut=None,
            kspacing=minimum_distance,
            base_incar=base_incar_with_encut
        )
        targets.append(target)
    
    return targets 

def execute_targets(
    targets: list[ConvergenceTarget],
    structure: Structure,
    potcar: Potcar,
    job_script: str | None = None
) -> None:
    """Execute convergence test targets by creating directories and writing files.
    
    Args:
        targets: List of ConvergenceTarget objects to execute.
        structure: Pymatgen Structure object.
        potcar: Potcar object.
        job_script: Optional path to job script file to copy.
        
    """
    for target in targets:
        os.mkdir(str(target.path))
        write_vasp_input_files(target.path, target.incar_params, structure, potcar, job_script)

def validate_inputs(
    poscar_path: str,
    incar_path: str,
    output_dir: str,
    job_script: str | None = None
) -> None:
    """Validate input files and output directory before starting work.
    
    Args:
        poscar_path: Path to POSCAR file.
        incar_path: Path to INCAR file.
        output_dir: Path to output directory.
        job_script: Optional path to job script file.
        
    Raises:
        FileNotFoundError: If POSCAR, INCAR, or job script file doesn't exist.
        FileExistsError: If output directory already exists.
        
    """
    # Check POSCAR exists
    if not os.path.isfile(poscar_path):
        raise FileNotFoundError(
            f"POSCAR file not found: {poscar_path}\n"
            f"Please provide a valid POSCAR file using the -p/--poscar argument."
        )
    
    # Check INCAR exists
    if not os.path.isfile(incar_path):
        raise FileNotFoundError(
            f"INCAR file not found: {incar_path}\n"
            f"Please provide a valid INCAR file using the -i/--incar argument."
        )
    
    # Check job script exists if provided
    if job_script is not None and not os.path.isfile(job_script):
        raise FileNotFoundError(
            f"Job script file not found: {job_script}\n"
            f"Please provide a valid job script file using the --job-script argument."
        )
    
    # Check output directory doesn't exist
    if os.path.exists(output_dir):
        raise FileExistsError(
            f"Output directory already exists: {output_dir}\n"
            f"Please choose a different directory using the -d/--directory argument."
        )
            
def load_inputs(poscar_path: str, incar_path: str) -> tuple[Structure, dict]:
    """Load structure and INCAR template.
    
    Args:
        poscar_path: Path to POSCAR file.
        incar_path: Path to INCAR file.
        
    Returns:
        Tuple of (Structure, INCAR dict).
        
    """
    print(f"Loading structure from {poscar_path}...")
    structure = Structure.from_file(poscar_path)
    print(f"Loading INCAR template from {incar_path}...")
    base_incar_dict = Incar.from_file(incar_path).as_dict()
    return structure, base_incar_dict

def print_dry_run_summary(
    directory: str,
    encut_targets: list[ConvergenceTarget],
    kspacing_targets: list[ConvergenceTarget],
    potcar: Potcar,
    job_script: str | None = None
) -> None:
    """Print summary of what would be created in dry-run mode.
    
    Args:
        directory: Base directory path.
        encut_targets: List of ENCUT test targets.
        kspacing_targets: List of KSPACING test targets.
        potcar: Potcar object.
        job_script: Optional path to job script file.
        
    """
    print("\n=== DRY RUN - No files will be created ===\n")
    print(f"Would create base directory: {directory}\n")
    print(f"Would create {len(encut_targets)} ENCUT test directories:")
    for target in encut_targets:
        print(f"  - {target.path} (ENCUT={target.encut})")
    print(f"\nWould create {len(kspacing_targets)} KSPACING test directories:")
    for target in kspacing_targets:
        print(f"  - {target.path} (KSPACING={target.kspacing})")
    print(f"\nTotal: {len(encut_targets) + len(kspacing_targets)} convergence tests")
    print(f"POTCARs will be included in all calculations")
    if job_script:
        print(f"Job script '{job_script}' will be copied to each directory")
                    
def execute_convergence_tests(
    directory: str,
    encut_targets: list[ConvergenceTarget],
    kspacing_targets: list[ConvergenceTarget],
    structure: Structure,
    potcar: Potcar,
    job_script: str | None = None
) -> None:
    """Execute convergence tests by creating directories and writing files.
    
    Args:
        directory: Base directory path.
        encut_targets: List of ENCUT test targets.
        kspacing_targets: List of KSPACING test targets.
        structure: Pymatgen Structure object.
        potcar: Potcar object.
        job_script: Optional path to job script file to copy.
        
    """
    print(f"\nCreating directory structure at {directory}...")
    create_directory_structure(directory)
    
    print(f"Generating {len(encut_targets)} ENCUT convergence tests...")
    execute_targets(encut_targets, structure, potcar, job_script)
    
    print(f"Generating {len(kspacing_targets)} KSPACING convergence tests...")
    execute_targets(kspacing_targets, structure, potcar, job_script)
    
    print(f"\n✓ Complete! Created {len(encut_targets) + len(kspacing_targets)} tests in {directory}")

def load_potcar(pseudopotentials: list[str] | None, potcar_file: str | None) -> Potcar:
    """Load or create Potcar object.
    
    Args:
        pseudopotentials: List of pseudopotential names to create Potcar from.
        potcar_file: Path to existing POTCAR file to load.
        
    Returns:
        Potcar object.
        
    Raises:
        ValueError: If neither pseudopotentials nor potcar_file is provided.
        
    """
    if pseudopotentials:
        potcar = Potcar(pseudopotentials)
    elif potcar_file:
        potcar = Potcar.from_file(potcar_file)
    else:
        raise ValueError("Either pseudopotentials or potcar_file must be provided")
    assert(isinstance(potcar, Potcar))
    return potcar
                    
def main() -> None:
    """Main entry point for convergence testing script."""
    args = parse_args()
    
    # Validate and load inputs
    print("Validating inputs...")
    validate_inputs(args.poscar, args.incar, args.directory, job_script=args.job_script)
    structure, base_incar_dict = load_inputs(args.poscar, args.incar)
    
    # Load potcar
    potcar = load_potcar(args.pseudopotentials, args.potcar_file)
    
    # Calculate parameter ranges
    print("Calculating convergence test parameters...")
    reciprocal_lattice_vectors = structure.lattice.reciprocal_lattice_crystallographic.matrix
    kspacing_min, kspacing_max, step = args.kspacing
    kspacing_values = get_convergence_testing_kspacing(
        reciprocal_lattice_vectors, 
        (kspacing_min, kspacing_max), 
        step
    )
    encut_min, encut_max, step = args.encut
    encut_values = np.arange(encut_min, encut_max + step, step)
    
    # Calculate targets
    encut_targets = calculate_encut_targets(
        args.directory, encut_values, args.base_kspacing, base_incar_dict
    )
    kspacing_targets = calculate_kspacing_targets(
        args.directory, kspacing_values, args.base_encut, base_incar_dict
    )
    
    # Execute or dry-run
    if args.dry_run:
        print_dry_run_summary(args.directory, encut_targets, kspacing_targets, potcar, args.job_script)
    else:
        execute_convergence_tests(args.directory, encut_targets, kspacing_targets, structure, potcar, args.job_script)


if __name__ == '__main__':
    main()