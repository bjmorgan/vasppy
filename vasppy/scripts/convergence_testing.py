"""Generate a series of VASP inputs for convergence testing."""
import argparse
import os
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Incar, Potcar
from vasppy.kpoints import get_convergence_testing_kspacing


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.
    
    Returns:
        Parsed command line arguments.
        
    """
    parser = argparse.ArgumentParser(
        description='Generate a series of VASP inputs for convergence testing. '
                    'At minimum, this requires a template INCAR and the geometry of the system (a POSCAR).'
    )
    parser.add_argument('-i', '--incar', required=True, help='Specify the template INCAR.')
    parser.add_argument('-p', '--poscar', required=True, help='Specify the geometry of the system (POSCAR).')
    parser.add_argument(
        '-e', '--encut', nargs=3, default=(100, 700, 50),
        help='Set the upper/lower bounds and step size for the basis set size (ENCUT). Defaults to 100 700 50.'
    )
    parser.add_argument(
        '-k', '--kspacing', nargs=3, default=(0.1, 0.8, 0.02),
        help='Set the upper/lower bounds and step size for the minimum allowed distance between k-points (KSPACING). '
             'Defaults to 0.1 0.8 0.02'
    )
    parser.add_argument(
        '-d', '--directory', nargs='?', default='./convergence_testing',
        help='Specify the directory in which to place the generated VASP inputs. Defaults to ./convergence_testing.'
    )
    parser.add_argument(
        '--base-encut', nargs='?', default=400,
        help='Set the value of ENCUT for the KSPACING convergence tests. Defaults to 400.'
    )
    parser.add_argument(
        '--base-kspacing', nargs='?', default=0.3,
        help='Set the value of KSPACING for the ENCUT convergence tests. Defaults to 0.3.'
    )
    parser.add_argument(
        '--pseudopotentials', nargs='+', default=None,
        help='Specify the pseudopotentials (POTCARs) to be used e.g. Li_sv Mn_pv O. '
             'Defaults to None i.e. no POTCARs will be written. '
             'This functionality requires Pymatgen to be set up for VASP POTCARs.'
    )
    args = parser.parse_args()
    return args


def create_directory_structure(base_dir: str) -> None:
    """Create the base directory structure for convergence testing.
    
    Args:
        base_dir: Base directory path for convergence tests.
        
    """
    os.mkdir(base_dir)
    os.mkdir(base_dir + '/ENCUT')
    os.mkdir(base_dir + '/KSPACING')


def write_vasp_input_files(
    directory: str,
    incar_params: dict,
    structure: Structure,
    potcar: Potcar | None = None
) -> None:
    """Write VASP input files to a directory.
    
    Args:
        directory: Directory path to write files to.
        incar_params: Dictionary of INCAR parameters.
        structure: Pymatgen Structure object.
        potcar: Optional Potcar object. If provided, POTCAR will be written.
        
    """
    incar = Incar.from_dict(incar_params)
    incar.write_file(directory + '/INCAR')
    structure.to(fmt='poscar', filename=directory + '/POSCAR')
    if potcar is not None:
        potcar.write_file(directory + '/POTCAR')


def generate_encut_inputs(
    base_dir: str,
    encut_values: np.ndarray,
    base_kspacing: float,
    incar_dict: dict,
    structure: Structure,
    potcar: Potcar | None = None
) -> None:
    """Generate VASP input files for ENCUT convergence testing.
    
    Args:
        base_dir: Base directory for convergence tests.
        encut_values: Array of ENCUT values to test.
        base_kspacing: Fixed KSPACING value to use for ENCUT tests.
        incar_dict: Base INCAR parameters dictionary.
        structure: Pymatgen Structure object.
        potcar: Optional Potcar object.
        
    """
    for energy_cutoff in encut_values:
        path = base_dir + '/ENCUT' + f'/{energy_cutoff}'
        os.mkdir(path)
        
        incar_params = incar_dict.copy()
        incar_params['ENCUT'] = energy_cutoff
        incar_params['KSPACING'] = base_kspacing
        
        write_vasp_input_files(path, incar_params, structure, potcar)


def generate_kspacing_inputs(
    base_dir: str,
    kspacing_values: tuple[float, ...],
    base_encut: float,
    incar_dict: dict,
    structure: Structure,
    potcar: Potcar | None = None
) -> None:
    """Generate VASP input files for KSPACING convergence testing.
    
    Args:
        base_dir: Base directory for convergence tests.
        kspacing_values: Tuple of KSPACING values to test.
        base_encut: Fixed ENCUT value to use for KSPACING tests.
        incar_dict: Base INCAR parameters dictionary.
        structure: Pymatgen Structure object.
        potcar: Optional Potcar object.
        
    """
    for minimum_distance in kspacing_values:
        path = base_dir + '/KSPACING' + f'/{minimum_distance}'
        os.mkdir(path)
        
        incar_params = incar_dict.copy()
        incar_params['KSPACING'] = minimum_distance
        incar_params['ENCUT'] = base_encut
        
        write_vasp_input_files(path, incar_params, structure, potcar)


def main() -> None:
    """Main entry point for convergence testing script."""
    args = parse_args()

    # Load structure and template INCAR
    structure = Structure.from_file(args.poscar)
    base_incar_dict = Incar.from_file(args.incar).as_dict()
    
    # Calculate parameter ranges
    reciprocal_lattice_vectors = structure.lattice.reciprocal_lattice_crystallographic.matrix
    kspacing_min, kspacing_max, step = args.kspacing
    kspacing_values = get_convergence_testing_kspacing(
        reciprocal_lattice_vectors, 
        (kspacing_min, kspacing_max), 
        step
    )
    
    encut_min, encut_max, step = args.encut
    encut_values = np.arange(encut_min, encut_max + step, step)
    
    # Create POTCAR if pseudopotentials specified
    potcar = Potcar(args.pseudopotentials) if args.pseudopotentials is not None else None
    
    # Generate all convergence test inputs
    create_directory_structure(args.directory)
    generate_encut_inputs(
        args.directory, 
        encut_values, 
        args.base_kspacing, 
        base_incar_dict, 
        structure, 
        potcar
    )
    generate_kspacing_inputs(
        args.directory, 
        kspacing_values, 
        args.base_encut, 
        base_incar_dict, 
        structure, 
        potcar
    )


if __name__ == '__main__':
    main()