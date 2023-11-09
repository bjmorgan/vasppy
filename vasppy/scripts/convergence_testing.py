import argparse
import os
import numpy as np
from pymatgen.core import Structure # type: ignore
from pymatgen.io.vasp import Incar, Potcar # type: ignore
from vasppy.kpoints import get_convergence_testing_kspacing

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Generate a series of VASP inputs for convergence testing. At minimum, this requires a template INCAR and the geometry of the system (a POSCAR).')
    parser.add_argument('-i', '--incar', required=True, help='Specify the template INCAR.')
    parser.add_argument('-p', '--poscar', required=True, help='Specify the geometry of the system (POSCAR).')
    parser.add_argument('-e', '--encut', nargs=3, default=(100, 700, 50), help='Set the upper/lower bounds and step size for the basis set size (ENCUT). Defaults to 100 700 50.')
    parser.add_argument('-k', '--kspacing', nargs=3, default=(0.1, 0.8, 0.02), help='Set the upper/lower bounds and step size for the minimum allowed distance between k-points (KSPACING). Defaults to 0.1 0.8 0.02')
    parser.add_argument('-d', '--directory', nargs='?', default='./convergence_testing', help='Specify the directory in which to place the generated VASP inputs. Defaults to ./convergence_testing.')
    parser.add_argument('--base-encut', nargs='?', default=400, help='Set the value of ENCUT for the KSPACING convergence tests. Defaults to 400.')
    parser.add_argument('--base-kspacing', nargs='?', default=0.3, help='Set the value of KSPACING for the ENCUT convergence tests. Defaults to 0.3.')
    parser.add_argument('--pseudopotentials', nargs='+', default=None, help='Specify the pseudopotentials (POTCARs) to be used e.g. Li_sv Mn_pv O. Defaults to None i.e. no POTCARs will be written. This functionality requires Pymatgen to be set up for VASP POTCARs.')
    args = parser.parse_args()

    return args

def main() -> None:
    args = parse_args()

    structure = Structure.from_file(args.poscar)
    reciprocal_lattice_vectors = structure.lattice.reciprocal_lattice_crystallographic.matrix
    kspacing_min, kspacing_max, step = args.kspacing
    kspacing = get_convergence_testing_kspacing(reciprocal_lattice_vectors, (kspacing_min, kspacing_max), step)

    encut_min, encut_max, step = args.encut
    encut = np.arange(encut_min, encut_max + step, step)

    base_incar_params = Incar.from_file(args.incar).as_dict()

    os.mkdir(args.directory)
    os.mkdir(args.directory + '/ENCUT')
    os.mkdir(args.directory + '/KSPACING')

    if args.pseudopotentials is not None:
        potcar = Potcar(args.pseudopotentials)

    for energy_cutoff in encut:
        path = args.directory + '/ENCUT' + f'/{energy_cutoff}'
        os.mkdir(path)

        incar_params = base_incar_params.copy()
        incar_params['ENCUT'] = energy_cutoff
        incar_params['KSPACING'] = args.base_kspacing
        incar = Incar.from_dict(incar_params)

        incar.write_file(path + '/INCAR')
        structure.to(fmt='poscar', filename=path + '/POSCAR')
        if args.pseudopotentials is not None:
            potcar.write_file(path + '/POTCAR')

    for minimum_distance in kspacing:
        path = args.directory + '/KSPACING' + f'/{minimum_distance}'
        os.mkdir(path)

        incar_params = base_incar_params.copy()
        incar_params['KSPACING'] = minimum_distance
        incar_params['ENCUT'] = args.base_encut
        incar = Incar.from_dict(incar_params)

        incar.write_file(path + '/INCAR')
        structure.to(fmt='poscar', filename=path + '/POSCAR')
        if args.pseudopotentials is not None:
            potcar.write_file(path + '/POTCAR')

if __name__ == '__main__':
    main()
