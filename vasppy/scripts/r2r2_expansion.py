#! /usr/bin/env python3

from vasppy.poscar import Poscar
import argparse
import numpy as np
from typing import Literal

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser(description="Generate a sqrt(2) x sqrt(2) supercell from a VASP POSCAR")
    parser.add_argument("poscar", help="filename of the VASP POSCAR to be processed")
    parser.add_argument(
        "-a",
        "--axis",
        choices=['x', 'y', 'z'],
        type=str,
        help="normal vector for sqrt(2) x sqrt(2) expansion",
        required=True,
    )
    args = parser.parse_args()
    return args

def sqrt2_by_sqrt2_expansion(
        poscar: Poscar,
        axis: Literal['x', 'y', 'z']
        ) -> Poscar:
        axis_vectors = {'x': [1, 0, 0],
                        'y': [0, 1, 0],
                        'z': [0, 0, 1]}
        poscar.cell.rotate(axis=axis_vectors[axis], theta=np.pi/4)
        poscar = poscar.replicate(2, 2, 1)
        cc = poscar.cartesian_coordinates()
        poscar.cell.matrix = np.diag(poscar.cell.matrix.diagonal())
        poscar.coordinates = cc.dot(np.linalg.inv(poscar.cell.matrix))
        return poscar

def main():
    args = parse_command_line_arguments()
    poscar = Poscar.from_file(args.poscar)
    sqrt2_by_sqrt2_expansion(
        poscar = poscar,
        axis = args.axis
    ).output()

if __name__ == '__main__':
    main()


