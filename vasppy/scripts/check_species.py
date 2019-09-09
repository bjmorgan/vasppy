#! /usr/bin/env python3

from vasppy.poscar import Poscar
from vasppy.summary import potcar_spec, potcar_sets
import argparse
import re

"""
A command line utility for testing species consistency between a VASP POSCAR and POTCAR pair of files. Species are considered consistent if the species labels used in the POSCAR file match the start of the pseudopotential labels in the POTCAR file, in order. e.g. a POSCAR that contains `Ti O` will match a POTCAR that contains `Ti_pv O`. If any species labels do not match the script raises an AttributeError.

The `-p` flag will check that all the pseudopotentials in the POTCAR file belong to a specific pseudopotential set.
"""

def parse_command_line_arguments():
    parser = argparse.ArgumentParser( description='Check species consistency between a VASP POSCAR file and a POTCAR file.' )
    parser.add_argument( 'poscar', help="filename of the VASP POSCAR to be processed", nargs='?', default='POSCAR' )
    parser.add_argument( 'potcar', help="filename of the VASP POTCAR to be processed", nargs='?', default='POTCAR' )
    parser.add_argument( '-p', '--ppset', help="check whether the POTCAR pseudopotentials belong to a specific pseudopotential set", choices=potcar_sets )
    return parser.parse_args()

def main():
    args = parse_command_line_arguments()
    poscar = Poscar.from_file( args.poscar )
    potcars = potcar_spec( args.potcar )
    for i, ( species, potcar ) in enumerate( zip( poscar.atoms, potcars ), 1 ):
        matching_potcar = potcar.startswith( species )
        if not matching_potcar:
            raise AttributeError( 'Species {} mismatch:\nPOSCAR contains {}\nPOTCAR contains {}'.format( i, species, potcar ) )
        if args.ppset:
            this_ppset = potcars[potcar]
            if args.ppset != this_ppset:
                raise AttributeError( 'Pseudopotential set mismatch: {}'.format( potcars.values() ) )

if __name__ == '__main__':
    main()
