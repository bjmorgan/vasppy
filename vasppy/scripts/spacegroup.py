#! /usr/bin/env python3 

from vasppy.poscar import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import argparse

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser( description='Finds the spacegroup for a VASP POSCAR file' )
    parser.add_argument( 'poscar', help="filename of the VASP POSCAR to be processed" )
    parser.add_argument( '-s', '--symprec', type=float, help="Precision for symmetry analuysis (defalut=1e-3)", default=1e-3 )
    args = parser.parse_args()
    return( args )

def main():
    args = parse_command_line_arguments()
    # initialise
    poscar = Poscar() # this doesn't really need vasppy. Could just use pymatgen to read the POSCAR
    # read POSCAR file
    poscar.read_from( args.poscar )
    structure = poscar.to_pymatgen_structure()
    symmetry_analyzer = SpacegroupAnalyzer( structure, symprec = args.symprec )
    print( symmetry_analyzer.get_space_group_symbol() )

if __name__ == "__main__":
    main()
