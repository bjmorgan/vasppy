#! /usr/bin/env python3 

from vasppy.poscar import Poscar
import argparse

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser( description='Converts a VASP POSCAR file to the .xtl file format' )
    parser.add_argument( 'poscar', help="filename of the VASP POSCAR to be processed" )
    parser.add_argument( '-s', '--symprec', type=float, help="Symmetry precision for a symmetrised .cif output" )
    args = parser.parse_args()
    return( args )

def main():
    args = parse_command_line_arguments()
    # initialise
    poscar = Poscar()
    # read POSCAR file
    poscar.read_from( args.poscar )
    poscar.output_as_cif( args.symprec )

if __name__ == "__main__":
    main()
