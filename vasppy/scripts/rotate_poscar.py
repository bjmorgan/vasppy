#! /usr/bin/env python3

from vasppy.poscar import Poscar
import math
import argparse

def parse_command_line_arguments():
    parser = argparse.ArgumentParser( description='Rotates the cell lattice in VASP POSCAR files' )
    parser.add_argument( 'poscar', help="filename of the VASP POSCAR to be processed" )
    parser.add_argument( '-a', '--axis', nargs=3, type=float, help="vector for rotation axis", required=True )
    parser.add_argument( '-d', '--degrees', type=int, help="rotation angle in degrees", required=True )
    args = parser.parse_args()
    return( args )

def main():
    args = parse_command_line_arguments()
    poscar = Poscar()
    poscar.read_from( args.poscar )
    theta = math.pi * args.degrees / 180.0
    poscar.cell.rotate( args.axis, theta )
    poscar.output()

if __name__ == "__main__":
    main()
