#! /usr/bin/env python3 

import poscar
import sys
import re
import string
import numpy as np
import argparse
import copy
import signal      

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser( description='Manipulates VASP POSCAR files' )
    parser.add_argument( 'poscar', help="filename of the VASP POSCAR to be processed" )
    parser.add_argument( '-l', '--label', type=int, choices=[ 1, 4 ], help="label coordinates with atom name at position {1,4}" )
    parser.add_argument( '-c', '--coordinates-only', help='only output coordinates', action='store_true' )
    parser.add_argument( '-t', '--coordinate-type', type=str, choices=[ 'c', 'cartesian', 'd', 'direct' ], default='direct', help="specify coordinate type for output {(c)artesian|(d)irect} [default = (d)irect]" )
    parser.add_argument( '-g', '--group', help='group atoms within supercell', action='store_true' )
    parser.add_argument( '-s', '--supercell', type=int, nargs=3, metavar=( 'h', 'k', 'l' ), help='construct supercell by replicating (h,k,l) times along [a b c]' )
    parser.add_argument( '-b', '--bohr', action='store_true', help='assumes the input file is in Angstrom, and converts everything to bohr')
    args = parser.parse_args()
    return( args )

if __name__ == "__main__":
    args = parse_command_line_arguments()
    coordinate_types = { 'd' : 'Direct',
                         'direct' : 'Direct',
                         'c' : 'Cartesian',
                         'cartesian' : 'Cartesian' }
    coordinate_type = coordinate_types[ args.coordinate_type ]
    # initialise
    poscar = Poscar()
    # read POSCAR file
    poscar.read_from( args.poscar )
    if args.supercell: # generate supercell
        if args.group:
            # check that if grouping is switched on, we are asking for a supercell that allows a "3D-chequerboard" pattern.
            for (i,axis) in zip( args.supercell, range(3) ):
                if i%2 == 1 and i > 1:
                    raise Exception( "odd supercell expansions != 1 are incompatible with automatic grouping" )
        poscar = poscar.replicate( *args.supercell, group=args.group )
    if args.bohr:
        poscar = poscar.in_bohr()
    # output to stdout
    if args.coordinates_only:
        poscar.output_coordinates_only( coordinate_type = coordinate_type, label = args.label )
    else:
        poscar.output( coordinate_type = coordinate_type, label = args.label )