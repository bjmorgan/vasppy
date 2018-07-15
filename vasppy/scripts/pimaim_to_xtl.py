#! /usr/bin/env python3

import numpy as np
import argparse
from vasppy.poscar import Poscar
from vasppy.cell import Cell
import vasppy.pimaim

def parse_command_line_arguments():
    parser = argparse.ArgumentParser( description = 'TODO' )
    parser.add_argument( '-f', '--filename', help="PIMAIM restart file filename" )
    parser.add_argument( '-l', '--labels', nargs='+', help="labels for each species", required = True )
    parser.add_argument( '-n', '--atom-numbers', nargs='+', help="atom numbers for each species", required = True, type = int )
    args = parser.parse_args()
    if len( args.labels ) != len( args.atom_numbers ):
        print( '\'labels\' and \'atom-numbers\' require matched argument numbers' )
        exit()
    return args 

def main():
    args = parse_command_line_arguments()
    if args.filename is None:
        args.filename = 'testout.rst'
    poscar = vasppy.pimaim.poscar_from_pimaim_restart( args.filename, args.atom_numbers, args.labels )
    poscar.output_as_xtl()

if __name__ == '__main__':
    main()
