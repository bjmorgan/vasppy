#! /usr/bin/env python3

from vasppy.summary import potcar_spec
import argparse

def parse_command_line_arguments():
    parser = argparse.ArgumentParser( description='Generate POTCAR specification based on hashing individual pseudopotential strings' )
    parser.add_argument( 'potcar', help="filename of the VASP POTCAR to be processed" )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_command_line_arguments()
    for p, ps in potcar_spec( args.potcar ).items():
        print( p, ps )
