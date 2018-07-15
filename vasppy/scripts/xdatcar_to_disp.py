#! /usr/bin/env python3

from vasppy.xdatcar import Xdatcar
import argparse
import copy

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument( 'xdatcar' )
    args = parser.parse_args()
    return( args )

def main():
    args = parse_command_line_arguments()
    xdatcar = Xdatcar()
    xdatcar.read_from( args.xdatcar )
    poscar1 = copy.deepcopy( xdatcar.poscar[0] )
    for p in xdatcar.poscar[1:]:
        poscar2 = p
        for i, c in enumerate( poscar1.coordinates ):
            print( ' '.join( [ str(e) for e in ( poscar1.cell.minimum_image( poscar1.coordinates[i], poscar2.coordinates[i] ).dot( poscar1.cell.matrix ) ) ] ) )
        poscar1 = copy.deepcopy( p )

if __name__ == "__main__":
    main()
