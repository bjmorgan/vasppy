#! /usr/bin/env python3

from vasppy.xdatcar import Xdatcar
from vasppy.rdf import Rdf
import argparse
import copy
import math

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument( 'xdatcar' )
    parser.add_argument( 'label', nargs = 2 )
    parser.add_argument( 'max_r', type = float )
    parser.add_argument( 'n_bins', type = int )
    args = parser.parse_args()
    return( args )

def main():
    args = parse_command_line_arguments()
    max_r = args.max_r
    number_of_bins = args.n_bins
    species_1 = args.label[ 0 ]
    species_2 = args.label[ 1 ]
    xdatcar = Xdatcar()
    xdatcar.read_from( args.xdatcar )
    rdf = Rdf( max_r = max_r, number_of_bins = number_of_bins )
    volume_scaling_factor = 4.0 * math.pi / ( 3.0 * xdatcar.poscar[0].cell.volume())
    for poscar in xdatcar.poscar:
        rdf += poscar.to_configuration().partial_rdf( species_1, species_2, max_r = max_r, number_of_bins = number_of_bins )
    print( rdf.normalised_data() )
    #[ print( dr, g_of_r / volume_scaling_factor ) for dr, g_of_r in rdf.normalised_data() ]

if __name__ == "__main__":
    main()
