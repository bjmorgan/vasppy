#! /usr/bin/env python3

from pymatgen.io.vasp import Xdatcar  # type: ignore
from vasppy.rdf import RadialDistributionFunction
import argparse
import copy
import math
import numpy as np  # type: ignore

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument( 'xdatcar' )
    parser.add_argument( 'label', nargs=2 )
    parser.add_argument( 'max_r', type=float )
    parser.add_argument( 'n_bins', type=int )
    args = parser.parse_args()
    return( args )

def main():
    args = parse_command_line_arguments()
    max_r = args.max_r
    number_of_bins = args.n_bins
    species_1 = args.label[ 0 ]
    species_2 = args.label[ 1 ]
    xdatcar = Xdatcar( args.xdatcar )
    indices_i = [ i for i, s in enumerate(xdatcar.structures[0]) 
                  if s.species_string == species_1 ]
    if not indices_i:
        raise ValueError( f'No species {species_1} found' )
    indices_j = [ i for i, s in enumerate(xdatcar.structures[0]) 
                  if s.species_string == species_2 ]
    if not indices_j:
        raise ValueError( f'No species {species_2} found' )
    rdf = RadialDistributionFunction( xdatcar.structures, indices_i, indices_j, 
                                      number_of_bins, 0.0, max_r )
    for x, y in zip( rdf.r, rdf.rdf ):
        print(x, y)

if __name__ == "__main__":
    main()
