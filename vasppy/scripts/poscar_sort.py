#! /usr/bin/env python3 

from vasppy.poscar import Poscar
import argparse
import copy
import numpy as np

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser( description='Reorders the atomic species in a VASP POSCAR file' )
    parser.add_argument( 'labels', nargs='+', help="atomic lables in the desired order" )
    parser.add_argument( 'poscar', help="filename of the VASP POSCAR to be processed" )
    args = parser.parse_args()
    return( args )

def main():
    args = parse_command_line_arguments()
    # initialise
    poscar = Poscar()
    # read POSCAR file
    poscar.read_from( args.poscar )
    # construct new Poscar instance with sorted atom types
    sorted_poscar = copy.deepcopy( poscar )
    sorted_poscar.atoms = []
    sorted_poscar.atom_numbers = []
    coordinate_list = []
    atoms = poscar.to_configuration().atoms
    for label in args.labels:
        if label in poscar.atoms:
            sorted_poscar.atoms.append( label )
            matched_atoms = [ atom for atom in atoms if atom.label == label ]
            sorted_poscar.atom_numbers.append( len( matched_atoms ) )
            coordinate_list.extend( [ atom.r for atom in matched_atoms ] )
        else:
            raise( ValueError( "'{}' atom label not found in {}".format( label, args.poscar ) ) )
    sorted_poscar.coordinates = np.array( coordinate_list )
    sorted_poscar.output( coordinate_type = poscar.coordinate_type )

if __name__ == "__main__":
    main()
