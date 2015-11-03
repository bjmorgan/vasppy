#! /usr/bin/env python3

import numpy as np
import argparse
from vasppy.poscar import Poscar
from vasppy.cell import Cell

def parse_command_line_arguments():
    parser = argparse.ArgumentParser( description = 'TODO' )
    parser.add_argument( '-l', '--labels', nargs='+', help="labels for each species", required = True )
    parser.add_argument( '-n', '--atom-numbers', nargs='+', help="atom numbers for each species", required = True, type = int )
    args = parser.parse_args()
    if len( args.labels ) != len( args.atom_numbers ):
        print( '\'labels\' and \'atom-numbers\' require matched argument numbers' )
        exit()
    return args 

def lines_to_numpy_array( data ):
    return np.array( [ [ float( s ) for s in line.split() ] for line in data ] )

def read_pimaim_restart( filename, number_of_atoms ):
    with open( filename, 'r' ) as f:
        file_data = f.readlines()

    cr_dump_log, vel_dump_log, chg_dump_log, full_dump_log = [ ( line.strip() == 'T' ) for line in file_data[:4] ]

    # this assumes coordinates, velocities, and dipoles are all present.
    # not sure what happens if atoms have qudrupoles, etc.
    coordinates = lines_to_numpy_array( file_data[ 4 : 4 + number_of_atoms ] ) 
    if vel_dump_log:
        velocities = lines_to_numpy_array( file_data[ 4 + number_of_atoms : 4 + number_of_atoms * 2 ] ) 
    else:
        velocities = None
    if chg_dump_log:
        dipoles = lines_to_numpy_array( file_data[ 4 + number_of_atoms * 2 : 4 + number_of_atoms * 3 ] ) 
    else:
        dipoles = None
    cell_matrix = lines_to_numpy_array( file_data[ -6: -3 ] )
    cell_lengths = lines_to_numpy_array( file_data[ -3: ] )
    full_cell_matrix = cell_matrix * cell_lengths
    # TODO! need to check this with a non-orthorhombic cell
    return( coordinates, velocities, dipoles, full_cell_matrix )

if __name__ == '__main__':
    filename = 'testout.rst'
    restart_file = True

    args = parse_command_line_arguments()
    number_of_atoms = sum( args.atom_numbers )
    coordinates, velocities, dipoles, full_cell_matrix = read_pimaim_restart( filename, number_of_atoms )

    poscar = Poscar()
    poscar.cell = Cell( full_cell_matrix ) # TODO: possibly this needs transposing?
    poscar.atoms = args.labels
    poscar.atom_numbers = args.atom_numbers
    poscar.coordinate_type = 'Cartesian'
    poscar.coordinates = coordinates

    poscar.output_as_xtl()
