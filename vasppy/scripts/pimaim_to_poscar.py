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

def read_pimaim_restart( filename ):
    with open( filename, 'r' ) as f:
        file_data = f.readlines()

    cr_dump_log, vel_dump_log, chg_dump_log, full_dump_log = [ bool( line.strip() ) for line in file_data[:4] ]

    # figure out how many atoms are in this calculation
    data_per_line = [ len( line.split() ) for line in file_data ] 

    number_of_atomic_data_lines = next( index for index, d in enumerate( data_per_line[4:] ) if d == 1 )
    number_of_atomic_data_types = sum( [ cr_dump_log, vel_dump_log, chg_dump_log ] )
    number_of_atoms = int( number_of_atomic_data_lines / number_of_atomic_data_types )
    # this assumes coordinates, velocities, and dipoles are all present.
    # not sure what happens if atoms have qudrupoles, etc.
    coordinates = lines_to_numpy_array( file_data[ 4 : 4 + number_of_atoms ] ) 
    velocities  = lines_to_numpy_array( file_data[ 4 + number_of_atoms : 4 + number_of_atoms * 2 ] ) 
    dipoles     = lines_to_numpy_array( file_data[ 4 + number_of_atoms * 2 : 4 + number_of_atoms * 3 ] ) 
    cell_matrix = lines_to_numpy_array( file_data[ -6: -3 ] )
    cell_lengths = lines_to_numpy_array( file_data[ -3: ] )
    full_cell_matrix = cell_matrix * cell_lengths
    # TODO! need to check this with a non-orthorhombic cell
    return( coordinates, velocities, dipoles, full_cell_matrix )

def main():
    filename = 'testout.rst'
    restart_file = True

    args = parse_command_line_arguments()
    coordinates, velocities, dipoles, full_cell_matrix = read_pimaim_restart( filename )
    assert( sum( args.atom_numbers ) == len( coordinates ) )
    poscar = Poscar()
    poscar.cell = Cell( full_cell_matrix ) # TODO: possibly this needs transposing?
    poscar.atoms = args.labels
    poscar.atom_numbers = args.atom_numbers
    poscar.coordinate_type = 'Cartesian'
    poscar.coordinates = coordinates
    poscar.output()

if __name__ == '__main__':
    main()
