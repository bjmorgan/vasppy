# routines for interacting with PIMAIM data / filetypes

import numpy as np
from vasppy.poscar import Poscar
from vasppy.cell import Cell

def lines_to_numpy_array( data ):
    return np.array( [ [ float( s ) for s in line.split() ] for line in data ] )

def read_restart_file( filename, number_of_atoms ):
    with open( filename, 'r' ) as f:
        file_data = f.readlines()

    cr_dump_log, vel_dump_log, chg_dump_log, full_dump_log = [ ( line.strip() == 'T' ) for line in file_data[:4] ]

    # this assumes coordinates, velocities, and dipoles are all present.
    # not sure what happens if atoms have qudrupoles, etc.
    coordinates = lines_to_numpy_array( file_data[ 4 : 4 + number_of_atoms ] ) * 0.52918 # convert bohr to Angstroms
    if vel_dump_log:
        velocities = lines_to_numpy_array( file_data[ 4 + number_of_atoms : 4 + number_of_atoms * 2 ] )
    else:
        velocities = None
    if chg_dump_log:
        dipoles = lines_to_numpy_array( file_data[ 4 + number_of_atoms * 2 : 4 + number_of_atoms * 3 ] )
    else:
        dipoles = None
    cell_matrix = lines_to_numpy_array( file_data[ -6: -3 ] )
    cell_lengths = lines_to_numpy_array( file_data[ -3: ] ) * 0.52918 # convert bohr to Angstroms
    full_cell_matrix = cell_matrix * cell_lengths
    # TODO! need to check this with a non-orthorhombic cell
    return( coordinates, velocities, dipoles, full_cell_matrix )

def poscar_from_pimaim_restart( filename, atom_numbers, atom_labels ):
    number_of_atoms = sum( atom_numbers )
    coordinates, velocities, dipoles, full_cell_matrix = read_restart_file( filename, number_of_atoms )

    poscar = Poscar()
    poscar.cell = Cell( full_cell_matrix ) # TODO: possibly this needs transposing
    poscar.atoms = atom_labels
    poscar.atom_numbers = atom_numbers
    poscar.coordinate_type = 'Direct'
    poscar.coordinates = poscar.cell.cartesian_to_fractional_coordinates( coordinates )
 
    return poscar
