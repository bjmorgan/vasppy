#! /usr/bin/env python3

from vasppy import procar
from vasppy.outcar import reciprocal_lattice_from_outcar
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument( '-k', '--k-points', help='index of k-points for calculating effective mass', nargs=2, type=int, required=True )
    parser.add_argument( '-b', '--band-index', help='index of band for calculating effective mass', type=int, required=True )
    parser.add_argument( '-f', '--procar', help='PROCAR filename (default PROCAR)', type=str, default='PROCAR' )
    args = parser.parse_args()

    reciprocal_lattice = reciprocal_lattice_from_outcar( 'OUTCAR' ) # Move reading the reciprocal lattice to procar.py

    pcar = procar.Procar()
    pcar.read_from_file( args.procar )
    pcar.effective_mass_calc( k_point_indices = args.k_points, 
                              band_index = args.band_index, 
                              reciprocal_lattice = reciprocal_lattice,
                              printing = True )

