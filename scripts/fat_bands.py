#! /usr/bin/env python3

from vasppy import procar
import argparse

def orbitals_with_l( l ):
    to_return = { 's' : [ 0 ],
                  'p' : [ 1, 2, 3 ],
                  'd' : [ 4, 5, 6, 7, 8 ],
                  'f' : [ 9, 10, 11, 12, 13 ],
                  'all' : None }
    return to_return[ l ]
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument( '-i', '--ions', help='ion indices for band projection (default: sum over all ions)', nargs='+', type=int )
    parser.add_argument( '-s', '--spins', help='spin indices for band projection (default 0)', nargs='+', type=int, default=[ 0 ] )
    parser.add_argument( '-o', '--orbitals', help='orbital indices for band projection (default: sum over all orbitals)', nargs='+', type=int )
    parser.add_argument( '-e', '--efermi', help='set fermi energy as reference for energy scale', type=float, default=0.0 )
    parser.add_argument( '-l', '--l-angular-momentum', help='select all orbitals with angular momentum L for band projection. This supercedes the --orbitals option', choices=[ 's', 'p', 'd', 'f', 'all' ] )
    parser.add_argument( '-f', '--procar', help='PROCAR filename (default PROCAR)', type=str, default='PROCAR' )
    parser.add_argument( '--scaling', help='Energy scaling for band widths (default 0.2 eV)', type=float, default=0.2 )

    args = parser.parse_args()

    if args.l_angular_momentum:
        args.orbitals = orbitals_with_l( args.l_angular_momentum )

    pcar = procar.Procar()
    pcar.read_from_file( args.procar )
    pcar.print_weighted_band_structure( spins = args.spins, ions = args.ions, orbitals = args.orbitals, scaling = args.scaling, e_fermi = args.efermi )
