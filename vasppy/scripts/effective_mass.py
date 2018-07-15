#! /usr/bin/env python3

from vasppy import procar
from vasppy.outcar import reciprocal_lattice_from_outcar
import argparse

def minimum_length( nmin ):
    class MinimumLength( argparse.Action ):
        def __call__( self, parser, args, values, option_string=None ):
            if not nmin <= len( values ):
                msg = 'argument "{f}" requires at least {nmin} arguments'.format( f = self.dest, nmin = nmin )
                raise argparse.ArgumentError( self, msg )
            setattr( args, self.dest, values )
    return MinimumLength
            
def main():
    parser = argparse.ArgumentParser( description='Calculate an effective mass from a VASP PROCAR using a fitted quadratic' )
    parser.add_argument( '-k', '--k-points', help='index of k-points for calculating effective mass', nargs='+', type=int, required=True, action=minimum_length( 2 ) )
    parser.add_argument( '-b', '--band-index', help='index of band for calculating effective mass', type=int, required=True )
    parser.add_argument( '-f', '--procar', help='PROCAR filename (default PROCAR)', type=str, default='PROCAR' )
    parser.add_argument( '-v', '--verbose', help='Verbose output', action='store_true' )
    parser.add_argument( '-o', '--outcar', help='OUTCAR filename (default OUTCAR)', type=str, default='OUTCAR' )
    parser.add_argument( '-s', '--spin', help='select spin channel (default 1 / non-spin-polarised)', type=int, default='1' )
    args = parser.parse_args()

    reciprocal_lattice = reciprocal_lattice_from_outcar( 'OUTCAR' ) # Move reading the reciprocal lattice to procar.py

    pcar = procar.Procar()
    pcar.read_from_file( args.procar )
    effective_mass = pcar.effective_mass_calc( k_point_indices = args.k_points, 
                                               band_index = args.band_index, 
                                               reciprocal_lattice = reciprocal_lattice,
                                               spin = args.spin,
                                               printing = args.verbose )
    print( effective_mass )

if __name__ == '__main__':
    main()
