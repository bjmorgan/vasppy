#! /usr/bin/env python3

from vasppy import grid
import argparse

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser( description='z-projection of a VASP (grid format) file' )
    parser.add_argument( 'gridfile', help="filename of the VASP (grid format) file to be processed" )
    parser.add_argument( '-p', '--projection', choices=[ 'x', 'y', 'z' ], help="output averaged projection perpendicular to [x,y,z]" )
    parser.add_argument( '-o', '--orthorhombic', help='map grid points onto an orthorhombic (non-space filling) grid', action = 'store_true' )
    args = parser.parse_args()
    return args

def main():
    args = parse_command_line_arguments()
    vgrid = grid.Grid()
    vgrid.read_from_filename( args.gridfile )
    if args.orthorhombic:
        vgrid = vgrid.interpolate_to_orthorhombic_grid( vgrid.dimensions )
    if args.projection:
        index = grid.Grid.projections[ args.projection ]
        grid_spacing = vgrid.poscar.cell_lengths()[ index ] / vgrid.dimensions[ index ]
        [ print( i * grid_spacing, av ) for i, av in enumerate( vgrid.average( normal_axis_label = args.projection ) ) ]

if __name__ == "__main__":
    main()
