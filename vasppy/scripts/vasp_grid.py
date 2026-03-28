#! /usr/bin/env python3
"""z-projection of a VASP grid format file."""

from vasppy import grid
import argparse


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace with ``gridfile``, ``projection``,
        and ``orthorhombic`` attributes.
    """
    parser = argparse.ArgumentParser(
        description="z-projection of a VASP (grid format) file"
    )
    parser.add_argument(
        "gridfile", help="filename of the VASP (grid format) file to be processed"
    )
    parser.add_argument(
        "-p",
        "--projection",
        choices=["x", "y", "z"],
        help="output averaged projection perpendicular to [x,y,z]",
    )
    parser.add_argument(
        "-o",
        "--orthorhombic",
        help="map grid points onto an orthorhombic (non-space filling) grid",
        action="store_true",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for the vasp_grid command-line script."""
    args = parse_command_line_arguments()
    vgrid = grid.Grid.from_file(args.gridfile)
    if args.orthorhombic:
        vgrid = vgrid.interpolate_to_orthorhombic_grid(vgrid.dimensions)
    if args.projection:
        index = grid.Grid.projections[args.projection]
        grid_spacing = vgrid.structure.lattice.lengths[index] / vgrid.dimensions[index]
        for i, av in enumerate(vgrid.average(normal_axis_label=args.projection)):
            print(i * grid_spacing, av)


if __name__ == "__main__":
    main()
