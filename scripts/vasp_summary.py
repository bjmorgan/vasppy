#! /usr/bin/env python3

"""
Script for collecting information about VASP calculations into YAML format, for further processing.
Expects a series of directories (listed in `result_dirs`) that each contain:
    vasprun.xml
    vaspmeta.yaml (additional metadata providing information about each calculation)
"""

import argparse
from vasppy.summary import Summary, find_vasp_calculations

def get_args():
    parser = argparse.ArgumentParser( description='Summarise a VASP calculation.' )
    parser.add_argument( '-r', '--recursive', help='Recursively analyse directories.', action='store_true' )
    args = parser.parse_args()
    return args

to_print=[ 'title', 'status', 'stoichiometry', 'potcar', 'plus_u', 'energy', 'k-points', 'functional', 'encut', 'ediffg', 'ibrion', 'converged', 'md5', 'directory' ]

if __name__ == "__main__":
    args = get_args()
    if args.recursive:
        path = sorted( find_vasp_calculations() )
    else:
        path = [ '.' ]
    for p in path:
        s = Summary( p )
        s.output( to_print=to_print )
