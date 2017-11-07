#! /usr/bin/env python3

import sys

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
    parser.add_argument( '-l', '--list', help="List supported data flags.", action='store_true' )
    parser.add_argument( '-p', '--print', help="Specify data to parse.", nargs='*' )
    parser.add_argument( '-f', '--file', help="Specify a file to read data flags from." )
    args = parser.parse_args()
    return args

# This should really be set in the vasppy.Summary code, so that it can be tested to be consistent with the supported print methods.
# In fact, ideally the key, print method, and description would all be collected in a single object, which suggests writing this as a simple class.
supported_flags = { 'status': 'Status',
                    'stoichiometry': 'Stoichiometry',
                    'potcar': 'POTCAR',
                    'plus_u': 'Dudarev +U parameters',
                    'energy': 'Energy',
                    'lreal': 'LREAL',
                    'k-points': 'k-points',
                    'functional': 'functional',
                    'encut': 'encut',
                    'ediffg': 'ediffg',
                    'ibrion': 'ibrion',
                    'converged': 'converged',
                    'md5': 'md5',
                    'directory': 'directory' }

to_print=[ 'title', 'status', 'stoichiometry', 'potcar', 'plus_u', 'energy', 'lreal', 'k-points', 'functional', 'encut', 'ediffg', 'ibrion', 'converged', 'version', 'md5', 'directory' ]

if __name__ == "__main__":
    args = get_args()
    if args.list:
        for k, v, in supported_flags.items():
            print( "{}: {}".format( k.ljust(15), v ) )
        sys.exit()
    if args.print:
        not_supported = [ p for p in args.print if p not in supported_flags ]
        if not_supported:
            raise ValueError( not_supported )
        else:
            to_print = args.print
    if args.recursive:
        path = sorted( find_vasp_calculations() )
    else:
        path = [ '.' ]
    for p in path:
        s = Summary( p )
        s.output( to_print=to_print )
