#! /usr/bin/env python3

import sys
from pathlib import Path
import yaml

"""
Script for collecting information about VASP calculations into YAML format, for further processing.
Expects a series of directories (listed in `result_dirs`) that each contain:
    vasprun.xml
    vaspmeta.yaml (additional metadata providing information about each calculation)
"""

import argparse
from vasppy.summary import Summary, find_vasp_calculations
from vasppy.vaspmeta import VASPMeta

def get_args():
    parser = argparse.ArgumentParser( description='Summarise a VASP calculation.' )
    parser.add_argument( '-r', '--recursive', help='Recursively analyse directories.', action='store_true' )
    parser.add_argument( '-l', '--list', help="List supported data flags.", action='store_true' )
    parser.add_argument( '-p', '--print', help="Specify data to parse.", nargs='*' )
    parser.add_argument( '-f', '--file', help="Specify a file to read data flags from." )
    parser.add_argument( '-c', '--check', help="Checks whether VASP directories contain vaspmeta.yaml and vasprun.xml files", action='store_true' )
    args = parser.parse_args()
    return args

# This should really be set in the vasppy.Summary code, so that it can be tested to be consistent with the supported print methods.
# In fact, ideally the key, print method, and description would all be collected in a single object, which suggests writing this as a simple class.

def main():
    supported_flags = Summary.supported_flags
    to_print=[ 'title', 'status', 'stoichiometry', 'potcar', 'plus_u', 'energy', 'lreal', 'k-points', 'functional', 'encut', 'ediffg', 'ibrion', 'converged', 'version', 'md5', 'directory' ]
    titles = None
    args = get_args()
    if args.list:
        for k, v, in supported_flags.items():
            print( "{}: {}".format( k.ljust(15), v ) )
        sys.exit()
    if args.file:
        with open( args.file, 'r' ) as stream:
            settings = yaml.load( stream, Loader=yaml.SafeLoader ) 
        if 'to_print' in settings:
            to_print = settings['to_print']
        if 'titles' in settings:
            titles = settings['titles']
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
    if args.check:
        for p in path:
           vaspmeta = Path( '{}/vaspmeta.yaml'.format( p ) )
           if not vaspmeta.is_file(): 
               print( '{} is missing vaspmeta.yaml'.format( p ) ) 
           vasprun  = Path( '{}/vasprun.xml'.format( p ) )
           if not vasprun.is_file():
               print( '{} is missing vasprun.xml'.format( p ) ) 
    else:
        if titles:
            # Only parse directories with matching vasp_meta titles
            matching_path = []
            for p in path:
                vm = VASPMeta.from_file( '{}/vaspmeta.yaml'.format(p) )
                if vm.title in titles:
                    matching_path.append( p )
            path = matching_path
        for p in path:
            s = Summary( p )
            s.output( to_print=to_print )

if __name__ == "__main__":
    main()
