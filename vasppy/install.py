#! /usr/bin/env python3

import argparse
import os
import sys

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser( description="Installs exectuable python scripts\nDefault location is $HOME/bin" )
    parser.add_argument( 'name', nargs='*', help='name(s) of scripts to install', default=['all'] )
    parser.add_argument( '--prefix', help='specify directory to install files to. (default $HOME/bin/)' )
    parser.add_argument( '--list', action='store_true', help='list installable scripts' )
    parser.add_argument( '-f', '--force', action='store_true', help='force files to be overwritten on installation' )
    parser.add_argument( '-r', '--remove', action='store_true', help='remove files' )
    args = parser.parse_args()
    return( args )

installable_scripts = [ 'super', 'poscar_to_xtl', 'poscar_to_pimaim', 'poscar_sort', 'vasp_grid', 'xdatcar_to_rdf', 'xdatcar_to_disp' ]
home = os.path.expanduser("~")
install_dir = os.path.join( home, "bin" )
origin_dir = os.path.dirname( os.path.abspath( __file__ ) )

if __name__ == "__main__":
    args = parse_command_line_arguments()
    if( args.prefix ):
        install_dir = args.prefix
    if( args.list ):
        print( "\n".join( installable_scripts ) )
        sys.exit( 0 )
    scripts_to_install = installable_scripts if ( args.name[0] == 'all' ) else args.name
    if args.remove:
        for script in scripts_to_install:
            to_remove = os.path.join( install_dir, script )
            if os.path.isfile( to_remove ):
                os.remove( to_remove )
                print( script + ' unlinked' )
            else:
                print( to_remove + ' not found' )
        sys.exit( 0 )
    for script in scripts_to_install:
        to_install = os.path.join( install_dir, script )
        if os.path.isfile( to_install ):
            if args.force:
                os.remove( to_install )
                print( 'deleting '+to_install )
        print( "linking "+ script +'.py --> '+ to_install )
        try:
            os.symlink( os.path.join( origin_dir, script+'.py' ), to_install )
        except FileExistsError:
            print( '  '+script+' already exists in '+install_dir )
            print( '  linking failed' )
            sys.exit( -1 )


