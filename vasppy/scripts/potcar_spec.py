#! /usr/bin/env python3

from vasppy.summary import potcar_spec
import argparse

def parse_command_line_arguments():
    parser = argparse.ArgumentParser( description='Generate POTCAR specification based on hashing individual pseudopotential strings' )
    parser.add_argument('potcar', help="filename of the VASP POTCAR to be processed", nargs='?', default='POTCAR' )
    parser.add_argument('--hash', help="return the md5 hashes of the individual pseudopotential strings", action='store_true') 
    args = parser.parse_args()
    return args

def main():
    args = parse_command_line_arguments()
    if args.hash:
        hashes = {}
        for p, md5hash in potcar_spec(args.potcar, return_hashes=True).items():
            hashes[p] = md5hash
    for p, ps in potcar_spec(args.potcar).items():
        if args.hash:
            print(p, ps, hashes[p])
        else:
            print(p, ps, hashes[p])

if __name__ == '__main__':
    main()
