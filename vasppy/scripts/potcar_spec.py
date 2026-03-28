#! /usr/bin/env python3
"""Generate a POTCAR specification based on hashing individual pseudopotential strings."""

import sys
if sys.platform != "win32":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

from vasppy.summary import potcar_spec
import argparse


def parse_command_line_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace with ``potcar`` filename and ``hash`` flag.
    """
    parser = argparse.ArgumentParser(
        description="Generate POTCAR specification based on hashing individual pseudopotential strings"
    )
    parser.add_argument(
        "potcar",
        help="filename of the VASP POTCAR to be processed",
        nargs="?",
        default="POTCAR",
    )
    parser.add_argument(
        "--hash",
        help="return the md5 hashes of the individual pseudopotential strings",
        action="store_true",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for the potcar_spec command-line script."""
    args = parse_command_line_arguments()
    names, datasets = potcar_spec(args.potcar)
    if args.hash:
        _, hashes = potcar_spec(args.potcar, return_hashes=True)
        for name, dataset, md5hash in zip(names, datasets, hashes):
            print(name, dataset, md5hash)
    else:
        for name, dataset in zip(names, datasets):
            print(name, dataset)


if __name__ == "__main__":
    main()
