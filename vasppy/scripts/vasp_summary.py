#! /usr/bin/env python3
"""Collect information about VASP calculations into YAML format for further processing.

Expects a series of directories (listed in ``result_dirs``) that each contain:
    vasprun.xml
    vaspmeta.yaml (additional metadata providing information about each calculation)
"""

import sys
if sys.platform != "win32":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)
from pathlib import Path
import yaml # type: ignore
import tqdm  # type: ignore
from multiprocessing import Pool
import argparse
from vasppy.summary import Summary, find_vasp_calculations
from vasppy.vaspmeta import VASPMeta


def get_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(description="Summarise a VASP calculation.")
    parser.add_argument(
        "-r",
        "--recursive",
        help="Recursively analyse directories.",
        action="store_true",
    )
    parser.add_argument(
        "-l", "--list", help="List supported data flags.", action="store_true"
    )
    parser.add_argument("-p", "--print", help="Specify data to parse.", nargs="*")
    parser.add_argument("-f", "--file", help="Specify a file to read data flags from.")
    parser.add_argument(
        "-c",
        "--check",
        help="Checks whether VASP directories contain vaspmeta.yaml and vasprun.xml files",
        action="store_true",
    )
    parser.add_argument(
        "-b",
        "--progress-bar",
        help="Show progress bar when parsing vasprun.xml files",
        action="store_true",
    )
    parser.add_argument(
        "-j",
        "--maxjobs",
        help="Maximum number of calculations to parse in parallel",
        type=int,
    )
    return parser.parse_args()


def get_summary(p: str) -> Summary:
    """Return a Summary object for the given path.

    Args:
        p: Path to a VASP calculation directory.

    Returns:
        Summary object for the calculation.
    """
    return Summary(p)


# This should really be set in the vasppy.Summary code, so that it can be tested to be consistent with the supported print methods.
# In fact, ideally the key, print method, and description would all be collected in a single object, which suggests writing this as a simple class.


def main() -> None:
    """Entry point for the vasp_summary command-line script."""
    supported_flags = Summary.supported_flags
    to_print = [
        "title",
        "status",
        "stoichiometry",
        "potcar",
        "plus_u",
        "energy",
        "lreal",
        "k-points",
        "functional",
        "encut",
        "ediffg",
        "ibrion",
        "converged",
        "version",
        "md5",
        "directory",
    ]
    titles = None
    args = get_args()
    if args.list:
        for k, v in supported_flags.items():
            print(f"{k.ljust(15)}: {v}")
        sys.exit()
    if args.file:
        with open(args.file, "r") as stream:
            settings = yaml.load(stream, Loader=yaml.SafeLoader)
        if "to_print" in settings:
            to_print = settings["to_print"]
        if "titles" in settings:
            titles = settings["titles"]
    if args.print:
        not_supported = [p for p in args.print if p not in supported_flags]
        if not_supported:
            raise ValueError(not_supported)
        else:
            to_print = args.print
    if args.recursive:
        path = sorted(find_vasp_calculations())
    else:
        path = ["."]
    if args.check:
        for p in path:
            vaspmeta = Path(f"{p}/vaspmeta.yaml")
            if not vaspmeta.is_file():
                print(f"{p} is missing vaspmeta.yaml")
            vasprun = Path(f"{p}/vasprun.xml")
            if not vasprun.is_file():
                print(f"{p} is missing vasprun.xml")
    else:
        if titles:
            # Only parse directories with matching vasp_meta titles
            matching_path = []
            for p in path:
                vm = VASPMeta.from_file(f"{p}/vaspmeta.yaml")
                if vm.title in titles:
                    matching_path.append(p)
            path = matching_path
        if args.maxjobs:
            with Pool(args.maxjobs) as pool:
                if args.progress_bar:
                    summaries = list(
                        tqdm.tqdm(pool.imap(get_summary, path), total=len(path))
                    )
                else:
                    summaries = pool.map(get_summary, path)
        else:
            path_iterator = tqdm.tqdm(path, unit="vasprun") if args.progress_bar else path
            summaries = [get_summary(p) for p in path_iterator]
        iterable = tqdm.tqdm(summaries, unit="records") if args.progress_bar else summaries
        for s in iterable:
            s.output(to_print=to_print)


if __name__ == "__main__":
    main()
