#! /usr/bin/env python3

"""
Script for collecting information about VASP calculations into YAML format, for further processing.
Expects a series of directories (listed in `result_dirs`) that each contain:
    vasprun.xml
    vaspmeta.yaml (additional metadata providing information about each calculation)
"""

from vasppy.summary import Summary, find_vasp_calculations

to_print=[ 'title', 'status', 'stoichiometry', 'potcar', 'plus_u', 'energy', 'k-points', 'functional', 'encut', 'ediffg', 'ibrion', 'converged', 'md5', 'directory' ]

if __name__ == "__main__":
    for path in sorted( find_vasp_calculations() ):
        s = Summary( path )
        s.output( to_print=to_print )
