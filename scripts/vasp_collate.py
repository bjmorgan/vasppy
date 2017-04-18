#! /usr/bin/env python3 

from vasppy.poscar import Poscar
from vasppy.outcar import final_energy_from_outcar
from vasppy.vaspmeta import VASPMeta

import argparse

if __name__ == "__main__":
#    args = parse_command_line_arguments()
    # initialise
    poscar = Poscar()
    # read POSCAR file
    poscar.read_from( 'POSCAR' )
    stoich = poscar.stoichiometry
    energy = final_energy_from_outcar( 'OUTCAR' )
    meta = VASPMeta.from_file( 'vaspmeta.yaml' )

    print( "title: {}".format( meta.title ) )
    print( "stoichiometry:" )
    for element in stoich:
        print( "    - {}: {}".format( element, stoich[ element ] ) )
    print( "status: {}".format( meta.status ) )
    print( "energy: {}".format( energy ) )
