import numpy as np
import re

def reciprocal_lattice_from_outcar( filename ): # from https://github.com/MaterialsDiscovery/PyChemia
    """Finds and return the reciprocal lattice vectors, if more than
    one set present, it just returns the last one.
    Args:
    -filename: the name of the outcar file  to be read
    """
    outcar = open(filename, "r").read()
    # just keeping the last component
    recLat = re.findall(r"reciprocal\s*lattice\s*vectors\s*([-.\s\d]*)",
                        outcar)[-1]
    recLat = recLat.split()
    recLat = np.array(recLat, dtype=float)
    # up to now I have, both direct and rec. lattices (3+3=6 columns)
    recLat.shape = (3, 6)
    recLat = recLat[:, 3:]
    return recLat
