import numpy as np
import re

def reciprocal_lattice_from_outcar( filename ): # from https://github.com/MaterialsDiscovery/PyChemia
    """
    Finds and returns the reciprocal lattice vectors, if more than
    one set present, it just returns the last one.
    Args:
        filename (Str): The name of the outcar file to be read

    Returns:
        List(Float): The reciprocal lattice vectors.
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

def final_energy_from_outcar( filename='OUTCAR' ):
    """
    Finds and returns the energy from a VASP OUTCAR file, by searching for the last `energy(sigma->0)` entry.

    Args:
        filename (Str, optional): OUTCAR filename. Defaults to 'OUTCAR'.

    Returns:
        (Float): The last energy read from the OUTCAR file.
    """
    with open( filename ) as f:
        outcar = f.read()
    energy_re = re.compile( "energy\(sigma->0\) =\s+([-\d\.]+)" )
    energy = float( energy_re.findall( outcar )[-1] )
    return energy

def vasp_version_from_outcar( filename='OUTCAR' ):
    """
    Returns the first line from a VASP OUTCAR file, to get the VASP source version string.

    Args:
        filename (Str, optional): OUTCAR filename. Defaults to 'OUTCAR'.

    Returns:
        (Str): The first line read from the OUTCAR file.
    """
    with open( filename ) as f:
        line = f.readline().strip()
    return line

def potcar_eatom_list_from_outcar( filename='OUTCAR' ):
    """
    Returns a list of EATOM values for the pseudopotentials used.

    Args:
        filename (Str, optional): OUTCAR filename. Defaults to 'OUTCAR'.

    Returns:
        (List(Float)): A list of EATOM values, in the order they appear in the OUTCAR.
    """
    with open( filename ) as f:
        outcar = f.read()
    eatom_re = re.compile( "energy of atom\s+\d+\s+EATOM=\s*([-\d\.]+)" )
    eatom = [ float( e ) for e in eatom_re.findall( outcar ) ]
    return eatom


def fermi_energy_from_outcar( filename='OUTCAR' ):
    """Finds and returns the fermi energy.
    Args:
    -filename: the name of the outcar file to be read

    Returns:
        (Float): The fermi energy as found in the OUTCAR 
    """
    outcar = open(filename, "r").read()
    # returns a match object
    fermi_energy = re.search(r"E-fermi\s*:\s*([-.\d]*)", outcar)
    # take the first group - group(0) contains entire match
    fermi_energy = float(fermi_energy.group(1))
    return fermi_energy
