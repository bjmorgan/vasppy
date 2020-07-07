import numpy as np  # type: ignore
import re
from pymatgen.io.vasp.outputs import Outcar  # type: ignore

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
    """Finds and returns the Fermi energy.

    Args:
        filename (:obj:'str', optional): the name of the ``OUTCAR`` file to be read. Default is `OUTCAR`.

    Returns:
        (Float): The Fermi energy as found in the ``OUTCAR`` file.

    """
    outcar = open(filename, "r").read()
    # returns a match object
    fermi_energy = re.search(r"E-fermi\s*:\s*([-.\d]*)", outcar)
    # take the first group - group(0) contains entire match
    fermi_energy = float(fermi_energy.group(1))
    return fermi_energy

def forces_from_outcar( filename='OUTCAR' ):
    """Finds and returns forces from the OUTCAR file.
      
    Args:
        filename (:obj:'str', optional): the name of the ``OUTCAR`` file to be read. Default is `OUTCAR`.

    Returns:
        (np.array): The force as found in the ``OUTCAR`` file, as a NSTEPS x NIONS x 3 numpy array.

    """
    outcar = Outcar(filename)
    forces = outcar.read_table_pattern(
        header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
        row_pattern=r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        footer_pattern=r"\s--+",
        postprocess=lambda x: float(x),
        last_one_only=False
    )
    return np.array( forces )

def coords_from_outcar( filename='OUTCAR' ):
    """Finds and returns Cartesian coordinates from the OUTCAR file.
      
    Args:
        filename (:obj:'str', optional): the name of the ``OUTCAR`` file to be read. Default is `OUTCAR`.

    Returns:
        (np.array): The Cartesian coordinates as found in the ``OUTCAR`` file, as a NSTEPS x NIONS x 3 numpy array.

    """
    outcar = Outcar(filename)
    coords = outcar.read_table_pattern(
        header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
        row_pattern=r"\s+[+-]?(\d+\.\d+)\s+[+-]?(\d+\.\d+)\s+[+-]?(\d+\.\d+)\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+",
        footer_pattern=r"\s--+",
        postprocess=lambda x: float(x),
        last_one_only=False
    )
    return np.array( coords )
      
