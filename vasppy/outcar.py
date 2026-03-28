import numpy as np
import re
from pymatgen.io.vasp.outputs import Outcar

_RECIP_LAT_RE = re.compile(r"reciprocal\s*lattice\s*vectors\s*([-.\s\d]*)")
_ENERGY_RE = re.compile(r"energy\(sigma->0\) =\s+([-\d\.]+)")
_EATOM_RE = re.compile(r"energy of atom\s+\d+\s+EATOM=\s*([-\d\.]+)")
_FERMI_RE = re.compile(r"E-fermi\s*:\s*([-.\d]*)")


def reciprocal_lattice_from_outcar(
    filename: str,
) -> np.ndarray:  # from https://github.com/MaterialsDiscovery/PyChemia
    """Find and return the reciprocal lattice vectors from an OUTCAR file.

    If more than one set is present, the last one is returned.

    Args:
        filename: The name of the OUTCAR file to be read.

    Returns:
        A 3x3 numpy array of the reciprocal lattice vectors.
    """
    with open(filename) as f:
        outcar = f.read()
    # just keeping the last component
    matches = _RECIP_LAT_RE.findall(outcar)
    if not matches:
        raise ValueError(f"Reciprocal lattice vectors not found in {filename}")
    rec_lat = matches[-1]
    rec_lat = rec_lat.split()
    rec_lat = np.array(rec_lat, dtype=float)
    # up to now we have both direct and reciprocal lattices (3+3=6 columns)
    rec_lat.shape = (3, 6)
    rec_lat = rec_lat[:, 3:]
    return rec_lat


def final_energy_from_outcar(filename: str = "OUTCAR") -> float:
    """Find and return the energy from a VASP OUTCAR file.

    Searches for the last ``energy(sigma->0)`` entry.

    Args:
        filename: OUTCAR filename. Defaults to ``'OUTCAR'``.

    Returns:
        The last energy read from the OUTCAR file.
    """
    with open(filename) as f:
        outcar = f.read()
    matches = _ENERGY_RE.findall(outcar)
    if not matches:
        raise ValueError(f"Energy not found in {filename}")
    energy = float(matches[-1])
    return energy


def vasp_version_from_outcar(filename: str = "OUTCAR") -> str:
    """Return the VASP source version string from an OUTCAR file.

    Reads the first line of the file.

    Args:
        filename: OUTCAR filename. Defaults to ``'OUTCAR'``.

    Returns:
        The first line of the OUTCAR file (stripped of surrounding whitespace).
    """
    with open(filename) as f:
        line = f.readline().strip()
    return line


def potcar_eatom_list_from_outcar(filename: str = "OUTCAR") -> list[float]:
    """Return a list of EATOM values for the pseudopotentials used.

    Args:
        filename: OUTCAR filename. Defaults to ``'OUTCAR'``.

    Returns:
        A list of EATOM values, in the order they appear in the OUTCAR.

    Raises:
        ValueError: If no EATOM values are found in the file.
    """
    with open(filename) as f:
        outcar = f.read()
    eatom = [float(e) for e in _EATOM_RE.findall(outcar)]
    if not eatom:
        raise ValueError(f"No EATOM values found in {filename}")
    return eatom


def fermi_energy_from_outcar(filename: str = "OUTCAR") -> float:
    """Find and return the Fermi energy from an OUTCAR file.

    Args:
        filename: The name of the ``OUTCAR`` file to be read. Defaults to ``'OUTCAR'``.

    Returns:
        The Fermi energy as found in the ``OUTCAR`` file.
    """
    with open(filename) as f:
        outcar = f.read()
    # returns a match object
    fermi_energy_match = _FERMI_RE.search(outcar)
    if fermi_energy_match is None:
        raise ValueError("Fermi energy not found in OUTCAR file.")
    # take the first group — group(0) contains the entire match
    return float(fermi_energy_match.group(1))


def forces_from_outcar(
    filename: str = "OUTCAR", last_one_only: bool = False
) -> np.ndarray:
    """Find and return forces from the OUTCAR file.

    Args:
        filename: The name of the ``OUTCAR`` file to be read. Defaults to ``'OUTCAR'``.
        last_one_only: If ``True``, return only the last ionic step. Defaults to
            ``False``.

    Returns:
        The forces as found in the ``OUTCAR`` file.
        If ``last_one_only`` is ``False``: an NSTEPS x NIONS x 3 numpy array.
        If ``last_one_only`` is ``True``: an NIONS x 3 numpy array.
    """
    outcar = Outcar(filename)
    forces = outcar.read_table_pattern(
        header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
        row_pattern=r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
        footer_pattern=r"\s--+",
        postprocess=float,
        last_one_only=last_one_only,
    )
    return np.array(forces)


def coords_from_outcar(filename: str = "OUTCAR") -> np.ndarray:
    """Find and return Cartesian coordinates from the OUTCAR file.

    Args:
        filename: The name of the ``OUTCAR`` file to be read. Defaults to ``'OUTCAR'``.

    Returns:
        The Cartesian coordinates as found in the ``OUTCAR`` file, as an
        NSTEPS x NIONS x 3 numpy array.
    """
    outcar = Outcar(filename)
    coords = outcar.read_table_pattern(
        header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
        row_pattern=r"\s+[+-]?(\d+\.\d+)\s+[+-]?(\d+\.\d+)\s+[+-]?(\d+\.\d+)\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+",
        footer_pattern=r"\s--+",
        postprocess=float,
        last_one_only=False,
    )
    return np.array(coords)
