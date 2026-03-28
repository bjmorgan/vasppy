"""Unit conversion constants derived from CODATA physical constants."""

from scipy.constants import physical_constants, angstrom  # type: ignore

angstrom_to_bohr: float = physical_constants["atomic unit of length"][0] / angstrom
ev_to_hartree: float = physical_constants["electron volt-hartree relationship"][0]
