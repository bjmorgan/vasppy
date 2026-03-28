from functools import cached_property
from lxml import etree  # type: ignore
from pymatgen.core import Structure  # type: ignore
from typing import TypedDict, cast
import numpy as np


class StructureData(TypedDict):
    """TypedDict for structure data returned by parse_structure."""

    lattice: list[list[float]]
    frac_coords: list[list[float]]
    selective_dynamics: list[list[bool]] | None


def parse_varray(
    varray: etree.Element,
) -> list[list[float]] | list[list[int]] | list[list[bool]]:
    """Parse ``<varray>`` data.

    Args:
        varray: xml ``<varray>`` element.

    Returns:
        A nested list of either float, int, or bool.
    """
    m: list[list[int]] | list[list[float]] | list[list[bool]]
    varray_type = varray.get("type", None)
    v_list = [v.text.split() for v in varray.findall("v") if v.text is not None]
    if varray_type == "int":
        m = [[int(number) for number in v] for v in v_list]
    elif varray_type == "logical":
        m = [[i == "T" for i in v] for v in v_list]
    else:
        m = [[float(number) for number in v] for v in v_list]
    return m


def parse_structure(structure: etree.Element) -> StructureData:
    """Parse ``<structure>`` data.

    Args:
        structure: xml ``<structure>`` element.

    Returns:
        Dictionary of structure data with the following keys:

        - ``lattice``: cell matrix (list[list[float]]).
        - ``frac_coords``: atom fractional coordinates (list[list[float]]).
        - ``selective_dynamics``: selective dynamics (list[list[bool]] or None).
    """
    crystal = structure.find("crystal")
    if crystal is None:
        raise ValueError("Truncated vasprun.xml: missing <crystal> element in <structure>")
    crystal_varray = crystal.find("varray")
    if crystal_varray is None:
        raise ValueError("Truncated vasprun.xml: missing <varray> in <crystal>")
    latt = parse_varray(crystal_varray)
    pos_element = structure.find("varray")
    if pos_element is None:
        raise ValueError("Truncated vasprun.xml: missing <varray> element in <structure>")
    pos = parse_varray(pos_element)
    sdyn = structure.find("varray/[@name='selective']")
    if sdyn is not None:
        sdyn = parse_varray(sdyn)
    structure_dict: StructureData = {
        "lattice": cast(list[list[float]], latt),
        "frac_coords": cast(list[list[float]], pos),
        "selective_dynamics": cast(list[list[bool]] | None, sdyn),
    }
    return structure_dict


def structure_from_structure_data(
    lattice: list[list[float]],
    atom_names: list[str],
    frac_coords: list[list[float]],
) -> Structure:
    """Generate a pymatgen Structure.

    Args:
        lattice: 3x3 cell matrix.
        atom_names: List of atom name strings.
        frac_coords: Nx3 list of fractional coordinates.

    Returns:
        A pymatgen Structure object.
    """
    return Structure(
        lattice=lattice,
        species=atom_names,
        coords=frac_coords,
        coords_are_cartesian=False,
    )


class Vasprun:
    """Object for parsing vasprun.xml data.

    Attributes:
        atom_names: List of atom name strings.
        structures: List of structures as pymatgen Structure objects.
        frac_coords: timesteps x atoms x 3 numpy array of fractional coordinates.
        cart_coords: timesteps x atoms x 3 numpy array of Cartesian coordinates.
        forces: timesteps x atoms x 3 numpy array of forces (if present).

    Examples:
        >>> vasprun = Vasprun('vasprun.xml')
        >>> cart_coords = vasprun.cart_coords
        >>> forces = vasprun.forces
    """

    def __init__(self, filename: str) -> None:
        """Initialise a Vasprun object from a vasprun.xml file.

        Args:
            filename: The vasprun.xml filename.
        """
        doc = etree.parse(filename)
        self.doc = doc.getroot()

    @cached_property
    def structures(self) -> list[Structure]:
        """Parse and return all structures from the vasprun.xml file.

        Returns:
            A list of pymatgen Structure objects.
        """
        return self.parse_structures()

    @cached_property
    def atom_names(self) -> list[str]:
        """Parse and return atom names from the vasprun.xml file.

        Returns:
            A list of atom name strings.
        """
        return self.parse_atom_names()

    def parse_atom_names(self) -> list[str]:
        """Return a list of atom names for the atoms in this calculation.

        Returns:
            A list of atom name strings.

        Raises:
            ValueError: If no atominfo element is found in the file.
            ValueError: If no atom names are found in the atominfo element.
        """
        atominfo = self.doc.find("atominfo")
        if atominfo is None:
            raise ValueError("No atominfo found in file")
        atom_names: list[str] = []
        for array in atominfo.findall("array"):
            if array.attrib["name"] == "atoms":
                atom_names = [rc.find("c").text.strip() for rc in array.find("set")]
        if not atom_names:
            raise ValueError("No atomname found in file")
        return atom_names

    def parse_structures(self) -> list[Structure]:
        """Return a list of pymatgen Structures for this calculation.

        Returns:
            A list of pymatgen Structure objects.
        """
        structures: list[Structure] = []
        for i, child in enumerate(self.doc.iterfind("calculation")):
            elem = child.find("structure")
            if elem is None:
                raise ValueError(
                    f"Truncated vasprun.xml: missing <structure> in <calculation> {i}"
                )
            structure_data = parse_structure(elem)
            structures.append(
                structure_from_structure_data(
                    lattice=structure_data["lattice"],
                    atom_names=self.atom_names,
                    frac_coords=structure_data["frac_coords"],
                )
            )
        return structures

    @property
    def frac_coords(self) -> np.ndarray:
        """Fractional coordinates from each calculation structure.

        Returns:
            timesteps x atoms x 3 numpy array of fractional coordinates.
        """
        return np.array([s.frac_coords for s in self.structures])

    @property
    def cart_coords(self) -> np.ndarray:
        """Cartesian coordinates from each calculation structure.

        Returns:
            timesteps x atoms x 3 numpy array of Cartesian coordinates.
        """
        return np.array([s.cart_coords for s in self.structures])

    @property
    def forces(self) -> np.ndarray | None:
        """Cartesian forces from each calculation structure (if present in the vasprun XML).

        Returns:
            timesteps x atoms x 3 numpy array of Cartesian forces if forces are
            included in the vasprun XML, otherwise ``None``.
        """
        forces = []
        for child in self.doc.iterfind("calculation"):
            elem = child.find("varray/[@name='forces']")
            if elem is not None:
                forces.append(parse_varray(elem))
        return np.array(forces) if forces else None
