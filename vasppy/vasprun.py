from lxml import etree # type: ignore
from typing import List, Union, Optional, Any, Dict
from pymatgen import Structure # type: ignore

def parse_varray(varray: etree.Element) -> Union[List[List[float]], 
                                                 List[List[int]],
                                                 List[List[bool]]]:
    """Parse <varray> data.

    Args:
        varray (etree.Element): xml <varray> element.

    Returns:
        (list(list): A nested list of either float, int, or bool.

    """
    m: Union[List[List[int]], List[List[float]], List[List[bool]]]
    varray_type = varray.get("type", None)
    v_list = [v.text.split() for v in varray.findall("v")] 
    if varray_type == 'int':
        m = [[int(number) for number in v] for v in v_list]
    elif varray_type == 'logical':
        m = [[i == "T" for i in v] for v in v_list]
    else:
        m = [[float(number) for number in v] for v in v_list]
    return m

def parse_structure(structure: etree.Element) -> Dict[str, Any]:
    """Parse <structure> data..

    Args:
        structure (etree.Element): xml <structure> element.

    Returns:
        (dict): Dictionary of structure data:
            `lattice`: cell matrix (list(list(float))..
            `frac_coords`: atom fractional coordinates (list(list(float)).
            `selective_dynamics`: selective dynamics (bool|None).

    """
    latt = parse_varray(structure.find("crystal").find("varray"))
    pos = parse_varray(structure.find("varray"))
    sdyn = structure.find("varray/[@name='selective']")
    structure_dict = {'lattice': latt,
                      'frac_coords': pos,
                      'selective_dynamics': sdyn}
    return structure_dict

def structure_from_structure_data(lattice: List[List[float]],
                                  atom_names: List[str],
                                  frac_coords: List[List[float]]) -> Structure:
    """Generate a pymatgen Structure.

    Args:
        lattice (list(list(float)): 3x3 cell matrix.
        atom_names (list(str)): list of atom name strings.
        frac_coords (list(list(float): Nx3 list of fractional coordinates.

    Returns:
        (pymatgen.Structure)

    """
    structure = Structure(lattice=lattice,
                          species=atom_names,
                          coords=frac_coords,
                          coords_are_cartesian=False)
    return structure

class Vasprun:
    """Object for parsing vasprun.xml data."""

    def __init__(self,
                 filename: str) -> None: 
        """Initialise a Vasprun object from a vasprun.xml file.

        Args:
            filename (str): The vasprun.xml filename.

        Returns:
            None

        """
        doc = etree.parse(filename)
        self.doc = doc.getroot()
        self._atom_names = None # type: Optional[List[str]]
        self._structures = None # type: Optional[List[Structure]]

    @property
    def structures(self) -> List[Structure]:
        """Getter for structures attribut.

        Returns:
            (list(pymatgen.Structure)): A list of pymatgen Structure objects.

        Notes:
            When first called this parses the vasprun XML data and
            caches the result.

        """
        if not self._structures:
            self._structures = self.parse_structures()
        return self._structures
 
    @property
    def atom_names(self) -> List[str]:
        """Getter for atom_names attribute.

        Returns:
            (list(str)): A list of atom name strings.

        Notes:
            When first called this parses the vasprun XML data and
            caches the result.

        """
        if not self._atom_names:
            self._atom_names = self.parse_atom_names()
        return self._atom_names

    def parse_atom_names(self) -> List[str]:
        """Return a list of atom names for the atoms in this calculation.

        Args:
            None

        Returns:
            (list(str))

        """
        atominfo = self.doc.find("atominfo")
        if atominfo is None:
            raise ValueError("No atominfo found in file")
        atom_names = []
        for array in atominfo.findall("array"):
            if array.attrib["name"] == "atoms":
                atom_names = [rc.find("c").text.strip() for rc in array.find("set")]
        if not atom_names:
            raise ValueError("No atomname found in file")
        return atom_names 

    def parse_structures(self) -> List[Structure]:
        """Returns a list of pymatgen Structures for this calculation.

        Args:
            None

        Returns:
            (list(pymatgen.Structure))

        """
        structures = []
        for elem in self.doc.iterfind("structure"):
            structure_data = parse_structure(elem)
            structures.append(
                structure_from_structure_data(
                    lattice=structure_data['lattice'],
                    atom_names=self.atom_names,
                    frac_coords=structure_data['frac_coords']
                )
            )
        return structures
