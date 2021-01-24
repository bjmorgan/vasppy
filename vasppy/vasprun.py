from lxml import etree # type: ignore
from typing import List, Union, TypeVar
TNum = TypeVar('TNum', int, float)

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

    def atom_names(self) -> List[str]:
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

