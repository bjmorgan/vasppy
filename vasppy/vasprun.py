from lxml import etree # type: ignore
from typing import List

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

