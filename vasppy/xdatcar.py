from vasppy.poscar import Poscar
import re
import copy
import numpy as np


class Xdatcar:
    """Class for parsing and working with VASP XDATCAR files."""

    lines_offset = 9

    def __init__(self):
        """
        Initialise a Xdatcar object.

        Args:
            None

        Returns:
            None
        """
        self.poscar = []
        self.poscar.append(Poscar())

    def read_from(self, filename: str) -> None:
        """Read XDATCAR data from a VASP XDATCAR file.

        Args:
            filename (str): The XDATCAR file to read.

        Returns:
            None

        """
        self.poscar[0].read_from(filename)
        with open(filename) as f:
            lines = f.read()
        frame_header = re.compile("\nDirect configuration=\s+\d+\n")
        frame_coordinates = [
            frame.split("\n") for frame in frame_header.split(lines)[2:]
        ]
        for frame in frame_coordinates:
            self.poscar.append(copy.deepcopy(self.poscar[0]))
            self.poscar[-1].coordinates = np.array(
                [
                    [float(e) for e in frame.pop(0).split()[0:3]]
                    for i in range(sum(self.poscar[0].atom_numbers))
                ]
            )
