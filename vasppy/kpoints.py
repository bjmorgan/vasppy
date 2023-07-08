import numpy as np  # type: ignore
from typing import Optional


class AutoKPoints:
    """Class for automatic k-point generation data in KPOINTS."""

    def __init__(
        self,
        title: str,
        subdivisions: np.ndarray,
        grid_centering: Optional[str] = "G",
        shift: Optional[np.ndarray] = np.array([0.0, 0.0, 0.0]),
    ) -> None:
        """Initialise an AutoKPoints object.

        Args:
            title (str): The first line of the file, treated as a comment by VASP.
            subdivisions: (np.ndarray(int, int, int)):
                Numbers of subdivisions along each reciprocal lattice vector.
            grid_centering (str, optional):
                Specify gamma-centered (G) or the original Monkhorst-Pack scheme (MP).
                Default is 'G'.
            shift: (np.ndarray(float, float, float), optional):
                Optional shift of the mesh (s_1, s_2, s_3). Default is ( [ 0., 0., 0. ] ).

        Returns:
            None

        Raises:
            ValueError: If an unrecognised grid-centering option is passed in.

        """
        accepted_grid_centerings = ["G", "MP"]
        if grid_centering not in accepted_grid_centerings:
            raise ValueError
        self.title = title
        self.grid_centering = grid_centering
        self.subdivisions = subdivisions
        self.shift = shift
