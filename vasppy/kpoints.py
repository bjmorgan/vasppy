import numpy as np  # type: ignore
from typing import Optional

class AutoKPoints:
    """Class for automatic k-point generation data in KPOINTS."""

    def __init__(
        self,
        title: str,
        subdivisions: np.ndarray,
        grid_centering: Optional[str] = "G",
        shift: Optional[np.ndarray] = None,
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
        if shift is None:
            self.shift = np.array([0.0, 0.0, 0.0])
        else:
            self.shift = shift

def get_subdivisions_from_kspacing(kspacing: float,
                                             reciprocal_lattice_vectors: np.ndarray) -> tuple[int, int, int]:
    """Calculate subdivisions along reciprocal lattice vectors from the miniumum allowed distance between k-points (KSPACING).

    Args:
        kspacing (float): The minimum allowed distance between k-points.
        reciprocal_lattice_vectors (np.ndarray): The reciprocal lattice vectors. These can be retrieved from ASE as atoms.cell.reciprocal() or from Pymatgen as
        structure.lattice.reciprocal_lattice_crystallographic.matrix.

    Returns:
        tuple[int, int, int]: The subdivisions along each reciprocal lattice vector.
    """
    subdivisions = []
    for row in reciprocal_lattice_vectors:
        magnitude = np.linalg.norm(row)
        subdivision = np.max([1, np.ceil((magnitude * 2 * np.pi) / kspacing)])

        subdivisions.append(int(subdivision))

    return tuple(subdivisions)

def get_convergence_testing_kspacing(reciprocal_lattice_vectors: np.ndarray,
                                     kspacing_range: tuple[float, float]=(0.1, 0.8),
                                     step: float=0.02) -> tuple[float, ...]:
    """Generate a range of minimum allowed distances between k-points (KSPACING) for convergence testing. This function ensures that no two values of KSPACING
    are generated that correspond to the same k-point mesh.

    Args:
        reciprocal_lattice_vectors (np.ndarray): The reciprocal lattice vectors. These can be retrieved from ASE as atoms.cell.reciprocal() or from Pymatgen as
        structure.lattice.reciprocal_lattice_crystallographic.matrix.
        kspacing_range (tuple[float, float]): The minimum and maximum KSPACING values. Defaults to (0.1, 0.8).
        step (float): The interval between KSPACING values to be tested. Defaults to 0.02.

    Returns:
        tuple[float, ...]: A range of KSPACING values which all correspond to distinct k-point grids.
    """
    allowed_kspacing = []
    highest_total = 0.0
    kspacing_min, kspacing_max = kspacing_range
    for kspacing in np.arange(kspacing_min, kspacing_max + step, step):
        subdivisions = get_subdivisions_from_kspacing(kspacing, reciprocal_lattice_vectors)
        total = 1 / sum(subdivisions)

        if total > highest_total:
            allowed_kspacing.append(round(kspacing, 3))
            highest_total = total

    return tuple(allowed_kspacing)
