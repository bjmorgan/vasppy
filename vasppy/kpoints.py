"""Classes and functions for VASP KPOINTS generation and convergence testing."""

import numpy as np


class AutoKPoints:
    """Class for automatic k-point generation data in KPOINTS."""

    def __init__(
        self,
        title: str,
        subdivisions: np.ndarray,
        grid_centering: str | None = "G",
        shift: np.ndarray | None = None,
    ) -> None:
        """Initialise an AutoKPoints object.

        Args:
            title: The first line of the file, treated as a comment by VASP.
            subdivisions: Numbers of subdivisions along each reciprocal
                lattice vector.
            grid_centering: Specify gamma-centred (``'G'``) or the original
                Monkhorst-Pack scheme (``'MP'``). Default is ``'G'``.
            shift: Optional shift of the mesh ``(s_1, s_2, s_3)``.
                Default is ``[0., 0., 0.]``.

        Raises:
            ValueError: If an unrecognised grid-centering option is passed in.

        """
        accepted_grid_centerings = ["G", "MP"]
        if grid_centering not in accepted_grid_centerings:
            raise ValueError(
                f"Unrecognised grid-centering option: '{grid_centering}'. "
                f"Expected one of {accepted_grid_centerings}."
            )
        self.title = title
        self.grid_centering = grid_centering
        self.subdivisions = subdivisions
        if shift is None:
            self.shift = np.array([0.0, 0.0, 0.0])
        else:
            self.shift = shift

def get_subdivisions_from_kspacing(
    kspacing: float,
    reciprocal_lattice_vectors: np.ndarray,
) -> tuple[int, ...]:
    """Calculate subdivisions from the minimum allowed distance between k-points.

    Args:
        kspacing: The minimum allowed distance between k-points (KSPACING).
        reciprocal_lattice_vectors: The reciprocal lattice vectors.

    Returns:
        The subdivisions along each reciprocal lattice vector.
    """
    return tuple(
        int(np.max([1, np.ceil(np.linalg.norm(row) * 2 * np.pi / kspacing)]))
        for row in reciprocal_lattice_vectors
    )

def get_convergence_testing_kspacing(
    reciprocal_lattice_vectors: np.ndarray,
    kspacing_range: tuple[float, float] = (0.1, 0.8),
    step: float = 0.02,
) -> tuple[float, ...]:
    """Generate KSPACING values for convergence testing.

    Produces a range of KSPACING values ensuring that no two values
    correspond to the same k-point mesh.

    Args:
        reciprocal_lattice_vectors: The reciprocal lattice vectors.
        kspacing_range: The minimum and maximum KSPACING values.
            Defaults to ``(0.1, 0.8)``.
        step: The interval between KSPACING values to be tested.
            Defaults to ``0.02``.

    Returns:
        KSPACING values that each correspond to a distinct k-point grid.
    """
    allowed_kspacing = []
    highest_total = 0.0
    kspacing_min, kspacing_max = kspacing_range
    for kspacing in np.arange(kspacing_min, kspacing_max + step, step):
        kspacing_float = float(kspacing)
        subdivisions = get_subdivisions_from_kspacing(kspacing_float, reciprocal_lattice_vectors)
        total = 1 / sum(subdivisions)
        if total > highest_total:
            allowed_kspacing.append(round(kspacing_float, 3))
            highest_total = total
    return tuple(allowed_kspacing)
