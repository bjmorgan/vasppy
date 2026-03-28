"""Cell geometry utilities: angle/rotation helpers and the Cell class."""

import math
from typing import cast

import numpy as np


def angle(x: np.ndarray, y: np.ndarray) -> float:
    """Calculate the angle between two vectors, in degrees.

    Args:
        x: One vector.
        y: The other vector.

    Returns:
        The angle between x and y in degrees.
    """
    dot = np.dot(x, y)
    x_mod = np.linalg.norm(x)
    y_mod = np.linalg.norm(y)
    cos_angle = dot / (x_mod * y_mod)
    return float(np.degrees(np.arccos(cos_angle)))


def rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """Return the 3D rotation matrix for a counter-clockwise rotation about an axis.

    Args:
        axis: Length-3 array defining the axis of rotation.
        theta: Rotation angle in radians.

    Returns:
        The corresponding 3x3 rotation matrix.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2)
    b, c, d = -axis * math.sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


class Cell:
    """Represents a periodic unit cell defined by a 3x3 lattice matrix."""

    def __init__(self, matrix: np.ndarray) -> None:
        """Initialise a Cell object.

        Args:
            matrix: 3x3 NumPy array containing the cell matrix.

        Raises:
            ValueError: If *matrix* is not a NumPy ndarray.
            ValueError: If *matrix* does not have shape (3, 3).
        """
        if not isinstance(matrix, np.ndarray):
            raise ValueError("matrix must be a numpy ndarray.")
        if matrix.shape != (3, 3):
            raise ValueError(
                f"matrix must have shape (3, 3); got {matrix.shape}."
            )
        self.matrix = matrix  # 3 x 3 numpy array
        self.inv_matrix = np.linalg.inv(matrix)

    def dr(self, r1: np.ndarray, r2: np.ndarray, cutoff: float | None = None) -> float | None:
        """Calculate the distance between two fractional coordinates in the cell.

        Args:
            r1: Fractional coordinates for position 1.
            r2: Fractional coordinates for position 2.
            cutoff: If set, returns None for distances greater than the cutoff.
                Defaults to None.

        Returns:
            The distance between r1 and r2, or None if it exceeds *cutoff*.
        """
        delta_r_cartesian = (r1 - r2).dot(self.matrix)
        delta_r_squared = sum(delta_r_cartesian**2)
        if cutoff is not None:
            cutoff_squared = cutoff**2
            if delta_r_squared > cutoff_squared:
                return None
        return math.sqrt(delta_r_squared)

    def nearest_image(self, origin: np.ndarray, point: np.ndarray) -> np.ndarray:
        """Find the fractional coordinates of the nearest periodic image to a point of origin.

        Args:
            origin: Fractional coordinates of the point of origin.
            point: Fractional coordinates of the other point.

        Returns:
            The fractional coordinates of the nearest image of *point* to *origin*.
        """
        return cast(np.ndarray, origin + self.minimum_image(origin, point))

    def minimum_image(self, r1: np.ndarray, r2: np.ndarray) -> np.ndarray:
        """Find the minimum image vector from point r1 to point r2.

        Args:
            r1: Fractional coordinates of point r1.
            r2: Fractional coordinates of point r2.

        Returns:
            The fractional coordinate vector from r1 to the nearest image of r2.
        """
        delta_r = r2 - r1
        delta_r = np.array(
            [x - math.copysign(1.0, x) if abs(x) > 0.5 else x for x in delta_r]
        )
        return delta_r

    def minimum_image_dr(self, r1: np.ndarray, r2: np.ndarray, cutoff: float | None = None) -> float | None:
        """Calculate the shortest distance between two points, accounting for periodic boundary conditions.

        Args:
            r1: Fractional coordinates of point r1.
            r2: Fractional coordinates of point r2.
            cutoff: If set, return None if the minimum distance exceeds *cutoff*.
                Defaults to None.

        Returns:
            The minimum image distance between r1 and r2, or None if it exceeds *cutoff*.
        """
        delta_r_vector = self.minimum_image(r1, r2)
        return self.dr(np.zeros(3), delta_r_vector, cutoff)

    def lengths(self) -> np.ndarray:
        """The cell lengths.

        Returns:
            Array of cell lengths (a, b, c).
        """
        return np.array([math.sqrt(sum(row**2)) for row in self.matrix])

    def angles(self) -> list[float]:
        """The cell angles in degrees.

        Returns:
            List of cell angles [alpha, beta, gamma].
        """
        a, b, c = self.matrix
        return [angle(b, c), angle(a, c), angle(a, b)]

    def cartesian_to_fractional_coordinates(self, coordinates: np.ndarray) -> np.ndarray:
        """Convert a set of Cartesian coordinates to fractional coordinates in the cell.

        Args:
            coordinates: Array of shape (N, 3) containing Cartesian coordinates.

        Returns:
            Array of shape (N, 3) containing the corresponding fractional coordinates.
        """
        return cast(np.ndarray, coordinates.dot(self.inv_matrix))

    def fractional_to_cartesian_coordinates(self, coordinates: np.ndarray) -> np.ndarray:
        """Convert a set of fractional coordinates in the cell to Cartesian coordinates.

        Args:
            coordinates: Array of shape (N, 3) containing fractional coordinates.

        Returns:
            Array of shape (N, 3) containing the corresponding Cartesian coordinates.
        """
        return cast(np.ndarray, coordinates.dot(self.matrix))

    def inside_cell(self, r: np.ndarray) -> np.ndarray:
        """Return the equivalent point inside the cell for a fractional coordinate.

        Args:
            r: Fractional coordinates of a point (may lie outside the cell boundaries).

        Returns:
            Fractional coordinates of an equivalent point inside the cell boundaries.
        """
        centre = np.array([0.5, 0.5, 0.5])
        new_r = self.nearest_image(centre, r)
        return new_r

    def volume(self) -> float:
        """The cell volume.

        Returns:
            The scalar cell volume.
        """
        return float(np.dot(self.matrix[0], np.cross(self.matrix[1], self.matrix[2])))

    def unit_vectors(self) -> np.ndarray:
        """The unit vectors for the cell lattice vectors.

        Returns:
            Array of shape (3, 3) containing the unit vectors of each lattice vector.
        """
        return cast(np.ndarray, (self.matrix.transpose() / self.lengths()).transpose())

    def rotate(self, axis: np.ndarray, theta: float) -> None:
        """Rotate the cell in place about the given axis by theta radians.

        Args:
            axis: Length-3 array defining the axis of rotation.
            theta: Rotation angle in radians.
        """
        self.matrix = np.array(
            [np.dot(rotation_matrix(axis, theta), v) for v in self.matrix]
        )
