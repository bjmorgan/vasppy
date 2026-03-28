"""Class representing a single atom with a label and coordinates."""

import numpy as np


class Atom:
    """Represents a single atom with a label and positional coordinates."""

    def __init__(self, label: str, r: np.ndarray) -> None:
        """Initialise an Atom instance.

        Args:
            label: A label for this atom (e.g. the element symbol).
            r: The atom coordinates as a NumPy array. Currently assumed to be
                fractional coordinates.
        """
        self.label = label
        self.r = r
