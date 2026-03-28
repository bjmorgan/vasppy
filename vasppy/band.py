import warnings


def handle_occupancy(
    occupancy: float,
    negative_occupancies: str = "warn",
) -> float:
    """Handle a potentially negative band occupancy.

    Args:
        occupancy: The band occupancy value to check.
        negative_occupancies: Policy for negative occupancies. One of:

            - ``'warn'``: issue a :class:`UserWarning` (default).
            - ``'raise'``: raise a :class:`ValueError`.
            - ``'ignore'``: return the value unchanged.
            - ``'zero'``: set the occupancy to ``0.0``.

    Returns:
        The (possibly corrected) occupancy value.

    Raises:
        ValueError: If ``negative_occupancies`` is not a recognised keyword.
        ValueError: If ``negative_occupancies`` is ``'raise'`` and the
            occupancy is negative.
    """
    valid_negative_occupancies = ["warn", "raise", "ignore", "zero"]
    if negative_occupancies not in valid_negative_occupancies:
        raise ValueError(
            f"valid options for negative_occupancies are {valid_negative_occupancies}"
        )
    if occupancy < 0:
        if negative_occupancies == "warn":
            warnings.warn(
                "One or more occupancies in your PROCAR file are negative.",
                stacklevel=2,
            )
        elif negative_occupancies == "raise":
            raise ValueError(
                "One or more occupancies in your PROCAR file are negative."
            )
        elif negative_occupancies == "ignore":
            pass
        elif negative_occupancies == "zero":
            occupancy = 0.0
    return occupancy


class Band:
    """Represents a single Kohn-Sham band.

    Attributes:
        index: Band index.
        energy: Band energy (eV).
        occupancy: Band occupancy.
    """

    def __init__(
        self,
        index: int,
        energy: float,
        occupancy: float,
        negative_occupancies: str = "warn",
    ) -> None:
        """Initialise a Band instance.

        Args:
            index: Band index.
            energy: Band energy (eV).
            occupancy: Band occupancy.
            negative_occupancies: Policy for negative occupancies passed to
                :func:`handle_occupancy`. Defaults to ``'warn'``.
        """
        self.index = index
        self.energy = energy
        self.occupancy = handle_occupancy(
            occupancy, negative_occupancies=negative_occupancies
        )

    def __eq__(self, other: object) -> bool:
        """Test equality between two Band instances.

        Args:
            other: Object to compare against.

        Returns:
            ``True`` if index, energy and occupancy are all equal.
        """
        if not isinstance(other, Band):
            return NotImplemented
        return (
            self.index == other.index
            and self.energy == other.energy
            and self.occupancy == other.occupancy
        )
