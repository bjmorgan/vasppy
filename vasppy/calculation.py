"""Classes and functions for representing and manipulating VASP calculations."""

import yaml
import re
from collections import Counter
from collections.abc import Mapping


class Calculation:
    """Represents a single VASP calculation with a title, energy, and stoichiometry."""

    def __init__(self, title: str, energy: float, stoichiometry: Mapping[str, int | float]) -> None:
        """Initialise a Calculation object.

        Args:
            title: The title string for this calculation.
            energy: Final energy in eV.
            stoichiometry: A Mapping describing the calculation stoichiometry,
                e.g. ``{'Ti': 1, 'O': 2}``.
        """
        self.title = title
        self.energy = energy
        self.stoichiometry = Counter(stoichiometry)

    def __mul__(self, scaling: float) -> "Calculation":
        """Return a new Calculation with scaled energy and stoichiometry.

        Args:
            scaling: The scaling factor.

        Returns:
            The scaled Calculation.
        """
        return Calculation(
            title=self.title,
            energy=self.energy * scaling,
            stoichiometry=self.scale_stoichiometry(scaling),
        )

    def __truediv__(self, scaling: float) -> "Calculation":
        """Return a new Calculation divided by a scaling factor.

        Args:
            scaling: The scaling factor.

        Returns:
            The scaled Calculation.
        """
        return self * (1 / scaling)

    def scale_stoichiometry(self, scaling: float) -> dict[str, float]:
        """Return the stoichiometry scaled by *scaling*.

        Args:
            scaling: The scaling factor.

        Returns:
            The scaled stoichiometry as a dict of ``{label: count}`` pairs.
        """
        return {k: v * scaling for k, v in self.stoichiometry.items()}


def delta_E(
    reactants: list["Calculation"],
    products: list["Calculation"],
    check_balance: bool = True,
) -> float:
    """Calculate the change in energy for reactants --> products.

    Args:
        reactants: A list of Calculation objects representing the initial state.
        products: A list of Calculation objects representing the final state.
        check_balance: Check that the reaction stoichiometry is balanced.
            Defaults to True.

    Returns:
        The change in energy in eV.

    Raises:
        ValueError: If *check_balance* is True and the reaction is not balanced.
    """
    if check_balance:
        imbalance = delta_stoichiometry(reactants, products)
        if imbalance:
            raise ValueError(
                f"reaction is not balanced: {imbalance}"
            )
    return sum(r.energy for r in products) - sum(r.energy for r in reactants)


def delta_stoichiometry(
    reactants: list["Calculation"],
    products: list["Calculation"],
) -> dict[str, float]:
    """Calculate the change in stoichiometry for reactants --> products.

    Args:
        reactants: A list of Calculation objects representing the initial state.
        products: A list of Calculation objects representing the final state.

    Returns:
        A dict of non-zero stoichiometry changes, keyed by species label.
    """
    totals: Counter = Counter()
    for r in reactants:
        totals.update((r * -1.0).stoichiometry)
    for p in products:
        totals.update(p.stoichiometry)
    return {c: totals[c] for c in totals if totals[c] != 0}


def energy_string_to_float(string: str) -> float:
    """Convert an energy string such as ``'-1.2345 eV'`` to a float.

    Args:
        string: The string to convert.

    Returns:
        The numeric energy value.

    Raises:
        ValueError: If the string does not contain a valid energy value.
    """
    energy_re = re.compile(r"(-?\d+\.\d+)")
    match = energy_re.match(string)
    if match is None:
        raise ValueError(f"Could not parse energy from string: {string!r}")
    return float(match.group(0))


def import_calculations_from_file(
    filename: str,
    skip_incomplete_records: bool = False,
) -> dict[str, "Calculation"]:
    """Construct a dict of Calculation objects by reading a YAML file.

    Each YAML document should include ``title``, ``stoichiometry``, and
    ``energy`` fields, e.g.::

        title: my calculation
        stoichiometry:
            - A: 1
            - B: 2
        energy: -0.1234 eV

    Separate calculations should be distinct YAML documents, separated by ``---``.

    Args:
        filename: Path of the YAML file to read.
        skip_incomplete_records: Skip YAML documents that are missing one or
            more of the required keys. Defaults to False.

    Returns:
        A dict mapping each calculation title to its Calculation object.

    Raises:
        ValueError: If a document lacks a ``stoichiometry`` field (and
            *skip_incomplete_records* is False).
        ValueError: If more than one calculation shares the same title.
    """
    calcs: dict[str, Calculation] = {}
    with open(filename, "r") as stream:
        docs = yaml.load_all(stream, Loader=yaml.SafeLoader)
        for d in docs:
            if skip_incomplete_records:
                if (
                    ("title" not in d)
                    or ("stoichiometry" not in d)
                    or ("energy" not in d)
                ):
                    continue
            if "stoichiometry" in d:
                stoichiometry: Counter = Counter()
                for s in d["stoichiometry"]:
                    stoichiometry.update(s)
            else:
                raise ValueError(f'stoichiometry not found for "{d["title"]}"')
            if d["title"] in calcs:
                raise ValueError(
                    f'More than one calculation has the same title: {d["title"]}'
                )
            calcs[d["title"]] = Calculation(
                title=d["title"],
                stoichiometry=stoichiometry,
                energy=energy_string_to_float(d["energy"]),
            )
    return calcs
