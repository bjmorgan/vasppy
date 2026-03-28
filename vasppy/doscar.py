import numpy as np
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.axes import Axes  # type: ignore
from matplotlib.figure import Figure  # type: ignore
import matplotlib._color_data as mcd  # type: ignore
from typing import ClassVar, Literal
from collections.abc import Iterable, Sequence


TABLEAU_GREY: str = "#bab0ac"


def pdos_column_names(lmax: int, ispin: int) -> list[str]:
    """Return column names for a projected DOS dataframe.

    Args:
        lmax: Maximum angular momentum quantum number (2 for d, 3 for f).
        ispin: Number of spin channels (1 or 2).

    Returns:
        List of column name strings beginning with ``'energy'``.

    Raises:
        ValueError: If ``lmax`` is not 2 or 3.
    """
    if lmax == 2:
        names = ["s", "p_y", "p_z", "p_x", "d_xy", "d_yz", "d_z2-r2", "d_xz", "d_x2-y2"]
    elif lmax == 3:
        names = [
            "s",
            "p_y",
            "p_z",
            "p_x",
            "d_xy",
            "d_yz",
            "d_z2-r2",
            "d_xz",
            "d_x2-y2",
            "f_y(3x2-y2)",
            "f_xyz",
            "f_yz2",
            "f_z3",
            "f_xz2",
            "f_z(x2-y2)",
            "f_x(x2-3y2)",
        ]
    else:
        raise ValueError("lmax value not supported")
    if ispin == 2:
        all_names: list[str] = [f"{n}_{s}" for n in names for s in ("up", "down")]
    else:
        all_names = names
    all_names.insert(0, "energy")
    return all_names


class Doscar:
    """Contains all the data in a VASP DOSCAR file and methods for manipulating it."""

    number_of_header_lines: int = 6
    _spin_map: ClassVar[dict[str, list[int]]] = {"up": [0], "down": [1], "both": [0, 1]}
    _valid_m_values: ClassVar[dict[str, list[str]]] = {
        "s": [],
        "p": ["y", "z", "x"],
        "d": ["xy", "yz", "z2-r2", "xz", "x2-y2"],
        "f": ["y(3x2-y2)", "xyz", "yz2", "z3", "xz2", "z(x2-y2)", "x(x2-3y2)"],
    }
    _l_offsets: ClassVar[dict[str, int]] = {"s": 0, "p": 1, "d": 4, "f": 9}
    _l_widths: ClassVar[dict[str, int]] = {"s": 1, "p": 3, "d": 5, "f": 7}

    def __init__(
        self,
        filename: str,
        ispin: int = 2,
        lmax: int = 2,
        lorbit: int = 11,
        spin_orbit_coupling: bool = False,
        read_pdos: bool = True,
        species: list[str] | None = None,
    ) -> None:
        """Create a Doscar object from a VASP DOSCAR file.

        Args:
            filename: Filename of the VASP DOSCAR file to read.
            ispin: ISPIN flag. Set to 1 for non-spin-polarised or 2 for
                spin-polarised calculations. Default is 2.
            lmax: Maximum l angular momentum (d=2, f=3). Default is 2.
            lorbit: The VASP LORBIT flag. Default is 11.
            spin_orbit_coupling: Spin-orbit coupling flag. Default is False.
            read_pdos: Set to True to read the atom-projected density of
                states. Default is True.
            species: List of atomic species strings, e.g.
                ``['Fe', 'Fe', 'O', 'O', 'O']``. Default is None.

        Raises:
            NotImplementedError: If ``spin_orbit_coupling`` is True.
        """
        self.filename = filename
        self.ispin = ispin
        self.lmax = lmax
        self.spin_orbit_coupling = spin_orbit_coupling
        if self.spin_orbit_coupling:
            raise NotImplementedError("Spin-orbit coupling is not yet implemented")
        self.lorbit = lorbit
        self.pdos: np.ndarray | None = None
        self.species = species
        self.read_header()
        self.read_total_dos()
        if read_pdos:
            self.read_projected_dos()

    @property
    def number_of_channels(self) -> int:
        """Number of lm-projection channels for the current lorbit and lmax.

        Returns:
            Integer channel count.

        Raises:
            NotImplementedError: If lorbit is not 11.
        """
        if self.lorbit == 11:
            return {2: 9, 3: 16}[self.lmax]
        raise NotImplementedError

    def read_header(self) -> None:
        """Read the six-line header from the DOSCAR file.

        Populates ``self.header`` and then calls :meth:`process_header`.
        """
        self.header: list[str] = []
        with open(self.filename, "r") as file_in:
            for _i in range(Doscar.number_of_header_lines):
                self.header.append(file_in.readline())
        self.process_header()

    def process_header(self) -> None:
        """Parse key values from the header lines.

        Sets ``self.number_of_atoms``, ``self.number_of_data_points``, and
        ``self.efermi`` from the raw header strings.
        """
        self.number_of_atoms = int(self.header[0].split()[0])
        self.number_of_data_points = int(self.header[5].split()[2])
        self.efermi = float(self.header[5].split()[3])

    def read_total_dos(self) -> None:
        """Read the total DOS block from the DOSCAR file.

        Populates ``self.energy`` and ``self.tdos``.
        """
        start_to_read: int = Doscar.number_of_header_lines
        df: pd.DataFrame = pd.read_csv(
            self.filename,
            skiprows=start_to_read,
            nrows=self.number_of_data_points,
            sep=r'\s+',
            names=["energy", "up", "down", "int_up", "int_down"],
            index_col=False,
        )
        self.energy: np.ndarray = df.energy.values
        df = df.drop("energy", axis=1)
        self.tdos = df
    def read_atomic_dos_as_df(self, atom_number: int) -> pd.DataFrame:
        """Read the projected DOS for a single atom as a dataframe.

        Args:
            atom_number: 1-based index of the atom to read.

        Returns:
            Dataframe of projected DOS values (energy column dropped).

        Raises:
            ValueError: If ``atom_number`` is outside the valid range.
        """
        if not (0 < atom_number <= self.number_of_atoms):
            raise ValueError(
                f"atom_number must be between 1 and {self.number_of_atoms}, got {atom_number}"
            )
        start_to_read = Doscar.number_of_header_lines + atom_number * (
            self.number_of_data_points + 1
        )
        df = pd.read_csv(
            self.filename,
            skiprows=start_to_read,
            nrows=self.number_of_data_points,
            sep=r'\s+',
            names=pdos_column_names(lmax=self.lmax, ispin=self.ispin),
            index_col=False,
        )
        return df.drop("energy", axis=1)

    def read_projected_dos(self) -> None:
        """Read the projected density of states data.

        Populates ``self.pdos`` as a 4D numpy array with dimensions
        ``[atom_no, energy_value, lm-projection, spin]``.
        """
        pdos_list = [self.read_atomic_dos_as_df(i + 1) for i in range(self.number_of_atoms)]
        if all(df.empty for df in pdos_list):
            raise ValueError("No projected DOS data found in file")
        self.pdos = np.vstack([np.array(df) for df in pdos_list]).reshape(
            self.number_of_atoms,
            self.number_of_data_points,
            self.number_of_channels,
            self.ispin,
        )

    def pdos_select(
        self,
        atoms: int | Sequence[int] | None = None,
        spin: str | None = None,
        l: str | None = None,
        m: list[str] | None = None,
    ) -> np.ndarray:
        """Return a subset of the projected density of states array.

        Args:
            atoms: Atom numbers to include in the selection. Atom numbers
                count from 0 (array index). Default selects all atoms.
            spin: Spin channel(s) to include. Accepted values are ``'up'``,
                ``'down'``, and ``'both'``. Default selects all available spin channels.
            l: Angular momentum to include. Accepted values are ``'s'``,
                ``'p'``, ``'d'``, and ``'f'``. Setting ``l`` without ``m``
                returns all projections for that angular momentum.
            m: One or more m-values. Requires ``l`` to be set. Valid values
                depend on ``l``:

                - ``l='s'``: no sub-selection (ignored).
                - ``l='p'``: one or more of ``['x', 'y', 'z']``.
                - ``l='d'``: one or more of
                  ``['xy', 'yz', 'z2-r2', 'xz', 'x2-y2']``.
                - ``l='f'``: one or more of
                  ``['y(3x2-y2)', 'xyz', 'yz2', 'z3', 'xz2',
                  'z(x2-y2)', 'x(x2-3y2)']``.

        Returns:
            4-dimensional numpy array of selected pDOS values with dimensions
            ``[atom_no, energy_value, lm-projection, spin]``.

        Raises:
            TypeError: If pDOS data has not been loaded.
            TypeError: If ``atoms`` is not a list.
            ValueError: If ``spin`` is not a recognised value.
            ValueError: If ``l`` is not a recognised value.
        """
        if self.pdos is None:
            raise TypeError("pdos data is not available; ensure read_pdos=True and the file contains pDOS data")
        atom_idx = self._resolve_atom_idx(atoms)
        spin_idx = self._resolve_spin_idx(spin)
        channel_idx = self._resolve_channel_idx(l, m)
        return self.pdos[atom_idx, :, :, :][:, :, channel_idx, :][:, :, :, spin_idx]

    def _resolve_atom_idx(self, atoms: int | Sequence[int] | None) -> list[int]:
        """Resolve the atoms argument to a list of atom indices.

        Args:
            atoms: A single atom index, a sequence of indices, or None to
                select all atoms.

        Returns:
            List of atom indices.
        """
        if atoms is None:
            atoms = range(self.number_of_atoms)
        elif isinstance(atoms, int):
            atoms = [atoms]
        return list(atoms)

    def _resolve_spin_idx(self, spin: str | None) -> list[int]:
        """Resolve the spin argument to a list of spin channel indices.

        Args:
            spin: One of ``'up'``, ``'down'``, ``'both'``, or None to select
                all available spin channels.

        Returns:
            List of spin channel indices.

        Raises:
            ValueError: If ``spin`` is specified for a non-spin-polarised calculation.
            ValueError: If ``spin`` is not a recognised value.
        """
        if spin is None:
            return list(range(self.ispin))
        if self.ispin == 1:
            raise ValueError("spin selection is not available for non-spin-polarised calculations")
        if spin not in self._spin_map:
            raise ValueError(f"'{spin}' is not a valid spin value; use 'up', 'down', or 'both'")
        return self._spin_map[spin]

    def _resolve_channel_idx(self, l: str | None, m: list[str] | None) -> list[int]:
        """Resolve the l and m arguments to a list of channel indices.

        Args:
            l: Angular momentum label (``'s'``, ``'p'``, ``'d'``, ``'f'``),
                or None to select all channels.
            m: List of m-value strings to sub-select within ``l``. Ignored
                for ``l='s'``. Pass None to select all projections for ``l``.

        Returns:
            List of channel indices.

        Raises:
            ValueError: If ``l`` is not a recognised angular momentum label.
        """
        if l is None:
            return list(range(self.number_of_channels))
        if l not in self._l_offsets:
            raise ValueError(f"'{l}' is not a valid angular momentum label; use 's', 'p', 'd', or 'f'")
        offset = self._l_offsets[l]
        if m is None or not self._valid_m_values[l]:
            return list(range(offset, offset + self._l_widths[l]))
        return [offset + i for i, v in enumerate(self._valid_m_values[l]) if v in m]

    def pdos_sum(
        self,
        atoms: int | Sequence[int] | None = None,
        spin: str | None = None,
        l: str | None = None,
        m: list[str] | None = None,
    ) -> np.ndarray:
        """Return the summed projected DOS over atoms, lm-projections, and spin.

        Args:
            atoms: Atom indices to include. Default is all atoms.
            spin: Spin channel(s) to include. Default is all spins.
            l: Angular momentum to include. Default is all.
            m: Sub-angular-momentum projections. Requires ``l`` to be set.

        Returns:
            1D numpy array of summed pDOS values indexed by energy.
        """
        return np.array(
            np.sum(self.pdos_select(atoms=atoms, spin=spin, l=l, m=m), axis=(0, 2, 3))
        )

    def plot_pdos(
        self,
        ax: Axes | None = None,
        to_plot: dict[str, list[str]] | None = None,
        colours: Iterable | None = None,
        plot_total_dos: bool | None = True,
        xrange: tuple[float, float] | None = None,
        ymax: float | None = None,
        scaling: dict[str, dict[str, float]] | None = None,
        split: bool = False,
        title: str | None = None,
        title_loc: Literal['left', 'center', 'right'] = "center",
        labels: bool = True,
        title_fontsize: int = 16,
        legend_pos: str = "outside",
    ) -> Figure | None:
        """Plot the projected density of states.

        Args:
            ax: Matplotlib axes object to plot on. If None, a new figure is
                created.
            to_plot: Dictionary mapping species labels to lists of orbital
                labels to plot, e.g. ``{'Fe': ['s', 'd'], 'O': ['p']}``.
                Default is to plot s, p, d (and f if ``lmax=3``) for every
                unique species in ``self.species``.
            colours: Iterable of colour values for successive traces. Defaults
                to the matplotlib Tableau colour set.
            plot_total_dos: Whether to plot the total DOS as a shaded region.
                Default is True.
            xrange: ``(x_min, x_max)`` energy window. Default plots all data.
            ymax: Maximum y-axis value. Default is auto-scaled.
            scaling: Nested dictionary of multiplicative scaling factors, e.g.
                ``{'Fe': {'d': 0.1}}``. Default applies no scaling.
            split: Not yet implemented. Reserved for future spin-split layouts.
            title: Plot title string. Default is no title.
            title_loc: Title alignment — ``'left'``, ``'center'``, or
                ``'right'``. Default is ``'center'``.
            labels: Whether to add axis labels. Default is True.
            title_fontsize: Font size for the title. Default is 16.
            legend_pos: Legend position. Use ``'outside'`` to place the legend
                to the right of the axes, or any valid Matplotlib legend
                ``loc`` string. Default is ``'outside'``.

        Returns:
            The :class:`matplotlib.figure.Figure` if a new figure was created,
            otherwise None.
        """
        if not ax:
            fig, ax = plt.subplots(1, 1, figsize=(8.0, 3.0))
        else:
            fig = None
        if not isinstance(ax, Axes):
            raise TypeError("ax must be a matplotlib Axes instance")
        if not colours:
            colours = mcd.TABLEAU_COLORS
        if not isinstance(colours, Iterable):
            raise TypeError("colours must be an iterable")
        color_iterator = (c for c in colours)

        if not scaling:
            scaling = {}

        if xrange:
            e_range = (self.energy >= xrange[0]) & (self.energy <= xrange[1])
        else:
            e_range = np.ones(self.energy.shape, dtype=bool)

        auto_ymax = 0.0

        if not to_plot:
            to_plot = {}
            if not isinstance(self.species, Iterable):
                raise TypeError("species must be set before calling plot_pdos without to_plot")
            for s in set(self.species):
                to_plot[s] = ["s", "p", "d"]
                if self.lmax == 3:
                    to_plot[s].append("f")

        for species in to_plot.keys():
            if not isinstance(self.species, Iterable):
                raise TypeError("species must be set before calling plot_pdos")
            index = [i for i, s in enumerate(self.species) if s == species]
            for state in to_plot[species]:
                if state not in ["s", "p", "d", "f"]:
                    raise ValueError(f"'{state}' is not a valid orbital label")
                color = next(color_iterator)
                label = f"{species} {state}"
                up_dos = self.pdos_sum(atoms=index, l=state, spin="up")[e_range]
                down_dos = self.pdos_sum(atoms=index, l=state, spin="down")[e_range]
                if species in scaling:
                    if state in scaling[species]:
                        up_dos *= scaling[species][state]
                        down_dos *= scaling[species][state]
                        label = rf"{species} {state} $\times${scaling[species][state]}"
                auto_ymax = max([auto_ymax, up_dos.max(), down_dos.max()])
                ax.plot(self.energy[e_range], up_dos, label=label, c=color)
                ax.plot(self.energy[e_range], down_dos * -1.0, c=color)
        if plot_total_dos:
            ax.fill_between(
                self.energy[e_range],
                self.tdos.up.values[e_range],
                self.tdos.down.values[e_range] * -1.0,
                facecolor=TABLEAU_GREY,
                alpha=0.2,
            )
            auto_ymax = max(
                [
                    auto_ymax,
                    self.tdos.up.values[e_range].max(),
                    self.tdos.down.values[e_range].max(),
                ]
            )

        if xrange:
            ax.set_xlim(xrange[0], xrange[1])

        if not ymax:
            ymax = 1.1 * auto_ymax
        ymax = float(ymax)
        ax.set_ylim(-ymax * 1.1, ymax * 1.1)
        if legend_pos == "outside":
            ax.legend(bbox_to_anchor=(1.01, 1.04), loc="upper left")
        else:
            ax.legend(loc=legend_pos)
        if labels:
            ax.set_xlabel("Energy [eV]")
        ax.axhline(y=0, c="lightgrey")
        ax.axes.grid(False, axis="y") # type: ignore

        ax.tick_params(
            axis="y",  # changes apply to the y-axis
            which="both",  # both major and minor ticks are affected
            left=False,  # ticks along the left edge are off
            right=False,  # ticks along the right edge are off
            labelleft=False,
        )  # labels along the left edge are off

        if title:
            ax.set_title(title, loc=title_loc, fontdict={"fontsize": title_fontsize})

        return fig
