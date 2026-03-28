# Summary class and helper methods
# Used for summarising VASP calculations as YAML

from collections.abc import Callable

from pymatgen.io.vasp.outputs import Vasprun  # type: ignore
from pymatgen.analysis.transition_state import NEBAnalysis  # type: ignore
from vasppy.vaspmeta import VASPMeta
from vasppy.outcar import (
    final_energy_from_outcar,
    vasp_version_from_outcar,
    potcar_eatom_list_from_outcar,
)
from vasppy.data.potcar_data import potcar_md5sum_data
from vasppy.utils import file_md5, md5sum, match_filename, cd
from xml.etree import ElementTree as ET
import glob
import re
import warnings

import yaml # type: ignore

potcar_sets = [
    "PBE",
    "PBE_52",
    "PBE_54",
    "PBE_54r",
    "LDA_54r",
    "LDA",
    "LDA_52",
    "LDA_54",
    "GGA",
    "USPP_GGA",
    "USPP_LDA",
]


def load_vasp_summary(filename: str) -> dict[str, dict]:
    """Read a ``vasp_summary.yaml`` file and return a dictionary of documents.

    Each YAML document in the file becomes a sub-dictionary keyed by its
    ``title`` value.

    Example:
        The file::

            ---
            title: foo
            data: foo_data
            ---
            title: bar
            data: bar_data

        is converted to::

            {'foo': {'title': 'foo', 'data': 'foo_data'},
             'bar': {'title': 'bar', 'data': 'bar_data'}}

    Args:
        filename: File path for the ``vasp_summary.yaml`` file.

    Returns:
        Dictionary mapping each document's title to its full document
        dictionary.
    """
    with open(filename, "r") as stream:
        docs = yaml.load_all(stream, Loader=yaml.SafeLoader)
        data = {d["title"]: d for d in docs}
    return data


def potcar_spec(
    filename: str,
    return_hashes: bool = False,
) -> tuple[list[str], list[str]]:
    """Return pseudopotential names and dataset labels from a POTCAR file.

    Parses a POTCAR file, splitting it into individual pseudopotential
    blocks and matching each against known md5 checksums. Returns a pair
    of aligned lists so that duplicate species are preserved.

    Args:
        filename: The name of the POTCAR file to process.
        return_hashes: If True the second list contains the md5 hashes of
            the component pseudopotential strings instead of dataset labels.

    Returns:
        A ``(names, values)`` tuple of two lists. *names* contains the
        pseudopotential labels (e.g. ``['Fe_pv', 'O']``) and *values*
        contains either the dataset labels (e.g. ``['PBE_54', 'PBE_54']``)
        or, when *return_hashes* is True, the md5 hashes.

    Raises:
        ValueError: If any pseudopotential block cannot be matched to a
            known md5 hash.
    """
    with open(filename, "r") as f:
        potcars = [s for s in re.split("(End of Dataset\n)", f.read()) if s]
    potcar_md5sums = [
        md5sum("".join(pair))
        for pair in zip(potcars[::2], potcars[1::2])
    ]
    names: list[str] = []
    values: list[str] = []
    for this_md5sum in potcar_md5sums:
        for ps in potcar_sets:
            for p, p_md5sum in potcar_md5sum_data[ps].items():
                if this_md5sum == p_md5sum:
                    names.append(p)
                    if return_hashes:
                        values.append(this_md5sum)
                    else:
                        values.append(ps)
    if len(names) != len(potcar_md5sums):
        raise ValueError("One or more POTCARs did not have matching md5 hashes")
    return names, values


def find_vasp_calculations() -> list[str]:
    """Return a list of all subdirectories that contain a ``vasprun.xml`` file.

    Searches recursively from the current directory for ``vasprun.xml`` and
    ``vasprun.xml.gz`` files.

    Returns:
        List of directory paths (relative, with leading ``./``) that contain
        VASP calculation output.
    """
    dir_list = [
        "./" + re.sub(r"vasprun\.xml", "", path)
        for path in glob.iglob("**/vasprun.xml", recursive=True)
    ]
    gz_dir_list = [
        "./" + re.sub(r"vasprun\.xml\.gz", "", path)
        for path in glob.iglob("**/vasprun.xml.gz", recursive=True)
    ]
    return dir_list + gz_dir_list


class Summary:
    """Summarise a VASP calculation directory as a structured YAML document.

    Reads a ``vaspmeta.yaml`` file and a ``vasprun.xml`` (or
    ``vasprun.xml.gz``) from a given directory and exposes methods to
    print individual fields in YAML format.  The :meth:`output` method
    drives the overall summary output.

    Attributes:
        directory: Path to the VASP calculation directory.
        meta: A :obj:`VASPMeta` instance parsed from ``vaspmeta.yaml``.
        vasprun: A pymatgen :obj:`Vasprun` instance, or ``None`` if the
            ``vasprun.xml`` could not be parsed.
        vasprun_filename: The matched filename for the vasprun file.
        print_methods: Mapping from flag names to bound print methods.
        supported_flags: Mapping from flag names to human-readable labels.
    """

    supported_flags = {
        "title": "Title",
        "description": "Description",
        "notes": "Notes",
        "type": "Type",
        "status": "Status",
        "stoichiometry": "Stoichiometry",
        "potcar": "POTCAR",
        "eatom": "POTCAR EATOM values",
        "plus_u": "Dudarev +U parameters",
        "energy": "Energy",
        "lreal": "LREAL",
        "k-points": "k-points",
        "functional": "functional",
        "encut": "encut",
        "ediffg": "ediffg",
        "ibrion": "ibrion",
        "converged": "converged",
        "md5": "md5",
        "directory": "directory",
        "vbm": "Vasprun valence band maximum",
        "cbm": "Vasprun conduction band minimum",
        "track": "tracking for files",
        "version": "VASP executable version",
        "nelect": "NELECT",
    }

    def __init__(self, directory: str = ".") -> None:
        """Initialise a Summary object for a VASP calculation directory.

        Args:
            directory: Path to the directory containing ``vaspmeta.yaml``
                and ``vasprun.xml``. Default is the current directory.

        Raises:
            FileNotFoundError: If ``vaspmeta.yaml`` is not found in
                ``directory``.
            FileNotFoundError: If no ``vasprun.xml`` or ``vasprun.xml.gz``
                file is found in ``directory``.
            ValueError: If the set of supported flags does not match the
                set of registered print methods (internal consistency check).
        """
        self.directory = directory
        self.vasprun: Vasprun | None = None
        with cd(directory):
            try:
                self.meta = VASPMeta.from_file("vaspmeta.yaml")
            except FileNotFoundError as exc:
                raise FileNotFoundError(
                    f"vaspmeta.yaml not found in {directory}"
                ) from exc
            self.parse_vasprun()
        self.print_methods: dict[str, Callable[[], None]] = {
            "title": self.print_title,
            "description": self.print_description,
            "notes": self.print_notes,
            "type": self.print_type,
            "status": self.print_status,
            "stoichiometry": self.print_stoichiometry,
            "potcar": self.print_potcar,
            "eatom": self.print_eatom,
            "energy": self.print_energy,
            "k-points": self.print_kpoints,
            "functional": self.print_functional,
            "encut": self.print_encut,
            "plus_u": self.print_plus_u,
            "ediffg": self.print_ediffg,
            "ibrion": self.print_ibrion,
            "converged": self.print_converged,
            "version": self.print_version,
            "md5": self.print_vasprun_md5,
            "directory": self.print_directory,
            "lreal": self.print_lreal,
            "vbm": self.print_vbm,
            "cbm": self.print_cbm,
            "track": self.print_file_tracking,
            "nelect": self.print_nelect,
        }
        if set(self.print_methods.keys()) != set(self.supported_flags):
            raise ValueError(
                f"print_methods keys do not match supported_flags: "
                f"{set(self.print_methods.keys()) ^ set(self.supported_flags.keys())}"
            )

    def parse_vasprun(self) -> None:
        """Read ``vasprun.xml`` as a pymatgen Vasprun object.

        Sets ``self.vasprun_filename`` and ``self.vasprun``. If the
        ``vasprun.xml`` is malformed, ``self.vasprun`` is set to ``None``
        rather than raising.

        Raises:
            FileNotFoundError: If no ``vasprun.xml`` or ``vasprun.xml.gz``
                file can be found.
        """
        self.vasprun_filename = match_filename("vasprun.xml")
        if not self.vasprun_filename:
            raise FileNotFoundError("Could not find vasprun.xml or vasprun.xml.gz file")
        try:
            self.vasprun = Vasprun(
                self.vasprun_filename, parse_potcar_file=False, parse_dos=False
            )
        except ET.ParseError:
            warnings.warn(
                f"Could not parse {self.vasprun_filename} in {self.directory}; "
                "summary output will be incomplete",
                stacklevel=2,
            )
            self.vasprun = None

    @property
    def _vasprun(self) -> Vasprun:
        """Return the parsed Vasprun object, raising if unavailable.

        Raises:
            RuntimeError: If ``vasprun.xml`` could not be parsed.
        """
        if self.vasprun is None:
            raise RuntimeError("vasprun.xml could not be parsed")
        return self.vasprun

    @property
    def stoich(self) -> dict:
        """Elemental stoichiometry of the final structure.

        Returns:
            Dictionary mapping element symbol strings to float amounts.

        Raises:
            RuntimeError: If ``vasprun.xml`` could not be parsed.
        """
        return self._vasprun.final_structure.composition.get_el_amt_dict()

    @property
    def functional(self) -> str:
        """String description of the calculation functional.

        Returns:
            String describing the DFT functional used.

        Raises:
            RuntimeError: If ``vasprun.xml`` could not be parsed.
        """
        return self._vasprun.run_type

    def potcars_are_pbe(self) -> bool:
        """Check whether all POTCARs are PBE type.

        Returns:
            True if all POTCAR symbols contain ``'PBE'``, otherwise False.

        Raises:
            RuntimeError: If ``vasprun.xml`` could not be parsed.
        """
        return all("PBE" in s for s in self._vasprun.potcar_symbols)

    def output(self, to_print: list[str]) -> None:
        """Write summary fields to stdout in YAML document format.

        If the vasprun is unavailable (None), only ``title``, ``type``,
        and ``status`` are printed.

        Args:
            to_print: List of flag names to print in order. Valid flag
                names are the keys of :attr:`supported_flags`.
        """
        if not self.vasprun:
            to_print = ["title", "type", "status"]
        print("---")
        for p in to_print:
            self.print_methods[p]()
        print("", flush=True)

    def print_type(self) -> None:
        """Print the calculation type if set."""
        if self.meta.type:
            print(f"type: {self.meta.type}")

    def print_title(self) -> None:
        """Print the calculation title."""
        print(f"title: {self.meta.title}")

    def print_description(self) -> None:
        """Print the calculation description."""
        print(f"description: {self.meta.description.strip()}")

    def print_notes(self) -> None:
        """Print notes, or a YAML null marker if notes are not set."""
        if self.meta.notes:
            print(f"notes: {self.meta.notes.strip()}")
        else:
            print("notes: ~")

    def print_status(self) -> None:
        """Print the calculation status."""
        print(f"status: {self.meta.status}")

    def print_lreal(self) -> None:
        """Print the LREAL INCAR parameter."""
        print(f"lreal: {self._vasprun.parameters['LREAL']}")

    def print_stoichiometry(self) -> None:
        """Print the elemental stoichiometry."""
        print("stoichiometry:")
        for element in self.stoich:
            print(f"    - {element}: {int(self.stoich[element])}")

    def print_potcar(self) -> None:
        """Print the POTCAR species and symbols."""
        print("potcar:")
        for e, p in zip(self.stoich, self._vasprun.potcar_symbols):
            print(f"    - {e}: {p}")

    def print_energy(self) -> None:
        """Print the final energy of the calculation.

        For NEB calculations, prints image energies via
        :meth:`print_neb_energy`. For standard calculations, prints the
        ``vasprun.final_energy`` value.

        Raises:
            ValueError: If ``meta.type`` is set to an unsupported value.
        """
        if not self.meta.type:
            print(f"energy: {self._vasprun.final_energy}")
        elif self.meta.type == "neb":
            self.print_neb_energy()
        else:
            raise ValueError(f"VASPMeta type not supported: {self.meta.type}")

    def print_neb_energy(self) -> None:
        """Print the NEB image energies."""
        image_00_energy = final_energy_from_outcar("00/OUTCAR")
        print(f"reference energy: {image_00_energy} eV")
        neb = NEBAnalysis.from_dir(".")
        print("neb image energies:")
        for i, e in enumerate(neb.energies):
            print(f"    - {i:02d}: {e:10.6f} eV")

    def print_version(self) -> None:
        """Print the VASP executable version string."""
        version_string = vasp_version_from_outcar(
            f"{self.directory}/OUTCAR"
        ).split()[0]
        print(f"version: {version_string}")

    def print_eatom(self) -> None:
        """Print the EATOM values from the OUTCAR for each species."""
        print("eatom:")
        for e, eatom in zip(
            self.stoich,
            potcar_eatom_list_from_outcar(f"{self.directory}/OUTCAR")
        ):
            print(f"    - {e}: {eatom} eV")

    def print_kpoints(self) -> None:
        """Print the k-point scheme and grid."""
        print("k-points:")
        print(f"    scheme: {self._vasprun.kpoints.style}")
        print(f"    grid: {' '.join(str(k) for k in self._vasprun.kpoints.kpts[0])}")

    def print_functional(self) -> None:
        """Print the DFT functional."""
        print(f"functional: {self.functional}")

    def print_ibrion(self) -> None:
        """Print the IBRION INCAR parameter."""
        print(f"ibrion: {self._vasprun.incar['IBRION']}")

    def print_ediffg(self) -> None:
        """Print the EDIFFG INCAR parameter."""
        print(f"ediffg: {self._vasprun.incar['EDIFFG']}")

    def print_encut(self) -> None:
        """Print the ENCUT (or ENMAX) INCAR parameter."""
        if "ENCUT" in self._vasprun.incar:
            print(f"encut: {self._vasprun.incar['ENCUT']}")
        elif "ENMAX" in self._vasprun.incar:
            print(f"encut: {self._vasprun.incar['ENMAX']}")

    def print_converged(self) -> None:
        """Print the convergence status."""
        print(f"converged: {self._vasprun.converged}")

    def print_vasprun_md5(self) -> None:
        """Print the md5 checksum of the vasprun.xml file."""
        print(
            f"vasprun md5: {file_md5(f'{self.directory}/{self.vasprun_filename}')}"
        )

    def print_file_tracking(self) -> None:
        """Print tracking information for any files listed in the metadata."""
        if self.meta.track:
            print("file tracking:")
            for f, new_filename in self.meta.track.items():
                print(f"    {f}:")
                if not new_filename:
                    new_filename = f
                print(f"        filename: {new_filename}")
                filename = match_filename(self.directory + f)
                if filename:
                    md5 = file_md5(filename)
                else:
                    md5 = "null"
                print(f"        md5: {md5}")

    def print_directory(self) -> None:
        """Print the calculation directory path."""
        print(f"directory: {self.directory}")

    def print_plus_u(self) -> None:
        """Print Dudarev DFT+U parameters if present."""
        if "LDAUU" in self._vasprun.incar:
            lqn = {0: "s", 1: "p", 2: "d", 3: "f"}
            ldauu = self._vasprun.incar["LDAUU"]
            ldauj = self._vasprun.incar["LDAUJ"]
            ldaul = self._vasprun.incar["LDAUL"]
            if any(v != 0 for v in ldauu):
                print("ldau:")
                for e, u, j, l in zip(self.stoich, ldauu, ldauj, ldaul):
                    if u != 0:
                        print(f"    - {e}: {lqn[l]} {u} {j}")

    def print_cbm(self) -> None:
        """Print the conduction band minimum from the vasprun."""
        print(f"cbm: {self._vasprun.eigenvalue_band_properties[1]}")

    def print_vbm(self) -> None:
        """Print the valence band maximum from the vasprun."""
        print(f"vbm: {self._vasprun.eigenvalue_band_properties[2]}")

    def print_nelect(self) -> None:
        """Print the NELECT INCAR parameter."""
        print(f"nelect: {self._vasprun.parameters['NELECT']}")
