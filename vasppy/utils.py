"""Utility functions for file checksums, coordinate transforms, and context helpers."""

import hashlib
from monty.io import zopen  # type: ignore
from pathlib import Path
import os
from contextlib import contextmanager
from pymatgen.core import Structure  # type: ignore
import numpy as np
from typing import Generator


@contextmanager
def cd(path: str) -> Generator[None, None, None]:
    """Context manager that temporarily changes the working directory.

    Args:
        path: The directory to change into for the duration of the context.

    Yields:
        None
    """
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)


def md5sum(string: str) -> str:
    """Generate the md5 checksum for a string.

    Args:
        string: The string to be checksummed.

    Returns:
        The hex checksum.
    """
    h = hashlib.new("md5")
    h.update(string.encode("utf-8"))
    return h.hexdigest()


def file_md5(filename: str) -> str:
    """Generate the md5 checksum for a file.

    Args:
        filename: The file to be checksummed.

    Returns:
        The hex checksum.

    Notes:
        If the file is gzipped, the md5 checksum returned is for the
        uncompressed ASCII file.
    """
    with zopen(filename, "r") as f:
        file_string = f.read()

    if isinstance(file_string, bytes):
        file_string = file_string.decode()

    return md5sum(file_string)


def match_filename(filename: str) -> str | None:
    """Check whether a file exists, either as named or as a gzipped variant.

    Args:
        filename: The root filename to look for.

    Returns:
        The actual filename (with or without ``.gz`` extension) if found,
        otherwise None.
    """
    return next(
        (
            f"{filename}{extension}"
            for extension in ["", ".gz"]
            if Path(f"{filename}{extension}").is_file()
        ),
        None,
    )


def validate_checksum(filename: str, md5sum: str) -> None:
    """Compare the md5 checksum of a file with an expected value.

    If the filename *foo* is not found, will try to read a gzipped file named
    *foo.gz*. In that case the checksum is calculated for the unzipped content.

    Args:
        filename: Path of the file to be checksummed.
        md5sum: The expected hex checksum.

    Raises:
        FileNotFoundError: If the file is not found.
        ValueError: If the calculated and expected checksums do not match.
    """
    actual_filename = match_filename(filename)
    if actual_filename is None:
        raise FileNotFoundError(f"File not found: {filename}")
    md5_hash = file_md5(filename=actual_filename)
    if md5_hash != md5sum:
        raise ValueError(f"md5 checksums are inconsistent: {actual_filename}")


def dr_ij(
    structure: Structure,
    indices_i: list[int] | None = None,
    indices_j: list[int] | None = None,
    self_reference: bool = False,
) -> np.ndarray:
    """Calculate all i-j interatomic distances for a single pymatgen Structure.

    Args:
        structure: A pymatgen Structure.
        indices_i: List of indices for species i. If not specified, distances
            are calculated between all pairs of atoms. Defaults to None.
        indices_j: List of indices for species j. If not specified,
            *indices_j* is set equal to *indices_i*. Defaults to None.
        self_reference: If computing distances for i==j, whether to include
            the i==j dr=0 terms. Defaults to False.

    Returns:
        N_i x N_j NumPy array of i-j minimum image distances.
    """
    if indices_i is None:
        indices_i = list(range(len(structure)))
    if indices_j is None:
        indices_j = indices_i
    lattice = structure.lattice
    i_frac_coords = structure.frac_coords[indices_i]
    j_frac_coords = structure.frac_coords[indices_j]
    dr_ij_array = lattice.get_all_distances(i_frac_coords, j_frac_coords)
    # When computing self-referencing distances (indices_i == indices_j),
    # mask out the i==j dr=0 diagonal terms.
    if list(indices_i) == list(indices_j) and not self_reference:
        mask = np.array(indices_i)[:, None] != np.array(indices_j)[None, :]
        to_return = dr_ij_array[mask].reshape(len(indices_i), -1)
    else:
        to_return = dr_ij_array
    return np.asarray(to_return)
