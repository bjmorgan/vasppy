import hashlib
from monty.io import zopen  # type: ignore
from pathlib import Path
import os
import yaml
from contextlib import contextmanager
from pymatgen.core import Structure
from typing import Optional, List
import numpy as np

@contextmanager
def cd( path ):
    old_dir = os.getcwd()
    os.chdir( path )
    try:
        yield
    finally:
        os.chdir( old_dir )

def md5sum( string ):
    """
    Generate the md5 checksum for a string

    Args:
        string (Str): The string to be checksummed.

    Returns:
        (Str): The hex checksum.
    """
    h = hashlib.new( 'md5' )
    h.update( string.encode( 'utf-8' ) )
    return h.hexdigest()

def file_md5( filename ):
    """
    Generate the md5 checksum for a file

    Args:
        filename (Str): The file to be checksummed.

    Returns:
        (Str): The hex checksum

    Notes:
        If the file is gzipped, the md5 checksum returned is
        for the uncompressed ASCII file.
    """
    with zopen( filename, 'r' ) as f:
        file_string = f.read()
    try: # attempt to decode byte object
        file_string = file_string.decode()
    except AttributeError:
        pass
    return( md5sum( file_string ) )

def match_filename( filename ):
    """
    Checks whether a file exists, either as named, or as a a gzippped file (filename.gz)

    Args:
        (Str): The root filename.

    Returns:
        (Str|None): if the file exists (either as the root filename, or gzipped), the return
            value will be the actual filename. If no matching filename is found the return
            value is set to None
    """
    f = next( ( '{}{}'.format( filename, extension ) for extension in [ '', '.gz' ]
        if Path( '{}{}'.format( filename, extension ) ).is_file() ), None ) 
    return f
   
def validate_checksum( filename, md5sum ):
    """
    Compares the md5 checksum of a file with an expected value.
    If the calculated and expected checksum values are not equal, 
    ValueError is raised.
    If the filename `foo` is not found, will try to read a gzipped file named
    `foo.gz`. In this case, the checksum is calculated for the unzipped file.

    Args:
        filename (str): Path for the file to be checksummed.
        md5sum (str):  The expected hex checksum.

    Returns:
        None
    """
    filename = match_filename( filename )
    md5_hash = file_md5( filename=filename )
    if md5_hash != md5sum:
        raise ValueError('md5 checksums are inconsistent: {}'.format( filename ))
        
def dr_ij(structure: Structure,
      indices_i: Optional[List[int]] = None,
      indices_j: Optional[List[int]] = None,
      self_reference: bool = False) -> np.ndarray:
    """
    Calculate all i-j interatomic distances for a single pymatgen Structure.

    Args:
        structure (:obj:`pymatgen.Structure`): A pymatgen Structure
        indices_i (:obj:`list(int)`, optional): List of indices for species i.
            Optional, default is `None`.
            If `indices_i` is not specified, then distances will be calculated between 
            all pairs of atoms.
        indices_j (:obj:`list(int)`, optional): List of indices for species j.
            Optional, default is `None`.
            If `indices_j` is not specified, then `indices_j` will be set equal
            to `indices_i`.
        self_reference (bool, optional): If computing distances for i==j, whether
            to include the i==j dr=0 terms. Optional, default is `False`.
            
    Returns:
        np.array: N_i x N_j numpy array of i-j minimum image distances.

    """
    if indices_i is None:
        indices_i = list(range(len(structure)))
    if indices_j is None:
        indices_j = indices_i
    lattice = structure.lattice
    i_frac_coords = structure.frac_coords[indices_i]
    j_frac_coords = structure.frac_coords[indices_j]
    dr_ij = lattice.get_all_distances(i_frac_coords, j_frac_coords)
    # If indices_i and indices_j contain common elements AND self_reference == False
    # then we want to mask dr_ij to remove the i==j dr=0 terms
    if (np.intersect1d(indices_i, indices_j).size > 0) and not self_reference:
        mask = np.ones_like(dr_ij, dtype=bool)
        for i_loc, i in enumerate(indices_i):
            for j_loc, j in enumerate(indices_j):
                if i == j:
                    mask[i_loc, j_loc] = 0
        to_return = dr_ij[mask].reshape(len(indices_i), -1)
    else:
        to_return = dr_ij
    return to_return
