import hashlib
from monty.io import zopen
from pathlib import Path
import os
import yaml
from contextlib import contextmanager

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
