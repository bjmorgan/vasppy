import hashlib
from monty.io import zopen

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
