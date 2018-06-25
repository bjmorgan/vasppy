import hashlib

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
    """
    with open( filename, 'r' ) as f:
        file_string = f.read()
    return( md5sum( file_string ) )
