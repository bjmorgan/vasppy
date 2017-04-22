# Summary class and helper methods
# Used for summarising VASP calculations as YAML

from pymatgen.io.vasp.outputs import Vasprun
from vasppy.vaspmeta import VASPMeta
from contextlib import contextmanager
import os
import yaml
import hashlib
import glob
import re

def find_vasp_calculations():
    """
    Returns a list of all subdirectories that contain a vasprun.xml file

    Args:
        None

    Returns:
        (List): generator of all VASP calculation subdirectories.
    """
    dir_list = [ './' + re.sub( r'vasprun\.xml', '', path ) for path in glob.iglob( '**/vasprun.xml', recursive=True ) ]
    return dir_list

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

def vasprun_md5( filename ):
    """
    Generate the md5 checksum for a file

    Args:
        filename (Str): The file to be checksummed.

    Returns:
        (Str): The hex checksum
    """
    with open( filename, 'r' ) as f:
        vasprun = f.read()
    return( md5sum( vasprun ) )

def functional( str ):
    """
    Identifies the calculation functional as PBE or PBEsol based on the `GGA` INCAR tag.

    Args:
        str (Str): The GGA INCAR tag

    Returns:
        (Str): PBE | PBEsol

    Notes:
        Not tested whether this works for PBE calculations, where the GGA tag might be missing.
    """
    if str == 'PS':
        f = 'PBEsol'
    else:
        f = 'PBE'
    return f
   
@contextmanager
def cd( path ):
    old_dir = os.getcwd()
    os.chdir( path )
    try:
        yield
    finally:
        os.chdir( old_dir )

class Summary:

    def __init__( self, directory='.' ):
        self.directory = directory
        with cd( directory ):
            self.meta = VASPMeta.from_file( 'vaspmeta.yaml' )
            self.vasprun = Vasprun( 'vasprun.xml', parse_potcar_file=False )
        self.print_methods = { 'title': self.print_title,
                               'status': self.print_status,
                               'stoichiometry': self.print_stoichiometry,
                               'potcar': self.print_potcar,
                               'energy': self.print_energy,
                               'k-points': self.print_kpoints,
                               'functional': self.print_functional,
                               'encut': self.print_encut,
                               'plus_u': self.print_plus_u,
                               'ediffg': self.print_ediffg,
                               'ibrion': self.print_ibrion,
                               'converged': self.print_converged,
                               'md5': self.print_vasprun_md5,
                               'directory': self.print_directory }
            
    @property
    def stoich( self ):
        return self.vasprun.final_structure.composition.get_el_amt_dict()

    def output( self, to_print ):
        print( "---" )
        for p in to_print:
            self.print_methods[ p ]()
        print( '', flush=True )
        
    def print_title( self ):
        print( "title: {}".format( self.meta.title ) )
       
    def print_status( self ):
        print( "status: {}".format( self.meta.status ) )

    def print_stoichiometry( self ):
        print( "stoichiometry:" )
        for element in self.stoich:
            print( "    - {}: {}".format( element, int( self.stoich[ element ] ) ) )

    def print_potcar( self ):
        print( "potcar:" )
        for e, p in zip( self.stoich, self.vasprun.potcar_symbols ):
            print( "    - {}: {}".format( e, p ) )
        
    def print_energy( self ):
        print( "energy: {}".format( self.vasprun.final_energy ) )    
 
    def print_kpoints( self ):
        print( "k-points:" )
        print( "    scheme: {}".format( self.vasprun.kpoints.style ) )
        print( "    grid: {}".format( " ".join( str( k ) for k in self.vasprun.kpoints.kpts[0] ) ) )

    def print_functional( self ):
        print( "functional: {}".format( functional( self.vasprun.parameters['GGA'] ) ) )

    def print_ibrion( self ):
        print( "ibrion: {}".format( self.vasprun.incar['IBRION'] ) )

    def print_ediffg( self ):
        print( "ediffg: {}".format( self.vasprun.incar['EDIFFG'] ) )

    def print_encut( self ):
        print( "encut: {}".format( self.vasprun.incar['ENCUT'] ) )

    def print_converged( self ):
        print( "converged: {}".format( self.vasprun.converged ) )

    def print_vasprun_md5( self ):
        print( "vasprun md5: {}".format( vasprun_md5( "{}/vasprun.xml".format( self.directory ) ) ) )

    def print_directory( self ):
        print( "directory: {}".format( self.directory ) )

    def print_plus_u( self ):
        if 'LDAUU' in self.vasprun.incar:
            lqn = { 0: 's', 1: 'p', 2: 'd', 3: 'f' }
            ldauu = self.vasprun.incar[ 'LDAUU' ]
            ldauj = self.vasprun.incar[ 'LDAUJ' ]
            ldaul = self.vasprun.incar[ 'LDAUL' ]
            if any( v != 0 for v in ldauu ):
                print( 'ldau:' )
                for e, u, j, l in zip( self.stoich, ldauu, ldauj, ldaul ):
                    if u != 0:
                        print( "    - {}: {} {} {}".format( e, lqn[l], u, j ) )
