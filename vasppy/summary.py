# Summary class and helper methods
# Used for summarising VASP calculations as YAML

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.transition_state import NEBAnalysis
from vasppy.vaspmeta import VASPMeta
from vasppy.outcar import final_energy_from_outcar, vasp_version_from_outcar, potcar_eatom_list_from_outcar
from vasppy.data.potcar_data import potcar_md5sum_data, potcar_nelect
from vasppy.utils import file_md5, md5sum, match_filename, cd
from xml.etree import ElementTree as ET 
import sys
import yaml
import glob
import re

potcar_sets = [ 'PBE', 'PBE_52', 'PBE_54', 'LDA', 'LDA_52', 'LDA_54', 'GGA', 'USPP_GGA', 'USPP_LDA' ]

def load_vasp_summary( filename ):
    """
    Reads a `vasp_summary.yaml` format YAML file and returns
    a dictionary of dictionaries. Each YAML document in the file
    corresponds to one sub-dictionary, with the corresponding
    top-level key given by the `title` value.

    Example:
        The file::
    
            ---
            title: foo
            data: foo_data
            ---
            title: bar
            data: bar_data

        is converted to the dictionary::

            { 'foo': { 'title': 'foo', 'data': 'foo_data' },
              'bar': { 'title': 'bar', 'data': 'bar_data' } }

    Args:
        filename (str): File path for the `vasp_summary.yaml` file.

    Returns:
        (dict(dict,dict,...)): A dictionary of separate YAML documents,
            each as dictionaries.a

    """
    with open( filename, 'r' ) as stream:
        docs = yaml.load_all( stream, Loader=yaml.SafeLoader )
        data = { d['title']: d for d in docs }
    return data

def potcar_spec( filename ):
    """
    Returns a dictionary specifying the pseudopotentials contained in a POTCAR file.

    Args:
        filename (Str): The name of the POTCAR file to process.

    Returns:
        (Dict): A dictionary of pseudopotential filename: dataset pairs, e.g.
                { 'Fe_pv': 'PBE_54', 'O', 'PBE_54' }
    """
    p_spec = {}
    with open( filename, 'r' ) as f:
        potcars = re.split('(End of Dataset\n)', f.read() )
    potcar_md5sums = [ md5sum( ''.join( pair ) ) for pair in zip( potcars[::2], potcars[1:-1:2] ) ]
    for this_md5sum in potcar_md5sums:
        for ps in potcar_sets:
            for p, p_md5sum in potcar_md5sum_data[ ps ].items():
                if this_md5sum == p_md5sum:
                    p_spec[ p ] = ps
    if len( p_spec ) != len( potcar_md5sums ):
        raise ValueError( 'One or more POTCARs did not have matching md5 hashes' )
    return p_spec
  
def find_vasp_calculations():
    """
    Returns a list of all subdirectories that contain either a vasprun.xml file
    or a compressed vasprun.xml.gz file.

    Args:
        None

    Returns:
        (List): list of all VASP calculation subdirectories.
    """
    dir_list = [ './' + re.sub( r'vasprun\.xml', '', path ) for path in glob.iglob( '**/vasprun.xml', recursive=True ) ]
    gz_dir_list = [ './' + re.sub( r'vasprun\.xml\.gz', '', path ) for path in glob.iglob( '**/vasprun.xml.gz', recursive=True ) ]
    return dir_list + gz_dir_list

class Summary:
    """
    TODO Document Summary class
    """

    supported_flags = { 'title': 'Title',
                        'description': 'Description',
                        'notes': 'Notes',
                        'type': 'Type',
                        'status': 'Status',
                        'stoichiometry': 'Stoichiometry',
                        'potcar': 'POTCAR',
                        'eatom': 'POTCAR EATOM values',
                        'plus_u': 'Dudarev +U parameters',
                        'energy': 'Energy',
                        'lreal': 'LREAL',
                        'k-points': 'k-points',
                        'functional': 'functional',
                        'encut': 'encut',
                        'ediffg': 'ediffg',
                        'ibrion': 'ibrion',
                        'converged': 'converged',
                        'md5': 'md5',
                        'directory': 'directory',
                        'vbm': 'Vasprun valence band maximum',
                        'cbm': 'Vasprun conduction band minimum',
                        'track': 'tracking for files',
                        'version': 'VASP executable version',
                        'nelect': 'NELECT' }

    def __init__( self, directory='.' ):
        self.directory = directory
        with cd( directory ):
            try:
                self.meta = VASPMeta.from_file( 'vaspmeta.yaml' )
            except FileNotFoundError as e:
                raise type(e)( str(e) + ' in {}'.format( directory )).with_traceback( sys.exc_info()[2] )
            self.parse_vasprun()
        self.print_methods = { 'title': self.print_title,
                               'description': self.print_description,
                               'notes': self.print_notes,
                               'type': self.print_type,
                               'status': self.print_status,
                               'stoichiometry': self.print_stoichiometry,
                               'potcar': self.print_potcar,
                               'eatom': self.print_eatom,
                               'energy': self.print_energy,
                               'k-points': self.print_kpoints,
                               'functional': self.print_functional,
                               'encut': self.print_encut,
                               'plus_u': self.print_plus_u,
                               'ediffg': self.print_ediffg,
                               'ibrion': self.print_ibrion,
                               'converged': self.print_converged,
                               'version': self.print_version,
                               'md5': self.print_vasprun_md5,
                               'directory': self.print_directory,
                               'lreal': self.print_lreal,
                               'vbm': self.print_vbm,
                               'cbm': self.print_cbm,
                               'track': self.print_file_tracking, 
                               'nelect': self.print_nelect }
        if not set( self.print_methods.keys() ) == set( self.supported_flags ):
            print( set( self.print_methods.keys() ) )
            print( '--------------' )
            print( set( self.supported_flags.keys() ) )
            raise( ValueError )

    def parse_vasprun( self ):
        """
        Read in `vasprun.xml` as a pymatgen Vasprun object.

        Args:
            None

        Returns:
            None

        None:
            If the vasprun.xml is not well formed this method will catch the ParseError
            and set self.vasprun = None.
        """            
        self.vasprun_filename = match_filename( 'vasprun.xml' )
        if not self.vasprun_filename:
            raise FileNotFoundError( 'Could not find vasprun.xml or vasprun.xml.gz file' )
        try:
            self.vasprun = Vasprun( self.vasprun_filename, parse_potcar_file=False )
        except ET.ParseError:
            self.vasprun = None
        except:
            raise

    @property
    def stoich( self ):
        return self.vasprun.final_structure.composition.get_el_amt_dict()

    @property
    def functional( self ):
        """
        String description of the calculation functional.
       
        Recognises:
            - PBE
            - PBEsol
            - PBE-based hybrids:
                - PBE0 (alpha=0.25, no screening)
                - HSE06 (alpha=0.25, mu=0.2)
                - generic hybrids (alpha=?, no screening)
                - generic screened hybrids (alpha=?, mu=?)
        
        Returns:
            (Str): String describing the calculation functional.

        """
        if self.potcars_are_pbe(): # PBE base funtional
            if 'LHFCALC' in self.vasprun.parameters:
                alpha = float( self.vasprun.parameters['AEXX'] )
            else:
                alpha = 0.0
            if 'HFSCREEN' in self.vasprun.parameters:
                mu = float( self.vasprun.parameters['HFSCREEN'] )
            else:
                mu = 0
            if alpha > 0:
                if mu > 0: # screened hybrid
                    if ( mu == 0.2 ) and ( alpha == 0.25 ):
                        f = 'HSE06'
                    else:
                        f = "screened hybrid. alpha={}, mu={}".format( alpha, mu )
                else: # unscreened hybrid
                    if alpha == 0.25:
                        f = 'PBE0'
                    else: 
                        f = "hybrid. alpha={}".format( alpha )
            else: # not hybrid. Plain PBE or some variant.
                pbe_list = { 'PS': 'PBEsol',
                             'PE': 'PBE',
                             '91': 'PW91',
                             'RP': 'rPBE',
                             'AM': 'AM05' }
                f = pbe_list[ self.vasprun.parameters['GGA'] ]
        else:
            f = 'not recognised'    
        return f

    def potcars_are_pbe( self ):
        return all( 'PBE' in s for s in self.vasprun.potcar_symbols )

    def output( self, to_print ):
        if not self.vasprun:
            to_print = [ 'title', 'type', 'status' ]
        print( "---" )
        for p in to_print:
            self.print_methods[ p ]()
        print( '', flush=True )
       
    def print_type( self ):
        if self.meta.type:
            print( "type: {}".format( self.meta.type ) )
 
    def print_title( self ):
        print( "title: {}".format( self.meta.title ) )
      
    def print_description( self ):
        print( "description: {}".format( self.meta.description.strip() ) )
 
    def print_notes( self ):
        if self.meta.notes:
            print( "notes: {}".format( self.meta.notes.strip() ) )
        else:
            print( "notes: ~" )
 
    def print_status( self ):
        print( "status: {}".format( self.meta.status ) )

    def print_lreal( self ):
        print( "lreal: {}".format( self.vasprun.parameters['LREAL'] ) )

    def print_stoichiometry( self ):
        print( "stoichiometry:" )
        for element in self.stoich:
            print( "    - {}: {}".format( element, int( self.stoich[ element ] ) ) )

    def print_potcar( self ):
        print( "potcar:" )
        for e, p in zip( self.stoich, self.vasprun.potcar_symbols ):
            print( "    - {}: {}".format( e, p ) )
        
    def print_energy( self ):
        # if this gets more options, it might be a good idea to set the
        # appropriate method using a dictionary?
        # or we could subclass Summary --> NEB_Summary ?
        if not self.meta.type:
            print( "energy: {}".format( self.vasprun.final_energy ) )    
        elif self.meta.type == 'neb':
            self.print_neb_energy()
        else:
            raise ValueError( "VASPMeta type not supported: {}".format( self.meta.type ) )

    def print_neb_energy( self ):
        image_00_energy = final_energy_from_outcar( '00/OUTCAR' )
        print( "reference energy: {} eV".format( image_00_energy ) )
        neb = NEBAnalysis.from_dir( '.' )
        print( "neb image energies:" )
        for i, e in enumerate( neb.energies ):
            print( "    - {:02d}: {:10.6f} eV".format( i, e ) )

    def print_version( self ):
        version_string = vasp_version_from_outcar( '{}/OUTCAR'.format( self.directory ) ).split()[0]
        print( "version: {}".format( version_string ) )

    def print_eatom( self ):
        # This is one way to try to uniquely identify the POTCARs used, because the
        # potcar_symbol (e.g. `Ti_pv 07Sep2000`) is not sufficient.
        print( "eatom:" )
        for e, eatom in zip( self.stoich, potcar_eatom_list_from_outcar( '{}/OUTCAR'.format( self.directory ) ) ):
            print( "    - {}: {} eV".format( e, eatom ) )
        
    def print_kpoints( self ):
        print( "k-points:" )
        print( "    scheme: {}".format( self.vasprun.kpoints.style ) )
        print( "    grid: {}".format( " ".join( str( k ) for k in self.vasprun.kpoints.kpts[0] ) ) )

    def print_functional( self ):
        print( "functional: {}".format( self.functional ) )

    def print_ibrion( self ):
        print( "ibrion: {}".format( self.vasprun.incar['IBRION'] ) )

    def print_ediffg( self ):
        print( "ediffg: {}".format( self.vasprun.incar['EDIFFG'] ) )

    def print_encut( self ):
        if 'ENCUT' in self.vasprun.incar:
            print( "encut: {}".format( self.vasprun.incar['ENCUT'] ) )
        elif 'ENMAX' in self.vasprun.incar:
            print( "encut: {}".format( self.vasprun.incar['ENMAX'] ) )

    def print_converged( self ):
        print( "converged: {}".format( self.vasprun.converged ) )

    def print_vasprun_md5( self ):
        print( "vasprun md5: {}".format( file_md5( "{}/{}".format( self.directory, self.vasprun_filename ) ) ) )

    def print_file_tracking( self ):
        if self.meta.track:
            print( "file tracking:" )
            for f, new_filename in self.meta.track.items():
                print( "    {}:".format( f ) )
                if not new_filename:
                    new_filename = f
                print( "        filename: {}".format( new_filename ) )
                filename = match_filename( self.directory + f )
                if filename: 
                    md5 = file_md5( filename )
                else:
                    md5 = 'null'
                print( "        md5: {}".format( md5 ) )
 
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

    def print_cbm( self ):
        print( 'cbm: {}'.format( self.vasprun.eigenvalue_band_properties[1] ) )

    def print_vbm( self ):
        print( 'vbm: {}'.format( self.vasprun.eigenvalue_band_properties[2] ) )

    def print_nelect( self ):
        print( 'nelect: {}'.format( self.vasprun.parameters['NELECT'] ) )
        
