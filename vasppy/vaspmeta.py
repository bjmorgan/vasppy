import yaml
from vasppy.utils import file_md5

class VASPMeta:
    """
    VASPMeta class for storing additional VASP calculation metadata
    """

    def __init__( self, title, description, status, notes=None, type=None, md5=None ):
        """
        Initialise a VASPMeta object.

        Args:
            title (Str): The title string for this calculation
            description (Str): Long description
            status (Str): Current status of the calculation. 
                Expected strings are (to-run, incomplete, finished, dropped)
            notes (:obj:Str, optional): Any additional notes. Defaults to None.
            type  (:obj:Str, optional): Can be used to describe the calculation type. 
                Defaults to None.
            md5 (:obj:list(str), optional): An optional list of filenames to calculate md5 hashes
                files to calculate hashes for when summarising the calculation output.
                Defaults to None.

        Returns:
            None
        """
        self.title = title
        self.description = description
        self.notes = notes
        expected_status = [ 'to-run', 'incomplete', 'finished', 'dropped' ]
        if status not in expected_status:
            raise ValueError( status )
        self.status = status
        expected_types = [ 'single-point', 'neb' ]
        if type:
            if type not in expected_types:
                raise ValueError( type )
            self.type = type 
        else:
            self.type = None
        self.md5 = md5

    @classmethod
    def from_file( cls, filename ):
        """
        Create a VASPMeta object by reading a `vaspmeta.yaml` file

        Args:
            filename (Str): filename to read in.

        Returns:
            (vasppy.VASPMeta): the VASPMeta object
        """
        with open( filename, 'r' ) as stream:
            data = yaml.load( stream )
            notes = data.get( 'notes' )
            v_type = data.get( 'type' )
            file_md5 = data.get( 'file_md5' )
            xargs = {}
            if file_md5:
                if type( file_md5 ) is str:
                    file_md5 = [ file_md5 ]
                xargs['md5'] = file_md5 
            vaspmeta = VASPMeta( data['title'], 
                                 data['description'], 
                                 data['status'], 
                                 notes=notes, 
                                 type=v_type,
                                 **xargs )
        return vaspmeta

