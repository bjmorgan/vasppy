import yaml
from vasppy.utils import file_md5

class VASPMeta:
    """
    VASPMeta class for storing additional VASP calculation metadata
    """

    def __init__( self, title, description, status, notes=None, type=None, track=None ):
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
            track (:obj:dict(str: str), optional): An optional dict of pairs of filenames
                for files to track. For each key: value pair, the key is the current filename
                in the directory. The value is the new filename that list of filenames to calculate md5 hashes
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
        self.track = track

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
            data = yaml.load( stream, Loader=yaml.SafeLoader )
            notes = data.get( 'notes' )
            v_type = data.get( 'type' )
            track = data.get( 'track' )
            xargs = {}
            if track:
                if type( track ) is str:
                    track = [ track ]
                xargs['track'] = track 
            vaspmeta = VASPMeta( data['title'], 
                                 data['description'], 
                                 data['status'], 
                                 notes=notes, 
                                 type=v_type,
                                 **xargs )
        return vaspmeta

