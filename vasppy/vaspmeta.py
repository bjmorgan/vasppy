import yaml

class VASPMeta:
    """
    VASPMeta class for storing additional VASP calculation metadata
    """

    def __init__( self, title, description, status, notes=None ):
        """
        Initialise a VASPMeta object.

        Args:
            title (Str): The title string for this calculation
            description (Str): Long description
            status (Str): Current status of the calculation. 
                Expected strings are (to-run, incomplete, finished, dropped)
            notes (:obj:Str, optional): Any additional notes. Defaults to None

        Returns:
            None
        """
        self.title = title
        self.description = description
        self.notes = notes
        expected_status = [ 'to-run', 'incomplete', 'finished', 'dropped' ]
        if status not in expected_status:
            raise ValueError
        self.status = status

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
            if 'notes' in data:
                vaspmeta = VASPMeta( data['title'], data['description'], data['status'], notes=data['notes'] )
            else:
                vaspmeta = VASPMeta( data['title'], data['description'], data['status'] )
            return vaspmeta
                
        
