import numpy as np

class AutoKPoints:
    """
    class for automatic k-point generation data in KPOINTS.
    """

    def __init__( self, title, subdivisions, grid_centering='G', shift=np.array( [ 0., 0., 0. ] ) ):
        """
        Initialise an AutoKPoints object

        Args:
            title (Str): The first line of the file, treated as a comment by VASP.
            grid_centering (Str, optional): Specify gamma-centered (G) or the original Monkhorst-Pack scheme (MP). Default is 'G'.
            subdivisions: (np.Array( Int, Int, Int )): Numbers of subdivisions along each reciprocal lattice vector.
            shift: (np.Array( Float, Float, Float ), optional): Optional shift of the mesh (s_1, s_2, s_3). Default is ( [ 0., 0., 0. ] ).

        Returns:
            None
        """
        accepted_grid_centerings = [ 'G', 'MP' ]
        if grid_centering not in accepted_grid_centerings:
            raise ValueError
        self.title = title
        self.grid_centering = grid_centering
        self.subdivisions = subdivisions
        self.shift = shift

